"""Versioned result reading, tabular views and transactional run bundles."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import tempfile
import uuid

from .jsonio import to_jsonable, write_json

FIELDS = ('agent', 'analysis', 'method', 'estimand', 'analysis_scope', 'n', 'n_groups',
          'n_observed', 'n_missing', 'n_positive', 'estimate_scope', 'batch',
          'estimate', 'lower', 'upper', 'log_e', 'interval_kind', 'status', 'reason',
          'unit', 'assumptions')
OWNED_FILES = frozenset({'clonal_share_result.json', 'input_qc.json', 'input_qc.md',
                         'report.md', 'report.html', 'results.csv', 'strata_results.csv', 'run_manifest.json'})


def read_result(path):
    """Read archived schema1 or current schema2 without changing its meaning."""
    result = json.loads(Path(path).read_text(encoding='utf-8'))
    if not isinstance(result, dict) or result.get('schema_version') not in ('1.0', '2.0'):
        raise ValueError('Unsupported result schema; expected 1.0 or 2.0')
    return result


def _finite(value):
    return value if isinstance(value, (int, float)) and math.isfinite(value) else None


def result_rows(record):
    """One row per agent and estimand, derived solely from the canonical record."""
    md = record.get('metadata_diagnostics') or {}
    rows = []

    def add(agent, analysis, target, r, point=None, low=None, high=None,
            interval_kind='95% interval', scope='observed and typed subset',
            n=None, unit='proportion', assumptions='Method-specific sampling assumptions apply.',
            estimate_scope=None):
        complete = r.get('complete', True)
        available = r.get('estimable', True) and r.get('status') not in ('empty', 'unavailable', 'no_retained_calls')
        if not complete:
            status = 'computation_incomplete'
        elif r.get('empty'):
            status = 'empty_confidence_set'
        elif not available:
            status = 'unavailable'
        elif _finite(low) == 0 and _finite(high) == 1:
            status = 'full_range'
        elif all(_finite(v) is None for v in (point, low, high, r.get('log_e'))):
            status = 'unavailable'
        else:
            status = 'computed'
        # Diagnostic estimates remain in the source record; unavailable values
        # are never exported as reportable estimates in the primary CSV.
        usable = complete and available and not r.get('empty')
        reason = r.get('reason') or r.get('failure_reason') or r.get('not_estimable_because') or ''
        if status == 'unavailable' and not reason:
            if r.get('n') == 0:
                reason = 'No readable outcome with a recorded lineage'
            elif r.get('n_groups') == 0:
                reason = 'No retained lineage labels'
            elif r.get('n_groups') == 1:
                reason = 'Only one retained lineage; no between-lineage comparison'
            elif r.get('n_splits') == 0:
                reason = 'No eligible training and evaluation split'
            elif r.get('prevalence') in (0., 1.):
                reason = 'Constant retained outcome; no variation to attribute'
            else:
                reason = 'No reportable result under the recorded method conditions'
        if not reason:
            reason = {
                'full_range': 'The completed interval or bounds span the entire admissible range [0, 1]',
                'computation_incomplete': 'The requested computation did not complete; no interval is reportable',
                'empty_confidence_set': 'No parameter value was retained by the completed confidence-set calculation',
            }.get(status, '')
        rows.append(dict(agent=str(agent), analysis=analysis, method=r.get('estimator_method', r.get('method', analysis)), estimand=target,
                         analysis_scope=scope, n=r.get('n', n), n_groups=r.get('n_groups'),
                         n_observed=r.get('observed_count'), n_missing=r.get('missing_count'),
                         n_positive=r.get('positive_count'), estimate_scope=estimate_scope or scope,
                         batch=r.get('batch'),
                         estimate=_finite(point) if usable else None,
                         lower=_finite(low) if usable else None,
                         upper=_finite(high) if usable else None,
                         log_e=_finite(r.get('log_e')) if usable else None,
                         interval_kind=interval_kind, status=status,
                         reason=reason,
                         unit=unit, assumptions=assumptions))

    for agent, r in (record.get('collection_bounds') or {}).items():
        add(agent, 'finite_collection_prevalence', 'prevalence in recorded collection', r,
            r.get('observed_prevalence'), r.get('lower_bound'), r.get('upper_bound'),
            interval_kind='identification bounds; not a confidence interval',
            scope='recorded finite collection', n=r.get('total_count'),
            estimate_scope='observed outcomes only',
            assumptions='Binary outcomes; recorded collection size known. Estimate is observed-only prevalence; bounds cover missing outcomes.')
    for agent, r in (md.get('clonal_share') or {}).items():
        add(agent, 'collection_membership', 'collection lineage-membership share', r,
            r.get('kappa_adj'), r.get('ci_low'), r.get('ci_high'))
    for agent, r in (md.get('realised_share') or {}).items():
        add(agent, 'realised_component', 'realised variance-component ratio', r,
            r.get('kappa'), r.get('ci_low'), r.get('ci_high'),
            assumptions='Gaussian residual model for interval interpretation; diagnostic checks do not verify this assumption.')
    for agent, r in (md.get('population_probit_icc') or {}).items():
        add(agent, 'population_probit', r.get('target', 'Gaussian population liability ICC'), r,
            r.get('rho_hat'), r.get('ci_low'), r.get('ci_high'),
            interval_kind=r.get('confidence_kind', '95% model interval'),
            scope='population model conditional on retained grouped data',
            assumptions='Gaussian independent lineage effects; binomial observations; ignorable selection and group sizes. Finite validation scope.')
    for agent, entry in ((md.get('censored_share') or {}).get('per_agent') or {}).items():
        r = entry.get('share') or {}
        add(agent, 'censored_mic', 'MIC latent population variance-component ratio', r,
            r.get('kappa'), r.get('ci_low'), r.get('ci_high'),
            n=entry.get('n'), unit='proportion',
            interval_kind='approximate nominal 95% F interval; not a coverage guarantee',
            scope='population model fitted to readable and typed MIC subset',
            assumptions='Empirical-Bayes moment estimate, not maximum likelihood; Gaussian variance-components; ignorable sampling and coarsening. Nominal coverage is not guaranteed. MIC units: ' + str(entry.get('mic_units') or 'not supplied'))
        c = entry.get('calibrated')
        if c is not None:
            add(agent, 'censored_mic_calibrated', 'MIC latent population variance-component ratio', c,
                c.get('estimate'), c.get('low'), c.get('high'),
                n=entry.get('n'), unit='proportion',
                interval_kind=('exact 95% generalized F interval' if c.get('method') == 'exact_generalized_F'
                               else '95% interval by null-wise parametric-bootstrap likelihood-ratio inversion'),
                scope='population model fitted to readable and typed MIC subset',
                assumptions='Maximum-likelihood estimate; Gaussian lineage effects and residuals; independent lineages; ignorable sampling and coarsening; panel as recorded (' + str(entry.get('panel_source')) + ')' + ('; fixed effects of ' + ', '.join(c.get('covariate_columns')) if c.get('covariate_columns') else '') + '. Coverage was checked by simulation within this model. MIC units: ' + str(entry.get('mic_units') or 'not supplied'))
    for agent, r in ((md.get('lineage_evidence') or {}).get('per_feature') or {}).items():
        add(agent, 'lineage_evidence', 'evidence against lineage-independent binary outcomes', dict(r, estimable=r.get('n_splits', 1) > 0),
            r.get('e_value'), interval_kind='e-value; not an effect size', unit='e-value',
            assumptions='Held-out Bernoulli likelihood under the stated null; training data separate from evaluation.')
    sequential = (md.get('lineage_evidence') or {}).get('sequential') or {}
    for agent, r in (sequential.get('per_feature') or {}).items():
        logs, values = r.get('log_e') or [], r.get('e_value') or []
        analyzed = r.get('n_analyzed_per_batch')
        s = dict(r, log_e=logs[-1] if logs else None,
                 n=sum(analyzed) if analyzed is not None else None,
                 batch=(sequential.get('batches') or [None])[-1])
        if 'n_scored_batches' in r and not r['n_scored_batches']:
            s.update(estimable=False, reason='No batch had eligible prior training and readable evaluation outcomes; evidence remains initialized')
        add(agent, 'sequential_lineage_evidence', 'final cumulative evidence for lineage dependence', s,
            values[-1] if values else None, unit='e-value', interval_kind='e-process final look; not an effect size',
            scope='readable typed outcomes in prespecified ordered batches',
            assumptions='Disjoint batches; common-probability conditionally independent Bernoulli null within each batch; past-only training. Panel stopping requires joint-filtration validity; repeated rejection unions are not controlled.')
    for agent, r in ((md.get('prevalence_decomposition') or {}).get('per_feature') or {}).items():
        for component in ('composition', 'within_lineage'):
            limits = r.get(component + '_ci95') or [None, None]
            sub = dict(r, estimable=r.get(component + '_estimable', False))
            add(agent, 'decomposition_' + component, component + ' component of prevalence difference', sub,
                r.get(component), limits[0], limits[1], n=(r.get('n_a', 0) + r.get('n_b', 0)),
                assumptions='Descriptive labelled/observed subset; bootstrap sampling assumptions; no causal interpretation or representativeness certification.')
    # Preserve numerical streams by leaving empty columns outside estimators,
    # while explicitly accounting for each requested but unavailable analysis.
    cfg = record.get('config') or {}
    requested = []
    if (cfg.get('attribution') or {}).get('enabled', True):
        requested += [('collection_membership', 'collection lineage-membership share'),
                      ('realised_component', 'realised variance-component ratio')]
    if (cfg.get('evidence') or {}).get('enabled', True):
        requested += [('lineage_evidence', 'evidence against lineage-independent binary outcomes')]
        if (cfg.get('dataset') or {}).get('batch_column'):
            requested += [('sequential_lineage_evidence', 'final cumulative evidence for lineage dependence')]
    if (cfg.get('population_probit') or {}).get('enabled', False):
        requested += [('population_probit', 'Gaussian population liability ICC')]
    if (cfg.get('surveillance') or {}).get('enabled', True) and (cfg.get('dataset') or {}).get('contrast_column'):
        requested += [('decomposition_composition', 'composition component of prevalence difference'),
                      ('decomposition_within_lineage', 'within-lineage component of prevalence difference')]
    existing = {(r['agent'], r['analysis']) for r in rows}
    for agent in (record.get('traits') or []) if record.get('schema_version') == '2.0' else []:
        for analysis, target in requested:
            if (str(agent), analysis) not in existing:
                add(agent, analysis, target, dict(n=0, n_groups=0, estimable=False,
                    reason='No readable call available for this requested analysis'), interval_kind='not computed',
                    unit={'lineage_evidence': 'e-value', 'sequential_lineage_evidence': 'e-value'}.get(analysis, 'proportion'))
    return to_jsonable(rows)


def status_guidance(row):
    """Explain an existing result status without changing it or its estimand."""
    status = row['status']
    reason = row.get('reason') or ''
    if status == 'computed':
        return reason or 'A reportable result was produced for the stated target and data subset.', 'Read the interval kind, assumptions and retained cohort before comparing results.'
    if status == 'full_range':
        action = ('Review missing-outcome counts; the recorded frame permits every prevalence from 0 to 1.'
                  if row['analysis'] == 'finite_collection_prevalence' else
                  'Report the full range as uninformative; inspect repeated-lineage support and outcome variation.')
    elif status == 'computation_incomplete':
        action = 'Inspect the recorded numerical or resource failure and method diagnostics before retrying in a new output directory.'
    elif status == 'empty_confidence_set':
        action = 'Inspect numerical diagnostics and model compatibility; do not replace the empty set with a zero estimate.'
    elif row['analysis'] in ('censored_mic', 'censored_mic_calibrated'):
        action = 'Check MIC parsing, panel wells, censoring and typed-isolate support in input QC and MIC diagnostics.'
    elif row['analysis'] == 'population_probit':
        action = 'Check typed outcomes, repeated groups and the selected method domain in population-model diagnostics.'
    elif row['analysis'].startswith('decomposition_'):
        action = 'Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope.'
    elif row['analysis'] == 'finite_collection_prevalence':
        action = 'Check the recorded collection frame and readable binary outcomes in input QC.'
    else:
        action = 'Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics.'
    return reason, action


def analysis_summary(record):
    groups = {}
    for row in result_rows(record):
        block = groups.setdefault(row['analysis'], {'n_agents': 0, 'n_computed': 0,
                                                   'n_full_range': 0, 'n_unavailable': 0})
        block['n_agents'] += 1
        if row['status'] == 'computed':
            block['n_computed'] += 1
        elif row['status'] == 'full_range':
            block['n_full_range'] += 1
        else:
            block['n_unavailable'] += 1
    return groups


def validate_destination(destination, *, overwrite=False):
    dest = Path(destination).expanduser().resolve()
    if dest.exists() and not dest.is_dir():
        raise FileExistsError(f'Output destination is not a directory: {dest}')
    if dest.exists() and any(dest.iterdir()) and not overwrite:
        raise FileExistsError(f'Output directory is not empty: {dest}. Choose a new directory or use --overwrite.')
    return dest


def _umask():
    mask = os.umask(0)
    os.umask(mask)
    return mask


def _publish(destination, writer, *, overwrite=False):
    dest = validate_destination(destination, overwrite=overwrite)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f'.{dest.name}.staging-', dir=dest.parent) as tmp:
        staging = Path(tmp).resolve()
        writer(staging)
        # A temporary directory is private to its owner; the published
        # directory gets the permissions an ordinary mkdir would give it.
        staging.chmod(0o777 & ~_umask())
        # Rendering may take time; a second check protects a result or user
        # file created at this destination since the initial preflight.
        validate_destination(dest, overwrite=overwrite)
        # Preserve user files in the new view. Previous complete runs remain
        # recoverable as sibling backups; no recursive user-directory deletion.
        if dest.exists():
            for p in dest.iterdir():
                if p.name not in OWNED_FILES:
                    if p.is_symlink():
                        raise ValueError(f'Cannot copy output-directory symlink: {p}')
                    if p.is_dir():
                        shutil.copytree(p, staging / p.name, symlinks=True)
                    else:
                        shutil.copy2(p, staging / p.name)
        backup = None
        if staging.parent != dest.parent:
            raise RuntimeError('Output staging escaped destination parent')
        if dest.exists():
            backup = dest.with_name(f'.{dest.name}.previous-{uuid.uuid4().hex}')
            dest.rename(backup)
        try:
            staging.rename(dest)
        except BaseException:
            if backup is not None and not dest.exists():
                backup.rename(dest)
            raise
    return dest


def _manifest(directory, kind):
    files = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
             for p in sorted(directory.iterdir()) if p.is_file()}
    write_json({'status': 'complete', 'kind': kind, 'files': files}, directory / 'run_manifest.json')


def publish_input_check(qc, destination, *, overwrite=False):
    from .qc import render_markdown
    def writer(directory):
        write_json(qc, directory / 'input_qc.json')
        (directory / 'input_qc.md').write_text(render_markdown(qc) if qc else '# Input check\n\nNo input QC was stored with this record.\n', encoding='utf-8')
        _manifest(directory, 'input_check_only')
    return _publish(destination, writer, overwrite=overwrite)


def publish_results(record, destination, *, overwrite=False):
    from .cli import _summary
    from .qc import render_markdown
    from .report import render_report
    from .report_html import render_html_report
    def writer(directory):
        write_json(record, directory / 'clonal_share_result.json')
        # Both reports and CSV see exactly the strict JSON representation.
        saved = read_result(directory / 'clonal_share_result.json')
        summary = _summary(saved)
        qc = saved.get('input_qc') or {}
        write_json(qc, directory / 'input_qc.json')
        (directory / 'input_qc.md').write_text(render_markdown(qc) if qc else '# Input check\n\nNo input QC was stored with this record.\n', encoding='utf-8')
        record_bytes = (directory / 'clonal_share_result.json').read_bytes()
        (directory / 'report.md').write_text(render_report(saved, summary, input_qc=qc, record_bytes=record_bytes), encoding='utf-8')
        (directory / 'report.html').write_text(render_html_report(saved, summary, input_qc=qc, record_bytes=record_bytes), encoding='utf-8')
        with (directory / 'results.csv').open('w', newline='', encoding='utf-8-sig') as stream:
            writer_csv = csv.DictWriter(stream, fieldnames=FIELDS)
            writer_csv.writeheader()
            for row in result_rows(saved):
                writer_csv.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v for k, v in row.items()})
        strata = saved.get('strata')
        if strata:
            with (directory / 'strata_results.csv').open('w', newline='', encoding='utf-8-sig') as stream:
                writer_csv = csv.DictWriter(stream, fieldnames=('stratum',) + FIELDS)
                writer_csv.writeheader()
                for level, block in strata['levels'].items():
                    for row in block['results']:
                        writer_csv.writerow({'stratum': level, **{k: json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v for k, v in row.items()}})
        _manifest(directory, 'analysis')
    return _publish(destination, writer, overwrite=overwrite)
