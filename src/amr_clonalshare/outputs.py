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
from .verdict import agent_verdict

FIELDS = ('agent', 'analysis', 'method', 'estimand', 'analysis_scope', 'n', 'n_groups',
          'n_observed', 'n_missing', 'n_positive', 'estimate_scope', 'batch',
          'estimate', 'lower', 'upper', 'log_e', 'p_value', 'q_value', 'selected',
          'interval_kind', 'status', 'reason', 'unit', 'assumptions', 'conclusion')
OWNED_FILES = frozenset({'clonal_share_result.json', 'input_qc.json', 'input_qc.md',
                         'report.md', 'report.html', 'results.csv', 'strata_results.csv', 'run_manifest.json',
                         'run.log', 'run_failed.log'})


def read_result(path):
    """Read a canonical result record without changing its meaning."""
    result = json.loads(Path(path).read_text(encoding='utf-8'))
    if not isinstance(result, dict) or result.get('schema_version') != '1.0':
        raise ValueError('Unsupported result schema; expected 1.0')
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
            estimate_scope=None, conclusion='', p=None, q=None, selected=None):
        available = r.get('estimable', True) and r.get('status') not in ('empty', 'unavailable', 'no_retained_calls')
        if not available:
            status = 'unavailable'
        elif _finite(low) == 0 and _finite(high) == 1:
            status = 'full_range'
        elif all(_finite(v) is None for v in (point, low, high, r.get('log_e'))):
            status = 'unavailable'
        else:
            status = 'computed'
        usable = available
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
        if not reason and status == 'full_range':
            reason = 'The completed interval or bounds span the entire admissible range [0, 1]'
        rows.append(dict(agent=str(agent), analysis=analysis, method=r.get('estimator_method', r.get('method', analysis)), estimand=target,
                         analysis_scope=scope, n=r.get('n', n), n_groups=r.get('n_groups'),
                         n_observed=r.get('observed_count'), n_missing=r.get('missing_count'),
                         n_positive=r.get('positive_count'), estimate_scope=estimate_scope or scope,
                         batch=r.get('batch'),
                         estimate=_finite(point) if usable else None,
                         lower=_finite(low) if usable else None,
                         upper=_finite(high) if usable else None,
                         log_e=_finite(r.get('log_e')) if usable else None,
                         p_value=_finite(p) if usable else None,
                         q_value=_finite(q) if usable else None,
                         selected=('yes' if selected else 'no') if usable and selected is not None else None,
                         interval_kind=interval_kind, status=status,
                         reason=reason,
                         unit=unit, assumptions=assumptions, conclusion=conclusion))

    for agent, r in (record.get('collection_bounds') or {}).items():
        add(agent, 'finite_collection_prevalence', 'prevalence in recorded collection', r,
            r.get('observed_prevalence'), r.get('lower_bound'), r.get('upper_bound'),
            interval_kind='identification bounds; not a confidence interval',
            scope='recorded finite collection', n=r.get('total_count'),
            estimate_scope='observed outcomes only',
            assumptions='Binary outcomes; recorded collection size known. Estimate is observed-only prevalence; bounds cover missing outcomes.')
    latent_assumptions = ('Every reading holds the latent MIC it stands for (a call: the side of '
                          'its cut-off); the latent MIC is continuous within a stratum; lineages '
                          'of one isolate set aside. No distributional form for the latent values. '
                          'The bounds are sharp where the distinct readings of each stratum are disjoint; overlapping readings are merged into one cell and the bounds then hold but may not be attained.')
    # The sampling unit, where declared, is read by the permutation test alone:
    # the estimates, intervals and limits treat the isolates of a lineage as
    # independent whatever they share, and say so.
    units_note = ((' Sampling units (' + str((record.get('metadata_diagnostics') or {}).get('sampling_units', {}).get('column'))
                   + ') are declared; only the p-value reads them, exchanging lineage labels within a unit.')
                  if (record.get('metadata_diagnostics') or {}).get('sampling_units') else '')
    independent = 'Isolates of a lineage drawn independently from its distribution, whatever their sampling units share'
    limit_assumptions = (latent_assumptions + ' The limit also assumes readings independent between '
                         'isolates, whatever their sampling units share.' + units_note)

    def latent_rows(agent, r, prefix, n=None, scope='observed and typed subset'):
        if 'latent_order_lower' not in r:
            return
        exact = bool(r.get('latent_order_upper_exact'))
        kind = ('bounds the readings allow; not a confidence interval; lower end established by the readings'
                + ('; upper end sharp' if exact else '; upper end a certified bound, not sharp'))
        found = r.get('latent_order_upper')
        note = ('' if exact or not isinstance(found, (int, float)) or not math.isfinite(found) else
                f'The largest share attained by an arrangement of the lineages inside the readings is {found:.3f}')
        add(agent, prefix + '_latent_bounds', 'lineage share of the latent ordering', dict(r, reason=note),
            None, r.get('latent_order_lower'), r.get('latent_order_upper_bound'),
            interval_kind=kind, n=n, scope=scope, assumptions=latent_assumptions)
        add(agent, prefix + '_latent_lower_limit', 'lineage share of the latent ordering', dict(r, reason=''),
            None, r.get('latent_order_lower_limit'), None,
            interval_kind='one-sided 95% lower confidence limit for the lower bound, and so for the share',
            n=n, scope=scope, assumptions=limit_assumptions)

    call_strata = md.get('call_strata') or {}
    call_within = (('each level of ' + str(call_strata.get('column'))) if call_strata.get('column')
                   else 'the collection')
    call_selection = md.get('lineage_selection') or {}
    call_q = call_selection.get('q_values') or {}
    call_selected = set(call_selection.get('rejected_features') or [])
    alpha = float(call_selection.get('alpha') or 0.05)
    for agent, r in (md.get('clonal_share') or {}).items():
        add(agent, 'collection_membership', 'lineage share of the call among the represented lineages', r,
            r.get('kappa_adj'), r.get('observed_low'), r.get('observed_high'),
            interval_kind='95% interval for the represented lineages (isolates drawn again, lineages fixed)',
            assumptions=independent + '; the p-value assumes lineage labels exchangeable within ' + call_within
                        + ' under the null.' + units_note,
            conclusion=agent_verdict(r, alpha=alpha, q=call_q.get(agent))['text'],
            p=r.get('p_value'), q=call_q.get(agent), selected=(agent in call_selected) if call_selection else None)
        latent_rows(agent, r, 'call')
    order_selection = (md.get('censored_share') or {}).get('order_selection') or {}
    order_q = order_selection.get('q_values') or {}
    order_selected = set(order_selection.get('rejected_features') or [])
    order_alpha = float(order_selection.get('alpha') or 0.05)
    for agent, entry in ((md.get('censored_share') or {}).get('per_agent') or {}).items():
        r = entry.get('order') or dict(n=entry.get('n'), estimable=False, reason='not computed')
        add(agent, 'mic_order', 'lineage share of the recorded within-stratum MIC ordering among the represented lineages', r,
            r.get('kappa_adj'), r.get('observed_low'), r.get('observed_high'),
            interval_kind='95% interval for the represented lineages (isolates drawn again, lineages fixed)',
            n=r.get('n_readings', entry.get('n')), unit='proportion',
            scope='readable and typed MIC subset',
            assumptions='No distributional form for the readings; readings scored within ' + (('each level of ' + str(r['stratified_by'])) if r.get('stratified_by') else 'the collection') + '; isolates of a lineage drawn independently from its distribution, whatever their sampling units share.' + units_note + ' MIC units: ' + str(entry.get('mic_units') or 'not supplied'),
            conclusion=agent_verdict(r, alpha=order_alpha, q=order_q.get(agent))['text'] if entry.get('order') else '',
            p=r.get('p_value'), q=order_q.get(agent), selected=(agent in order_selected) if order_selection else None)
        latent_rows(agent, r, 'mic', n=r.get('n_readings', entry.get('n')), scope='readable and typed MIC subset')
    e_bh = (md.get('lineage_evidence') or {}).get('e_bh') or {}
    e_selected = set(e_bh.get('rejected_features') or [])
    for agent, r in ((md.get('lineage_evidence') or {}).get('per_feature') or {}).items():
        add(agent, 'lineage_evidence', 'evidence against lineage-independent binary outcomes', dict(r, estimable=r.get('n_splits', 1) > 0),
            r.get('e_value'), interval_kind='e-value; not an effect size', unit='e-value',
            selected=(agent in e_selected) if e_bh else None,
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
                p=r.get(component + '_p'), q=r.get(component + '_q'),
                selected=r.get(component + '_discovery') if (component + '_discovery') in r else None,
                assumptions='Descriptive labelled/observed subset; bootstrap sampling assumptions; no causal interpretation or representativeness certification.')
    # Preserve numerical streams by leaving empty columns outside estimators,
    # while explicitly accounting for each requested but unavailable analysis.
    cfg = record.get('config') or {}
    requested = []
    if (cfg.get('attribution') or {}).get('enabled', True):
        requested += [('collection_membership', 'lineage share of the call among the represented lineages')]
    if (cfg.get('evidence') or {}).get('enabled', True):
        requested += [('lineage_evidence', 'evidence against lineage-independent binary outcomes')]
        if (cfg.get('dataset') or {}).get('batch_column'):
            requested += [('sequential_lineage_evidence', 'final cumulative evidence for lineage dependence')]
    if (cfg.get('surveillance') or {}).get('enabled', True) and (cfg.get('dataset') or {}).get('contrast_column'):
        requested += [('decomposition_composition', 'composition component of prevalence difference'),
                      ('decomposition_within_lineage', 'within-lineage component of prevalence difference')]
    existing = {(r['agent'], r['analysis']) for r in rows}
    for agent in (record.get('traits') or []) if record.get('schema_version') == '1.0' else []:
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
        return reason or 'A reportable result was produced for the stated target and data subset.', 'Read the interval kind, assumptions and retained collection before comparing results.'
    if status == 'full_range':
        action = ('Review missing-outcome counts; the recorded frame permits every prevalence from 0 to 1.'
                  if row['analysis'] == 'finite_collection_prevalence' else
                  'Report the full range as uninformative; inspect repeated-lineage support and outcome variation.')
    elif row['analysis'].startswith('mic_'):
        action = 'Check MIC parsing, panel wells, censoring and typed-isolate support in input QC and MIC diagnostics.'
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


def _cell(value):
    """A CSV cell a spreadsheet will not run as a formula.

    A text cell that opens with =, +, -, @, a tab or a carriage return is
    prefixed with an apostrophe unless it is a number; every other value is
    written as it is.
    """
    if isinstance(value, (list, dict)):
        value = json.dumps(value, ensure_ascii=False)
    if isinstance(value, str) and value[:1] in ('=', '+', '-', '@', '\t', '\r'):
        try:
            float(value)
        except ValueError:
            return "'" + value
    return value


def validate_destination(destination, *, overwrite=False):
    """Refuse a destination before any computation, not after it.

    ``overwrite`` replaces a results directory; it is refused for a directory
    holding subdirectories or the working directory, which is a project folder
    rather than a results folder, and the parent must be writable now.
    """
    dest = Path(destination).expanduser().resolve()
    if dest.exists() and not dest.is_dir():
        raise FileExistsError(f'Output destination is not a directory: {dest}')
    if dest.exists() and any(dest.iterdir()):
        if not overwrite:
            raise FileExistsError(f'Output directory is not empty: {dest}. Choose a new directory or use --overwrite.')
        cwd = Path.cwd().resolve()
        if cwd == dest or dest in cwd.parents:
            raise FileExistsError(f'--overwrite would replace the working directory {dest}; choose a results directory inside it.')
        for p in dest.iterdir():
            if p.is_symlink():
                raise FileExistsError(f'Output directory holds a symbolic link, which is not copied: {p}')
            if p.is_dir():
                raise FileExistsError(f'Output directory holds the subdirectory {p.name}, so it is not a results '
                                      f'directory; --overwrite replaces only results directories.')
    parent = dest.parent
    while not parent.exists():
        parent = parent.parent
    try:
        Path(tempfile.mkdtemp(prefix=f'.{dest.name}.check-', dir=parent)).rmdir()
    except OSError as exc:
        raise PermissionError(f'Cannot write in {parent}: {exc.strerror or exc}') from exc
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
        # Preserve user files in the new view; the directory holds files only.
        if dest.exists():
            for p in dest.iterdir():
                if p.name not in OWNED_FILES:
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
    if backup is not None:
        # The replaced directory held files only, each now copied or replaced.
        for p in backup.iterdir():
            p.unlink()
        backup.rmdir()
    return dest


def write_failure_log(destination, lines, error: str) -> None:
    """The log of a run that did not complete, with the error at its end,
    written as the only file of the destination so that a reader can send
    it when asking for help; a later run replaces it. A destination that
    already holds a finished run keeps it as it was, and the log goes beside
    it as run_failed.log."""
    dest = Path(destination)
    name = 'run_failed.log' if (dest / 'run_manifest.json').is_file() else 'run.log'
    try:
        dest.mkdir(parents=True, exist_ok=True)
        (dest / name).write_text('\n'.join([*lines, f'error: {error}']) + '\n', encoding='utf-8')
    except OSError:
        pass


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


def publish_results(record, destination, *, overwrite=False, log=None):
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
                writer_csv.writerow({k: _cell(v) for k, v in row.items()})
        strata = saved.get('strata')
        if strata:
            with (directory / 'strata_results.csv').open('w', newline='', encoding='utf-8-sig') as stream:
                writer_csv = csv.DictWriter(stream, fieldnames=('stratum',) + FIELDS)
                writer_csv.writeheader()
                for level, block in strata['levels'].items():
                    for row in block['results']:
                        writer_csv.writerow({'stratum': _cell(level), **{k: _cell(v) for k, v in row.items()}})
        if log:
            (directory / 'run.log').write_text('\n'.join(log) + '\n', encoding='utf-8')
        _manifest(directory, 'analysis')
    return _publish(destination, writer, overwrite=overwrite)
