import json

import pytest

from amr_clonalshare import cli


def mic_record():
    return {'schema_version': '1.0', 'seed': 1, 'config': {}, 'n_isolates': 8,
            'n_traits': 0, 'traits': [], 'input_qc': {}, 'metadata_diagnostics': {
                'lineage_column': 'lineage', 'censored_share': {'per_agent': {
                    'A': {'n': 8, 'order': {'estimable': False, 'reason': 'limited panel'}}}}}}


def test_canonical_csv_preserves_unavailable_mic_and_reads_back(tmp_path):
    from amr_clonalshare.outputs import result_rows, read_result
    r = mic_record()
    rows = result_rows(r)
    assert len(rows) == 1
    assert rows[0]['agent'] == 'A'
    assert rows[0]['analysis'] == 'mic_order'
    assert rows[0]['status'] == 'unavailable'
    assert rows[0]['reason'] == 'limited panel'
    assert rows[0]['estimate'] is None
    p = tmp_path / 'record.json'
    p.write_text(json.dumps(r))
    assert read_result(p)['schema_version'] == '1.0'
    assert result_rows(read_result(p)) == rows
    r['schema_version'] = '0.1'
    p.write_text(json.dumps(r))
    with pytest.raises(ValueError, match='expected 1.0'):
        read_result(p)


@pytest.mark.parametrize('block, profile, expected_status, expected_reason, expected_action', [
    ('clonal_share', {'kappa_adj': .5, 'observed_low': 0., 'observed_high': 1.}, 'full_range',
     'entire admissible range', 'uninformative'),
    ('clonal_share', {'estimable': False, 'reason': 'no repeated groups'}, 'unavailable',
     'no repeated groups', 'recorded labels'),
])
def test_status_reasons_and_next_checks_match_both_reports(block, profile, expected_status,
                                                           expected_reason, expected_action):
    from amr_clonalshare.outputs import result_rows
    from amr_clonalshare.report import render_report
    from amr_clonalshare.report_html import render_html_report
    record = mic_record()
    record['metadata_diagnostics'] = {block: {'A': profile}}
    row = result_rows(record)[0]
    assert row['status'] == expected_status
    assert expected_reason in row['reason']
    if expected_status != 'full_range':
        assert row['estimate'] is None and row['lower'] is None and row['upper'] is None
    for render in (render_report, render_html_report):
        text = render(record, cli._summary(record))
        assert expected_status in text
        assert expected_reason in text
        assert expected_action in text


def test_bundle_failure_leaves_previous_directory_untouched(tmp_path, monkeypatch):
    from amr_clonalshare.outputs import publish_results
    from amr_clonalshare import report_html
    dest = tmp_path / 'out'
    dest.mkdir()
    marker = dest / 'report.md'
    marker.write_text('previous complete run')
    with pytest.raises(FileExistsError):
        publish_results(mic_record(), dest)
    def fail(*args, **kwargs):
        raise RuntimeError('deliberate rendering failure')
    monkeypatch.setattr(report_html, 'render_html_report', fail)
    with pytest.raises(RuntimeError, match='deliberate'):
        publish_results(mic_record(), dest, overwrite=True)
    assert marker.read_text() == 'previous complete run'
    assert sorted(p.name for p in dest.iterdir()) == ['report.md']


def test_bundle_records_own_hashes_and_protects_unrelated_files(tmp_path):
    import hashlib
    from amr_clonalshare.outputs import publish_results
    dest = tmp_path / 'out'
    publish_results(mic_record(), dest)
    manifest = json.loads((dest / 'run_manifest.json').read_text())
    assert manifest['status'] == 'complete'
    assert 'results.csv' in manifest['files']
    for name, sha in manifest['files'].items():
        assert hashlib.sha256((dest / name).read_bytes()).hexdigest() == sha
    (dest / 'my_notes.txt').write_text('keep')
    publish_results(mic_record(), dest, overwrite=True)
    assert (dest / 'my_notes.txt').read_text() == 'keep'


def test_mic_only_summary_and_report_do_not_claim_no_antimicrobial():
    from amr_clonalshare.report import render_report
    r = mic_record()
    summary = cli._summary(r)
    assert summary['analyses']['mic_order']['n_agents'] == 1
    rendered = render_report(r, summary).lower()
    assert 'no antimicrobial read' not in rendered
    assert 'no antimicrobial carried a readable call' not in rendered
    assert 'recorded mic readings' in rendered


def test_cli_protects_existing_outputs_before_loading_inputs(tmp_path, capsys):
    dest = tmp_path / 'out'
    dest.mkdir()
    marker = dest / 'report.md'
    marker.write_text('old')
    # A real config is not needed to prove no output mutation on bad input.
    assert cli.main(['--config', str(tmp_path / 'absent.yaml'),
                     '--results-dir', str(dest)]) == 2
    assert marker.read_text() == 'old'


def test_run_record_carries_bounds_and_canonical_rows(share_cfg):
    from amr_clonalshare import run
    _, cfg = share_cfg
    r = run(cfg, seed=7)
    assert r['schema_version'] == '1.0'
    assert set(r['collection_bounds']) == set(r['traits'])
    assert r['results']
    for b in r['collection_bounds'].values():
        assert b['lower_bound'] == b['upper_bound']


def test_new_data_arriving_during_render_is_not_overwritten(tmp_path, monkeypatch):
    from amr_clonalshare.outputs import publish_results
    from amr_clonalshare import report_html
    dest = tmp_path / 'out'
    render = report_html.render_html_report
    def concurrent_writer(*args, **kwargs):
        dest.mkdir()
        (dest / 'report.md').write_text('concurrent user result')
        return render(*args, **kwargs)
    monkeypatch.setattr(report_html, 'render_html_report', concurrent_writer)
    with pytest.raises(FileExistsError):
        publish_results(mic_record(), dest)
    assert (dest / 'report.md').read_text() == 'concurrent user result'


def test_canonical_rows_keep_denominators_scopes_and_log_evidence():
    from amr_clonalshare.outputs import result_rows
    r = mic_record()
    r['collection_bounds'] = {'B': dict(total_count=10, observed_count=6,
        missing_count=4, positive_count=3, observed_prevalence=.5, lower_bound=.3, upper_bound=.7)}
    r['metadata_diagnostics']['lineage_evidence'] = {'per_feature': {'B': dict(
        n=1000, n_groups=50, n_splits=5, e_value=None, log_e=900.)}}
    rows = {x['analysis']: x for x in result_rows(r)}
    b = rows['finite_collection_prevalence']
    assert (b['n'], b['n_observed'], b['n_missing'], b['n_positive']) == (10, 6, 4, 3)
    assert b['estimate_scope'] == 'observed outcomes only'
    e = rows['lineage_evidence']
    assert e['log_e'] == 900. and e['status'] == 'computed'
    assert e['estimate'] is None
    assert all(x['method'] for x in rows.values())


@pytest.mark.parametrize('calls,labels', [('R', 'A'), ('R\nS\nR', 'A\nA\nB'),
                                        ('\n\n', 'A\nB\nB'), ('R\nS', 'NULL\nNA')])
def test_partial_inputs_publish_usable_inventory(tmp_path, calls, labels):
    from amr_clonalshare.config import from_dict
    from amr_clonalshare import run
    outcomes, lineages = calls.split('\n'), labels.split('\n')
    (tmp_path / 'calls.csv').write_text('id,agent,call\n' + ''.join(
        f'{i},A,{v}\n' for i, v in enumerate(outcomes)))
    (tmp_path / 'meta.csv').write_text('id,lineage\n' + ''.join(
        f'{i},{v}\n' for i, v in enumerate(lineages)))
    cfg = from_dict(dict(dataset=dict(name='small', data_dir=str(tmp_path),
        metadata='meta.csv', strain_id_column='id', lineage_column='lineage',
        phenotype='calls.csv', phenotype_id_column='id', phenotype_antibiotic_column='agent',
        phenotype_call_column='call', phenotype_kind='clinical_sir'),
        attribution=dict(n_boot=0, n_perm=2, repeats=1), evidence=dict(repeats=1),
        surveillance=dict(n_boot=0)))
    r = run(cfg, results_dir=tmp_path / 'out', seed=7)
    assert r['n_isolates'] == len(outcomes)
    assert r['collection_bounds']['A']['total_count'] == len(outcomes)
    analyses = {x['analysis']: x for x in r['results']}
    assert {'collection_membership', 'lineage_evidence'} <= analyses.keys()
    if not any(outcomes):
        assert analyses['collection_membership']['status'] == 'unavailable'
    text = (tmp_path / 'out' / 'report.md').read_text(encoding='utf-8')
    assert 'R (resistant)' in text


def test_mic_export_leads_with_the_ordering_share():
    from amr_clonalshare.outputs import result_rows
    r = mic_record()
    agent = r['metadata_diagnostics']['censored_share']['per_agent']['A']
    agent['order'] = dict(estimable=True, kappa_adj=.3, observed_low=.15,
                          observed_high=.45, n_readings=8, stratified_by='laboratory',
                          latent_order_lower=.05, latent_order_upper=.7, latent_order_upper_bound=.8,
                          latent_order_upper_exact=False, latent_order_lower_limit=.02)
    rows = {row['analysis']: row for row in result_rows(r)}
    assert set(rows) == {'mic_order', 'mic_latent_bounds', 'mic_latent_lower_limit'}
    order = rows['mic_order']
    assert (order['estimate'], order['lower'], order['upper']) == (.3, .15, .45)
    bounds, limit = rows['mic_latent_bounds'], rows['mic_latent_lower_limit']
    # the upper end reported is the certified bound; the largest share found is named beside it
    assert (bounds['estimate'], bounds['lower'], bounds['upper']) == (None, .05, .8)
    assert 'not a confidence interval' in bounds['interval_kind']
    assert 'certified bound, not sharp' in bounds['interval_kind'] and '0.700' in bounds['reason']
    assert (limit['lower'], limit['upper']) == (.02, None) and 'lower confidence limit' in limit['interval_kind']
    assert 'each level of laboratory' in order['assumptions']
    agent['order'].update(latent_order_upper=.8, latent_order_upper_exact=True)
    bounds = {row['analysis']: row for row in result_rows(r)}['mic_latent_bounds']
    assert 'upper end sharp' in bounds['interval_kind'] and bounds['reason'] == ''


def test_overwrite_refuses_a_project_folder_and_the_working_directory(tmp_path, monkeypatch):
    from amr_clonalshare.outputs import validate_destination
    project = tmp_path / 'project'
    (project / 'data').mkdir(parents=True)
    (project / 'config.yaml').write_text('x')
    with pytest.raises(FileExistsError, match='subdirectory data'):
        validate_destination(project, overwrite=True)
    flat = tmp_path / 'flat'
    flat.mkdir()
    (flat / 'calls.csv').write_text('x')
    monkeypatch.chdir(flat)
    with pytest.raises(FileExistsError, match='working directory'):
        validate_destination(flat, overwrite=True)


def test_overwrite_leaves_no_hidden_backup(tmp_path):
    from amr_clonalshare.outputs import publish_results
    dest = tmp_path / 'out'
    publish_results(mic_record(), dest)
    publish_results(mic_record(), dest, overwrite=True)
    assert [p.name for p in tmp_path.iterdir()] == ['out']


def test_an_unwritable_parent_is_refused_before_the_run(tmp_path):
    import os
    from amr_clonalshare.outputs import validate_destination
    locked = tmp_path / 'locked'
    locked.mkdir()
    locked.chmod(0o500)
    try:
        if os.access(locked, os.W_OK):
            pytest.skip('running with privileges that ignore directory permissions')
        with pytest.raises(PermissionError):
            validate_destination(locked / 'out')
    finally:
        locked.chmod(0o700)


def test_csv_cells_that_open_like_a_formula_are_not_run():
    from amr_clonalshare.outputs import _cell
    assert _cell('=HYPERLINK("http://x")') == '\'=HYPERLINK("http://x")'
    assert _cell('@SUM(A1)') == "'@SUM(A1)"
    assert _cell('-0.25') == '-0.25'
    assert _cell(-0.25) == -0.25
    assert _cell('tetracycline') == 'tetracycline'
