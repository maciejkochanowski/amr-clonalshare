"""Reporting contracts for the numerical implementation."""
from pathlib import Path
from amr_clonalshare.outputs import result_rows


def test_mic_primary_export_labels_the_approximate_interval():
    record = {'metadata_diagnostics': {'censored_share': {'per_agent': {
        'a': {'share': {'kappa': .5, 'ci_low': .2, 'ci_high': .8,
                        'n': 40, 'n_groups': 4, 'estimable': True,
                        'estimator_method': 'empirical_Bayes_moment_iteration'}}}}}}
    row = next(r for r in result_rows(record) if r['analysis'] == 'censored_mic')
    assert row['method'] == 'empirical_Bayes_moment_iteration'
    assert 'approximate' in row['interval_kind'].lower()
    assert 'not guaranteed' in row['assumptions'].lower()


def test_ci_output_contract_uses_current_schema():
    text = (Path(__file__).resolve().parents[1]/'.github/workflows/ci.yml').read_text()
    assert "assert d['schema_version'] == '2.0'" in text


def test_report_renders_a_sequential_product_of_zero(tmp_path, monkeypatch):
    from amr_clonalshare import core, evalues
    from conftest import planted_cohort, write_config
    original = evalues.sequential_e_process

    class Zeroed:
        def __init__(self, result):
            self.record = result.as_dict()
            self.record['log_e'] = list(self.record['log_e'])
            self.record['log_e'][-1] = float('-inf')

        def as_dict(self):
            return self.record

    monkeypatch.setattr(evalues, 'sequential_e_process',
                        lambda *a, **k: Zeroed(original(*a, **k)))
    raw = planted_cohort(tmp_path)
    _, cfg = write_config(tmp_path, raw)
    core.run(cfg, results_dir=tmp_path / 'out', seed=5)
    assert 'Table 4b' in (tmp_path / 'out' / 'report.md').read_text(encoding='utf-8')
