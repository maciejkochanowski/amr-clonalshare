"""Reporting contracts for the numerical implementation."""
from pathlib import Path
from amr_clonalshare.outputs import result_rows


def test_mic_primary_export_is_the_ordering_share():
    record = {'metadata_diagnostics': {'censored_share': {'per_agent': {
        'a': {'order': {'kappa_adj': .3, 'observed_low': .1, 'observed_high': .5, 'n_readings': 40,
                        'n_groups': 4, 'estimable': True,
                        'method': 'lineage share of the within-stratum mid-distribution score'}}}}}}
    rows = result_rows(record)
    assert [r['analysis'] for r in rows] == ['mic_order']
    assert (rows[0]['estimate'], rows[0]['lower'], rows[0]['upper']) == (.3, .1, .5)
    assert rows[0]['method'] == 'lineage share of the within-stratum mid-distribution score'
    assert 'no distributional form' in rows[0]['assumptions'].lower()


def test_ci_output_contract_uses_current_schema():
    text = (Path(__file__).resolve().parents[1]/'.github/workflows/ci.yml').read_text()
    assert "assert d['schema_version'] == '1.0'" in text


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
