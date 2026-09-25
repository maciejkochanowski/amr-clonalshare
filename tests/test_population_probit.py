"""A separate grouped-population target must not alter classical results."""
from dataclasses import replace
import importlib
import json
import subprocess
import sys
from types import SimpleNamespace
from pathlib import Path
import hashlib

import numpy as np
import pandas as pd
import pytest

from amr_clonalshare import core
from amr_clonalshare.config import ConfigError, from_dict
from amr_clonalshare.cli import _summary
from amr_clonalshare.report import render_report
from amr_clonalshare.report_html import render_html_report


def api():
    return importlib.import_module("amr_clonalshare.population_probit")


def test_public_lazy_api_and_separate_target():
    import amr_clonalshare
    assert "population_probit_icc" in dir(amr_clonalshare)
    result = amr_clonalshare.population_probit_icc([0, 1, 3, 4], [4] * 4)
    assert isinstance(result, amr_clonalshare.PopulationProbitResult)
    rec = result.as_dict()
    assert rec["target"] == "gaussian_population_liability_icc"
    assert 0 <= result.ci_low <= result.rho_hat <= result.ci_high <= 1
    assert result.critical_value > 3.841458820694124
    assert rec["calibration_sha256"] and rec["validation_id"]
    assert rec["profile_topology_check"] == "17_point_profile_scan_passed"
    assert "clonal_share" not in rec
    json.dumps(rec, allow_nan=False)
    output = subprocess.check_output([sys.executable, "-c", "import amr_clonalshare,sys; assert 'numpy' not in sys.modules; assert 'scipy' not in sys.modules"], text=True)
    assert output == ""


@pytest.mark.parametrize("counts,sizes,status", [
    ([0, 0], [10, 10], "constant_outcome_unidentified"),
    ([0, 1], [1, 1], "insufficient_repeated_groups"),
    ([0, 1], [1, 3], "insufficient_repeated_groups"),
])
def test_unidentified_sets_remain_explicit(counts, sizes, status):
    r = api().population_probit_icc(counts, sizes)
    assert r.status == status and not r.identified and not r.informative
    assert r.rho_hat is None and (r.ci_low, r.ci_high) == (0, 1)


@pytest.mark.parametrize("k,m", [([], []), ([1], [0]), ([2], [1]), ([.5], [2]), ([0], [1e20]), ([np.nan], [1]), ([0], [10**1000])])
def test_invalid_counts_are_rejected(k, m):
    with pytest.raises(ValueError):
        api().population_probit_icc(k, m)


@pytest.mark.parametrize("error_type", [RuntimeError, ValueError, FloatingPointError, OverflowError])
def test_numerical_failure_has_no_invented_interval(monkeypatch, error_type):
    module = api()
    def fail(*args, **kwargs):
        raise error_type("optimizer failed")
    monkeypatch.setattr(module._numerics, "fit_profile", fail)
    r = module.population_probit_icc([1, 2], [3, 3])
    assert r.status == "numerical_failure"
    assert r.rho_hat is None and r.ci_low is None and r.ci_high is None
    assert r.failure_reason == "optimizer failed"
    json.dumps(r.as_dict(), allow_nan=False)


def test_domain_diagnostics_do_not_claim_a_continuous_validation_region():
    r = api().population_probit_icc([0, 0], [1001, 1001])
    assert r.domain_diagnostics["extrapolation_warnings"]
    assert r.domain_diagnostics["tested_group_count_range"] == [15, 200]
    assert r.domain_diagnostics["tested_max_group_size"] == 161
    assert r.domain_diagnostics["finite_grid_only"] is True


def test_config_is_disabled_by_default_and_no_cutoff_override(share_cfg):
    _, cfg = share_cfg
    assert cfg.population_probit.enabled is False
    raw = {"dataset": {"name": "a", "metadata": "m.csv", "phenotype": "p.csv", "lineage_column": "g"},
           "population_probit": {"enabled": True}}
    assert from_dict(raw).population_probit.enabled
    raw["population_probit"]["critical_value"] = 1
    with pytest.raises(ConfigError, match="unknown key"):
        from_dict(raw)
    raw["population_probit"] = {"enabled": "false"}
    assert from_dict(raw).population_probit.enabled is False
    raw["population_probit"] = {"enabled": "not-a-boolean"}
    with pytest.raises(ConfigError):
        from_dict(raw)


def test_enable_does_not_change_any_classical_result_or_rng_stream(share_cfg):
    _, cfg = share_cfg
    off = core.run(cfg, seed=81)
    on_cfg = replace(cfg, population_probit=replace(cfg.population_probit, enabled=True))
    on = core.run(on_cfg, seed=81)
    extra = on["metadata_diagnostics"].pop("population_probit_icc")
    assert extra and on["metadata_diagnostics"] == off["metadata_diagnostics"]
    assert off["config"]["population_probit"] == {"enabled": False, "interval_method": "fixed_cutoff"}
    assert on["config"]["population_probit"] == {"enabled": True, "interval_method": "fixed_cutoff"}


def test_exact_subset_and_failures_are_per_agent(share_cfg, monkeypatch):
    _, cfg = share_cfg
    cfg = replace(cfg, population_probit=replace(cfg.population_probit, enabled=True),
                  attribution=replace(cfg.attribution, enabled=False),
                  evidence=replace(cfg.evidence, enabled=False),
                  surveillance=replace(cfg.surveillance, enabled=False))
    panel = pd.DataFrame({"a": [0, 1, np.nan, 1, 0, 1], "b": [np.nan]*4 + [0, 1],
                          "invalid": [2, 0, 0, 1, 0, 1]})
    raw = ["x", "x", "y", "y", None, "NA"]
    ds = SimpleNamespace(metadata=pd.DataFrame({cfg.dataset.lineage_column: raw}))
    rec = core._metadata_diagnostics(cfg, ds, panel.index, np.random.default_rng(2), X_df=panel)
    rows = rec["population_probit_icc"]
    a = rows["a"]
    assert (a["n"], a["n_dropped_non_finite"], a["n_dropped_untyped"]) == (3, 1, 2)
    assert a["n_groups"] == 2 and a["n_repeated_groups"] == 1
    assert rows["b"]["status"] == "no_retained_calls"
    assert rows["invalid"]["status"] == "invalid_binary_input"
    assert rows["invalid"]["ci_low"] is None


def test_profile_only_reports_have_distinct_summary_and_safe_labels(share_cfg):
    _, cfg = share_cfg
    cfg = replace(cfg, population_probit=replace(cfg.population_probit, enabled=True),
                  attribution=replace(cfg.attribution, enabled=False))
    rec = core.run(cfg, seed=13)
    profile = rec["metadata_diagnostics"]["population_probit_icc"]
    first = next(iter(profile))
    profile['<img src=x onerror="alert(1)">'] = profile.pop(first)
    for output in (render_report(rec, _summary(rec)), render_html_report(rec, _summary(rec))):
        assert "Population liability ICC" in output
        assert "Classical attribution was disabled" in output
        assert "no antimicrobial read" not in output
        assert "finite" in output.lower() and "not a universal" in output
        assert "<img src=x" not in output
        assert "population" in output.lower()


def test_report_never_fills_failure_interval():
    module = api()
    r = module._unavailable_result(4, 2, 2, 1.0, "numerical_failure", "diagnostic failed").as_dict()
    rec = {"n_isolates": 4, "seed": 1, "config": {"attribution": {"enabled": False}},
           "metadata_diagnostics": {"population_probit_icc": {"bad": r}}}
    for output in (render_report(rec, _summary(rec)), render_html_report(rec, _summary(rec))):
        assert "numerical_failure" in output
        assert "diagnostic failed" in output


def test_profile_only_global_terms_respect_target_and_testing_scope():
    r = api().population_probit_icc([1, 2, 3, 4], [5]*4).as_dict()
    rec = {"n_isolates": 20, "seed": 1, "config": {"attribution": {"enabled": False}},
           "metadata_diagnostics": {"population_probit_icc": {"a": r}}}
    for output in (render_report(rec, _summary(rec)), render_html_report(rec, _summary(rec))):
        assert "not for a fresh draw" not in output
        assert "Population-model interval" in output
        assert "fixed-look" in output and "sequential" in output
        assert "null-model assumptions" in output
        assert "cannot identify the quantity" not in output


def test_profile_only_summary_uses_profile_agent_and_group_counts():
    from amr_clonalshare.report_model import build_report
    a = api().population_probit_icc([0, 0], [5, 5]).as_dict()
    b = api().population_probit_icc([0]*4, [5]*4).as_dict()
    rec = {"n_isolates": 20, "seed": 1, "config": {"attribution": {"enabled": False}},
           "metadata_diagnostics": {"lineage_column": "g", "population_probit_icc": {"a": a, "b": b}}}
    report = build_report(rec, _summary(rec))
    rows = next(block.rows for block in report.sections[0].blocks if block.kind == "kv")
    pairs = dict(rows)
    assert pairs["Antimicrobials read"] == "2"
    assert "2 to 4" in pairs["Cohort"] and "agent-specific" in pairs["Cohort"]


def test_recorded_numerical_source_and_packaged_evidence_match():
    module = api()
    record = module._evidence()
    assert hashlib.sha256(Path(module._numerics.__file__).read_bytes()).hexdigest() == record["numerical_source_sha256"]
    assert record["calibration_id"] != record["validation_id"]
    assert record["calibration_seed"] != record["validation_seed"]
    returned = module.population_probit_icc([0, 0], [2, 2]).as_dict()
    returned["domain_diagnostics"]["tested_group_count_range"][0] = -100
    assert module._evidence()["tested_group_count_range"] == [15, 200]


def test_direct_config_rejects_truthy_non_boolean():
    from amr_clonalshare.config import PopulationProbitConfig
    with pytest.raises(ConfigError):
        PopulationProbitConfig(enabled="false").validate()


def test_no_interval_native_path_and_lazy_large_boundary_rules():
    module = api()._numerics
    fit = module.fit_profile([1, 2, 4], [5, 5, 5], truth_rho=.2, compute_interval=False)
    assert fit["ci_low"] is None and fit["lr_at_truth"] >= 0
    module._normal_rule.cache_clear()
    module._conditional_rule.cache_clear()
    for rho in (0., 1.):
        module.selected_count_probabilities(1000000, 0., rho, (0, 500000, 1000000))
    assert module._normal_rule.cache_info().currsize == 0
    assert module._conditional_rule.cache_info().currsize == 0
    for args in ((0, 0., .5), (10001, 0., .5)):
        with pytest.raises(ValueError):
            module.count_probabilities(*args)
    for indices in ((-1,), (.5,), (np.nan,)):
        with pytest.raises(ValueError):
            module.selected_count_probabilities(2, 0., .5, indices)
    for cutoff in (0., np.inf):
        with pytest.raises(ValueError):
            module.fit_profile([0, 1], [1, 1], critical_value=cutoff)


def test_optimizer_and_disconnected_set_failures_are_exposed(monkeypatch):
    module = api()._numerics
    with monkeypatch.context() as patch:
        patch.setattr(module, "minimize_scalar", lambda *a, **k: SimpleNamespace(success=False, fun=0., x=.5))
        with pytest.raises(RuntimeError, match="Nuisance"):
            module.ProfileLikelihood([1, 2], [3, 3]).profile(.5)
        ll = module.ProfileLikelihood([1, 2], [3, 3])
        patch.setattr(ll, "profile", lambda r: ((r-.4)**2, 0.))
        with pytest.raises(RuntimeError, match="ICC interval"):
            ll.mle()
    class Disconnected:
        def __init__(self, *args):
            self.n_eval = 0
        def mle(self):
            return .2, 0., 0.
        def profile(self, r):
            return 100*min((r-.2)**2, (r-.8)**2), 0.
    with monkeypatch.context() as patch:
        patch.setattr(module, "ProfileLikelihood", Disconnected)
        with pytest.raises(RuntimeError, match="Disconnected"):
            module.fit_profile([1, 2], [3, 3])
    with monkeypatch.context() as patch:
        patch.setattr(module, "brentq", lambda *a, **k: .5)
        with pytest.raises(RuntimeError, match="endpoint"):
            module.fit_profile([1, 2], [3, 3])
