"""The permutation p-value under a null it claims to control.

A rejection rate at one alpha is a single point of the p-value distribution.
The check here is of the whole distribution: on cohorts drawn with no lineage
effect the p-value must be uniform on the permutation grid, and on cohorts
with a planted effect it must fall. The same generator is shipped as
``benchmarks/null_uniformity.py``, where the full run and its record live.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "benchmarks"))

from null_uniformity import run


def test_null_pvalue_is_uniform_on_the_permutation_grid():
    r = run(n_rep=300, n_power=30, n=160, n_lineages=30, seed=7)
    null = r["null"]
    assert null["uniform_on_the_grid"], null["grid_gap_p_value"]
    assert null["level"]["0.05"]["holds"], null["level"]["0.05"]
    assert abs(null["kappa_adj_mean"]) < 0.02
    power = r["power_at_0.05"]
    assert power["0.30"]["power"] > power["0.00"]["power"] + 0.5
