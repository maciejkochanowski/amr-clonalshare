"""Numerical support and explicit no-resampling contracts."""
from amr_clonalshare.config import from_dict


def test_zero_surveillance_bootstrap_budget_retains_point_only_route():
    cfg = from_dict({'dataset': {'name': 'point-only', 'metadata': 'm.csv',
        'lineage_column': 'lineage', 'phenotype': 'p.csv'},
        'surveillance': {'n_boot': 0}})
    assert cfg.surveillance.n_boot == 0

