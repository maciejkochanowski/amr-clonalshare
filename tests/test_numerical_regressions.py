# Numerical and public-API regression tests.
import numpy as np
import pytest
from amr_clonalshare import censored as c
from amr_clonalshare.config import ConfigError, from_dict

def test_explicit_wells_never_snap_api_values():
    with pytest.raises(ValueError,match="well|panel|concentration"):
        c.intervals_from_mic([1.3],wells=[.5,1.,2.],operators=["="])

@pytest.mark.parametrize("bad", [-1.,.25,2.,np.inf])
def test_binary_interval_api_rejects_nonbinary_values(bad):
    with pytest.raises(ValueError,match="binary|zero|one|0/1"):
        c.intervals_from_binary([0.,1.,bad])

def test_binary_interval_api_preserves_missing():
    lo,hi=c.intervals_from_binary([0.,1.,np.nan])
    assert np.isnan(lo[2]) and np.isnan(hi[2])

def test_from_dict_rejects_invalid_budget_without_separate_validation():
    with pytest.raises(ConfigError,match="n_boot"):
        from_dict({"dataset":{"name":"audit","metadata":"metadata.csv","lineage_column":"lineage","phenotype":"calls.csv"},"attribution":{"n_boot":-1}})
