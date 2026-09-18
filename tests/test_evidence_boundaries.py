# Independent boundary tests for public evidence and label APIs.
import numpy as np
import pandas as pd
import pytest
from amr_clonalshare.evalues import combine_independent, combine_within_cohort, e_bh, e_bh_log, sequential_e_process
from amr_clonalshare.attribution import _codes, _is_untyped


def test_independent_product_recovers_after_intermediate_overflow():
    values=[1e200,1e200,1e-300,1e-150]
    assert combine_independent(values) == pytest.approx(1e-50,rel=1e-10)


def test_dependent_split_average_does_not_overflow_a_finite_convex_average():
    assert combine_within_cohort([1e308,1e308]) == pytest.approx(1e308)


def test_sequence_preserves_finite_evalues_above_log_700():
    train=np.repeat([0.,1.],200)
    held=np.repeat([0.,1.],506)
    result=sequential_e_process([train,held],[np.repeat(['a','b'],200),np.repeat(['a','b'],506)])
    assert 700 < result.log_e[-1] < np.log(np.finfo(float).max)
    assert result.e_value[-1] == pytest.approx(np.exp(result.log_e[-1]),rel=1e-12)


@pytest.mark.parametrize('alpha',[0.,-0.1,1.,1.1,np.nan,np.inf])
@pytest.mark.parametrize('function',[e_bh,e_bh_log])
def test_fdr_api_refuses_invalid_alpha(function,alpha):
    with pytest.raises(ValueError,match='alpha'):
        function([1.],alpha=alpha)


def test_within_cohort_combination_rejects_missing_evalues():
    with pytest.raises(ValueError,match='NaN|missing|finite|non-negative'):
        combine_within_cohort([1.,np.nan])


@pytest.mark.parametrize('missing',[pd.NA,pd.NaT,np.nan,None,' NULL ',' N/A '])
def test_nullable_pandas_and_text_missing_labels_are_consistent(missing):
    assert _is_untyped(missing)
    codes=_codes(['a',missing,'b',missing])
    assert codes[1] == codes[3]
    assert len(set(codes))==3

def test_by_resolution_warning_agrees_with_possible_joint_discoveries():
    from amr_clonalshare.clonality import decompose_panel
    labels=np.repeat(np.arange(4),10)
    a=pd.DataFrame(np.ones((40,10)),columns=[str(i) for i in range(10)])
    b=a*0
    result=decompose_panel(a,labels,b,labels,n_boot=199,rng=np.random.default_rng(3))
    family=result['family']
    assert family['n_within_lineage_discoveries']==10
    assert family['smallest_attainable_q'] == pytest.approx(sum(1/i for i in range(1,11))*.01)
    assert 'warning' not in family


def test_sequential_product_survives_a_first_batch_of_identical_outcomes():
    batches=[np.array([1.,1.]),np.array([0.,1.,1.]),np.array([0.,0.,1.,1.])]
    labels=[np.array(['a','b'],dtype=object),np.array(['c','a','b'],dtype=object),np.array(['a','b','c','a'],dtype=object)]
    result=sequential_e_process(batches,labels)
    assert all(np.isfinite(result.log_e))
