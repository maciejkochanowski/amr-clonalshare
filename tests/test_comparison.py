"""Matched-frame comparison must not confuse labels with selected records."""
import json
import numpy as np
import pytest
from amr_clonalshare.comparison import compare_lineage_definitions, main


def fixture():
    ids=np.array([f'I{i:03d}' for i in range(72)])
    a=np.repeat(['A','B','C','D','E','F'],12).astype(object)
    b=np.repeat([f'S{i}' for i in range(12)],6).astype(object)
    y=np.random.default_rng(21).binomial(1,np.repeat([.1,.2,.4,.6,.8,.9],12)).astype(float)
    b[-9:]=None
    return ids,y,a,b


def run(ids,y,a,b):
    return compare_lineage_definitions(y,a,b,ids=ids,seed=72,folds=3,repeats=3,n_perm=19)


def test_common_frame_is_identical_for_both_label_definitions():
    ids,y,a,b=fixture(); r=run(ids,y,a,b)
    assert r['counts']=={'supplied':72,'observed':72,'available_a':72,'available_b':63,'common':63}
    assert r['arms']['a_common']['id_sha256']==r['arms']['b_common']['id_sha256']
    assert r['arms']['a_common']['n_positive']==r['arms']['b_common']['n_positive']
    d=r['difference_decomposition']
    assert d['total_difference']==pytest.approx(d['selection_from_a']+d['label_difference_on_common']+d['selection_to_b'],abs=1e-14)
    assert d['selection_to_b']==pytest.approx(0.,abs=1e-14)
    assert 'not causal' in r['interpretation']
    assert r['paired_uncertainty']=='not computed; individual intervals are not intervals for differences'


def test_unique_ids_make_result_invariant_to_input_row_order():
    ids,y,a,b=fixture(); r=run(ids,y,a,b)
    order=np.random.default_rng(43).permutation(len(ids))
    s=run(ids[order],y[order],a[order],b[order])
    assert r==s


def test_missing_outcome_leaves_all_four_arms():
    ids,y,a,b=fixture(); y[:4]=np.nan; r=run(ids,y,a,b)
    assert r['counts']['observed']==68 and r['counts']['common']==59
    assert r['arms']['a_available']['share']['n']==68


def test_duplicate_ids_and_nonbinary_values_are_refused():
    ids,y,a,b=fixture(); ids[1]=ids[0]
    with pytest.raises(ValueError,match='unique'): run(ids,y,a,b)
    ids,y,a,b=fixture(); y[0]=.3
    with pytest.raises(ValueError,match='binary'): run(ids,y,a,b)


def test_empty_common_frame_is_explicit_and_strict_json():
    r=compare_lineage_definitions([0.,1.],['A',None],[None,'B'],ids=['I1','I2'],n_perm=2)
    assert r['counts']['common']==0
    assert r['difference_decomposition']['total_difference'] is None
    json.dumps(r,allow_nan=False)


def test_cli_writes_a_complete_bundle_and_refuses_overwrite(tmp_path):
    import pandas as pd
    ids,y,a,b=fixture(); p=tmp_path/'in.csv'; out=tmp_path/'out'
    pd.DataFrame({'id':ids,'y':y,'a':a,'b':b}).to_csv(p,index=False)
    args=['--input',str(p),'--id-column','id','--outcome','y','--lineage-a','a','--lineage-b','b','--output',str(out),'--permutations','9','--repeats','2']
    assert main(args)==0
    r=json.loads((out/'comparison.json').read_text())
    assert r['counts']['common']==63
    assert json.loads((out/'run_manifest.json').read_text())['status']=='complete'
    assert main(args)==2


def test_cli_refuses_settings_that_cannot_produce_a_result(tmp_path, capsys):
    import pandas as pd
    ids,y,a,b=fixture(); p=tmp_path/'in.csv'
    pd.DataFrame({'id':ids,'y':y,'a':a,'b':b}).to_csv(p,index=False)
    args=['--input',str(p),'--id-column','id','--outcome','y','--lineage-a','a','--lineage-b','b','--output',str(tmp_path/'out')]
    assert main(args+['--permutations','0'])==2
    assert '--permutations' in capsys.readouterr().err
    assert not (tmp_path/'out').exists()
