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
    assert r['counts']=={'supplied':72,'observed':72,'available_a':72,'available_b':63,'common':63,
                         'common_scorable':63}
    assert r['arms']['a_common']['id_sha256']==r['arms']['b_common']['id_sha256']
    assert r['arms']['a_common']['n_positive']==r['arms']['b_common']['n_positive']
    d=r['difference_decomposition']
    terms=['selection_from_a','singletons_under_a','label_difference_on_scorable','singletons_under_b','selection_to_b']
    assert d['total_difference']==pytest.approx(sum(d[k] for k in terms),abs=1e-14)
    assert d['singletons_under_a']==pytest.approx(0.,abs=1e-14) and d['singletons_under_b']==pytest.approx(0.,abs=1e-14)
    assert d['selection_to_b']==pytest.approx(0.,abs=1e-14)
    assert 'not causal' in r['interpretation']
    assert r['paired_uncertainty'].startswith('95% bootstrap-t intervals, records and labels fixed')
    p=r['paired_intervals']
    assert set(p)==set(terms)|{'definition_difference_on_common','total_difference'}
    # a term that is zero by construction in every draw has an interval of zero width
    for k in ('singletons_under_a','singletons_under_b','selection_to_b'):
        assert p[k]==[0.,0.]
    for k,(lo,hi) in p.items():
        assert lo<=hi
    assert r['paired_n_boot_used']==999


def test_no_draws_leave_the_paired_intervals_uncomputed():
    ids,y,a,b=fixture()
    r=compare_lineage_definitions(y,a,b,ids=ids,seed=72,folds=3,repeats=3,n_perm=19,n_boot=0)
    assert r['paired_uncertainty']=='not computed'
    assert all(v is None for v in r['paired_intervals'].values())


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


def test_isolates_left_in_singletons_by_either_definition_leave_the_relabelling_term():
    """The relabelling term is taken on the isolates both definitions score.
    Here the second definition leaves six isolates alone, and removing them
    strands the partner of another under the first definition, which the
    repeated pruning also removes."""
    ids,y,a,b=fixture()
    b=b.copy(); a=a.copy()
    b[:6]=[f'U{i}' for i in range(6)]           # six singletons under the second definition
    a[0]='Z'; a[10]='Z'                          # a pair under the first whose one member is a singleton under the second
    r=run(ids,y,a,b)
    arms=r['arms']
    assert r['counts']['common_scorable']==arms['a_scorable']['n']==arms['b_scorable']['n']==63-7
    assert arms['a_scorable']['id_sha256']==arms['b_scorable']['id_sha256']
    for key in ('a_scorable','b_scorable'):
        assert arms[key]['share']['n_singletons_set_aside']==0
    d=r['difference_decomposition']
    middle=d['singletons_under_a']+d['label_difference_on_scorable']+d['singletons_under_b']
    assert d['definition_difference_on_common']==pytest.approx(middle,abs=1e-14)
    assert abs(d['identity_residual'])<1e-14


def test_every_arm_carries_the_digest_of_its_identifiers():
    """id_sha256 is the SHA-256 of the arm's identifiers in string order,
    written as compact JSON that keeps non-ASCII characters; the arms on one
    set of isolates share it and the others do not."""
    import hashlib
    ids,y,a,b=fixture()
    ids=np.array([f'Ł{i:03d}' for i in range(72)])
    a=a.copy(); a[:3]=None                      # three isolates without a label under the first definition
    r=run(ids,y,a,b)
    typed_a=np.array([x is not None for x in a]); typed_b=np.array([x is not None for x in b])

    def digest(mask):
        names=sorted(str(x) for x in ids[mask])
        return hashlib.sha256(json.dumps(names,ensure_ascii=False,separators=(',',':')).encode()).hexdigest()
    common=typed_a&typed_b
    expected={'a_available':digest(typed_a),'a_common':digest(common),'a_scorable':digest(common),
              'b_scorable':digest(common),'b_common':digest(common),'b_available':digest(typed_b)}
    assert {k:v['id_sha256'] for k,v in r['arms'].items()}==expected
    assert len(set(expected.values()))==3


def test_the_outer_terms_change_the_records_and_the_middle_terms_the_labels():
    ids,y,a,b=fixture(); a=a.copy(); a[:3]=None
    r=run(ids,y,a,b)
    k={key:arm['share']['kappa_adj'] for key,arm in r['arms'].items()}
    d=r['difference_decomposition']
    assert d['definition_difference_on_common']==k['b_common']-k['a_common']
    assert d['total_difference']==k['b_available']-k['a_available']
    assert d['selection_to_b']==pytest.approx(k['b_available']-k['b_common'],abs=1e-15)
    assert d['selection_to_b']!=0. and abs(d['identity_residual'])<1e-14


def test_a_lineage_of_two_under_both_definitions_is_scored():
    ids,y,a,b=fixture(); a=a.copy(); b=b.copy()
    a[[0,1]]='P'; b[[0,1]]='Q'
    assert run(ids,y,a,b)['counts']['common_scorable']==63


def test_an_arm_that_cannot_be_estimated_leaves_the_comparison_a_diagnostic():
    ids,y,a,b=fixture(); b=b.copy(); b[:63]='S0'  # one lineage under the second definition
    r=run(ids,y,a,b)
    assert r['status']=='diagnostic_only: one or more arms not estimable'
    assert not r['arms']['b_scorable']['share']['estimable']


def test_cli_reports_a_comparison_without_a_difference_as_a_diagnostic(tmp_path):
    import pandas as pd
    ids,y,a,b=fixture(); b=b.copy(); b[:63]='S0'
    p=tmp_path/'in.csv'; out=tmp_path/'out'
    pd.DataFrame({'id':ids,'y':y,'a':a,'b':b}).to_csv(p,index=False)
    assert main(['--input',str(p),'--id-column','id','--outcome','y','--lineage-a','a','--lineage-b','b',
                 '--output',str(out),'--permutations','9','--repeats','2'])==0
    report=(out/'report.md').read_text(encoding='utf-8')
    assert 'Status: diagnostic_only: one or more arms not estimable' in report
    assert '| Total difference | not available |' in report
    assert json.loads((out/'run_manifest.json').read_text())['kind']=='matched_lineage_comparison'


def test_the_centre_of_the_drawn_world_is_the_share_its_readings_estimate():
    """A lineage of the second definition that merges lineages of the first
    with different rates carries the spread of those rates within it: the
    share the paired intervals are centred on is the share the estimate
    estimates in the drawn world, not one that leaves that spread out."""
    from amr_clonalshare.attribution import clonal_share
    from amr_clonalshare.comparison import _arm_truth
    a = np.repeat(np.arange(4), 60).astype(object)
    b = np.repeat(['M0', 'M0', 'M1', 'M1'], 60).astype(object)
    p = np.repeat([0., 1., 0., 0.], 60)
    everyone = np.ones(p.size, dtype=bool)
    # merged lineage M0 reads 0 and 1 in equal parts, M1 only 0: a third of
    # the variance lies between them
    assert _arm_truth(p, b, everyone) == pytest.approx(1 / 3)
    assert _arm_truth(p, a, everyone) == pytest.approx(1.)
    estimate = clonal_share(p, b, folds=5, repeats=10, n_perm=99, n_boot=0, seed=3).kappa_adj
    assert estimate == pytest.approx(1 / 3, abs=.02)


def test_arms_that_score_the_same_records_are_one_analysis():
    """Setting aside singletons that no arm scores leaves the same analysis:
    its term is zero in the data and in every draw, and so is its interval."""
    ids = np.array([f'I{i:03d}' for i in range(90)])
    rng = np.random.default_rng(5)
    a = np.repeat([f'A{i}' for i in range(9)], 10).astype(object)
    b = np.repeat([f'S{i}' for i in range(18)], 5).astype(object)
    b[::10] = [f'solo{i}' for i in range(9)]
    y = rng.binomial(1, np.repeat(rng.uniform(.1, .9, 9), 10)).astype(float)
    r = compare_lineage_definitions(y, a, b, ids=ids, seed=4, folds=3, repeats=3, n_perm=19, n_boot=99)
    assert r['arms']['b_common']['share']['kappa_adj'] == r['arms']['b_scorable']['share']['kappa_adj']
    assert r['difference_decomposition']['singletons_under_b'] == 0.
    assert r['paired_intervals']['singletons_under_b'] == [0., 0.]
