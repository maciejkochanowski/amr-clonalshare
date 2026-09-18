import csv
import importlib
import json
import numpy as np


def api():
    spec=importlib.util.find_spec('amr_clonalshare.mic_inference_cli')
    assert spec is not None, 'The MIC interval command is not implemented'
    return importlib.import_module('amr_clonalshare.mic_inference_cli')


def write_exact(path):
    rng=np.random.default_rng(711)
    with path.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['lo','hi','lineage'])
        for g in range(8):
            effect=rng.normal()
            for value in effect+rng.normal(size=10):w.writerow([value,value,'L'+str(g)])


def test_exact_command_writes_strict_json_and_refuses_overwrite(tmp_path):
    source=tmp_path/'in.csv';output=tmp_path/'result.json';write_exact(source)
    args=[str(source),'--method','exact','--output',str(output)]
    assert api().main(args)==0
    result=json.loads(output.read_text())
    assert result['result']['method']=='exact_generalized_F'
    assert result['result']['low']<=result['result']['high']
    assert api().main(args)==2


def test_panel_method_requires_declared_edges(tmp_path):
    source=tmp_path/'in.csv';write_exact(source)
    assert api().main([str(source),'--method','bootstrap','--seed','33','--output',str(tmp_path/'r.json')])==2


def test_exact_method_refuses_interval_readings(tmp_path):
    source=tmp_path/'in.csv';source.write_text('lo,hi,lineage\n0,1,A\n0,1,A\n1,2,B\n1,2,B\n')
    assert api().main([str(source),'--method','exact','--output',str(tmp_path/'r.json')])==2


def test_written_result_follows_the_umask(tmp_path):
    import os
    source=tmp_path/'in.csv';output=tmp_path/'result.json';write_exact(source)
    mask=os.umask(0o022)
    try:
        assert api().main([str(source),'--method','exact','--output',str(output)])==0
    finally:
        os.umask(mask)
    assert output.stat().st_mode & 0o777 == 0o644


def test_panel_column_names_the_panel_of_each_record(tmp_path):
    from amr_clonalshare._mic_panels import observe_panels, _panels
    rng=np.random.default_rng(5)
    labels=np.repeat(['L'+str(g) for g in range(8)],10)
    y=rng.normal(size=8)[np.arange(8).repeat(10)]+rng.normal(size=80)
    panel=np.where(np.arange(80)%2==0,'A','B')
    edges={'A':[-1.,0.,1.],'B':[-2.,-1.,0.,1.,2.]}
    lo,hi=observe_panels(y,_panels(edges,panel,80))
    source=tmp_path/'in.csv'
    with source.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['lo','hi','lineage','laboratory'])
        for row in zip(lo,hi,labels,panel):w.writerow(row)
    common=[str(source),'--method','null-test','--rho','0.3','--seed','3','--bootstrap','19',
            '--panel-edges',json.dumps(edges)]
    assert api().main(common+['--panel-column','laboratory','--output',str(tmp_path/'r.json')])==0
    result=json.loads((tmp_path/'r.json').read_text())
    assert result['result']['bootstrap_attempted']==19
    assert api().main(common+['--panel-column','site','--output',str(tmp_path/'r2.json')])==2
    assert api().main(common+['--output',str(tmp_path/'r3.json')])==2


def test_covariate_column_enters_the_bootstrap_and_not_the_exact_pivot(tmp_path):
    rng=np.random.default_rng(9)
    labels=np.repeat(['L'+str(g) for g in range(8)],10); y=rng.normal(size=8)[np.arange(8).repeat(10)]+rng.normal(size=80)
    laboratory=np.where(np.arange(80)%2==0,'A','B')
    from amr_clonalshare._mic_panels import observe_panel
    lo,hi=observe_panel(y,[-2.,-1.,0.,1.,2.])
    source=tmp_path/'in.csv'
    with source.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['lo','hi','lineage','laboratory'])
        for row in zip(lo,hi,labels,laboratory):w.writerow(row)
    common=[str(source),'--method','null-test','--rho','0.3','--seed','3','--bootstrap','19',
            '--panel-edges','[-2,-1,0,1,2]','--covariate-column','laboratory']
    assert api().main(common+['--output',str(tmp_path/'r.json')])==0
    result=json.loads((tmp_path/'r.json').read_text())
    assert len(result['result']['unrestricted_fit']['coefficients'])==1
    assert api().main(common+['--covariate-column','lineage','--output',str(tmp_path/'r2.json')])==0
    assert len(json.loads((tmp_path/'r2.json').read_text())['result']['unrestricted_fit']['coefficients'])==8
    exact=tmp_path/'exact.csv';write_exact(exact)
    assert api().main([str(exact),'--method','exact','--covariate-column','lineage','--output',str(tmp_path/'e.json')])==2
