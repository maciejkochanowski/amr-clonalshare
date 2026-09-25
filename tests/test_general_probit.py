"""Selected general method: staged confidence contracts and preserved defaults."""
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from amr_clonalshare import core
from amr_clonalshare.config import ConfigError, PopulationProbitConfig, from_dict
from amr_clonalshare.cli import _summary
from amr_clonalshare import general_probit as api
from amr_clonalshare._general_probit_compute import (
    GeneralComputeOptions, GeneralComputeSession, JournalEngine, NodeJournal, kernel, memory_plan,
)
from amr_clonalshare.report import render_report


def test_recorded_sources_and_method_contract():
    evidence = api._protocol()
    assert evidence['accepted'] is True
    assert evidence['internal_alpha'] == .04 and evidence['B'] == 4999
    assert evidence['validation_status'] == 'validation_accepted'
    assert evidence['validation_evidence']['validation_summary_sha256']
    assert hashlib.sha256(Path(kernel.__file__).read_bytes()).hexdigest() == '558414c2bb68bec4b329f2b01671690d1da4fa08387a440f0e613e537493ea9d'


@pytest.mark.parametrize('k,m,status', [([0,1],[1,1],'structurally_unidentified'),([0,0],[2,2],'constant_outcome')])
def test_complete_full_ranges_without_engine_or_rng(k,m,status,monkeypatch):
    monkeypatch.setattr(GeneralComputeSession,'engine',lambda *a: pytest.fail('unnecessary engine'))
    before = np.random.get_state()
    result = api.general_probit_icc(k,m)
    after = np.random.get_state()
    assert all(np.array_equal(a,b) for a,b in zip(before,after))
    assert result.complete and not result.empty and not result.informative
    assert result.components == [[0.,1.]] and result.status == status
    assert result.rho_hat is None
    json.dumps(result.as_dict(),allow_nan=False)


def test_upper_boundary_fill_is_not_bootstrap_empty(monkeypatch):
    monkeypatch.setattr(api,'mc_upper_bound',lambda *a,**k: dict(upper=0.,status='empty_grid_acceptance_reported_as_zero',grid_acceptance=[False]*101,monte_carlo_rejects_at_zero=True))
    upper = api.general_probit_icc([1,0],[2,1])
    assert upper.complete and not upper.empty and upper.boundary_fill_applied
    assert upper.raw_acceptance_empty and upper.monte_carlo_rejects_at_zero
    assert upper.components == [[0.,0.]] and upper.confidence_kind == 'one_sided_upper_limit'


@pytest.mark.parametrize('usable,components,empty',[(True,[],True),(True,[[.1,.2],[.7,.9]],False),(False,[[0.,1.]],False)])
def test_empty_disconnected_and_unresolved_objects(monkeypatch,usable,components,empty):
    monkeypatch.setattr(GeneralComputeSession,'engine',lambda *a: object())
    def infer(*args,**kwargs):
        assert kwargs['B']==4999 and kwargs['alphas']==(.04,)
        return dict(status='ok' if usable else 'numerical_incomplete',usable=usable,
                    sets=[dict(components=components,empty=not components)],numerical_errors=[])
    monkeypatch.setattr(kernel,'infer_profile_bootstrap',infer)
    result=api.general_probit_icc([1,0],[2,2])
    assert result.empty == empty and result.complete == usable
    if not usable:
        assert result.ci_low is None and result.components == []
        assert result.computation['provisional_components']==components
    elif components:
        assert result.ci_low == .1 and result.ci_high == .9
        assert result.as_dict()['component_length']==pytest.approx(.3)
    else:
        assert result.as_dict()['hull_width']==0 and result.status=='empty_confidence_set'


def test_memory_failure_precedes_engine_allocation(monkeypatch):
    monkeypatch.setattr(GeneralComputeSession,'engine',lambda *a: pytest.fail('budget bypass'))
    options=GeneralComputeOptions(memory_budget_mb=2,table_cache_mb=1)
    result=api.general_probit_icc([1,0],[2,2],compute=options)
    assert result.status=='resource_incomplete' and not result.complete
    assert result.ci_low is None and not result.empty
    plan=memory_plan([20]*30,GeneralComputeOptions())
    assert plan['reference_score_matrix_bytes']==4999*2091*8
    assert plan['planned_array_bytes'] > plan['retained_table_cap_bytes']+plan['reference_score_matrix_bytes']


@pytest.mark.parametrize('kwargs',[dict(seed=True),dict(seed=-1),dict(case_key=''),dict(compute={})])
def test_arguments_are_checked_even_on_shortcuts(kwargs):
    with pytest.raises(ValueError): api.general_probit_icc([0],[1],**kwargs)


def test_query_journal_resume_identity_and_tamper(tmp_path):
    identity=dict(B=4999,sizes=[2,2],counts=[1,0],seed=713)
    journal=NodeJournal(tmp_path,identity)
    node=dict(rho=.5,status='ok',B=4999,exceedances=249,p_value=.05)
    calls=[]
    engine=SimpleNamespace(input_sizes=[2,2],alternative=None,
                           evaluate_node=lambda *a,**kw: calls.append(kw) or node)
    with journal.locked():
        delegated=JournalEngine(engine,journal)
        assert delegated.evaluate_node([1,0],.5,B=4999,seed=713,case_key='fixture')==node
    journal=NodeJournal(tmp_path,identity)
    with journal.locked():
        assert JournalEngine(engine,journal).evaluate_node([1,0],.5,B=4999,seed=713,case_key='fixture')==node
        assert journal.reused==1
    assert len(calls)==1 and journal.lock_path.exists()
    assert NodeJournal(tmp_path,{**identity,'seed':714}).key != journal.key
    stored=json.loads(journal.path.read_text());stored['payload']['nodes'][float(.5).hex()]['p_value']=.1
    journal.path.write_text(json.dumps(stored))
    with pytest.raises(RuntimeError,match='integrity'):
        with NodeJournal(tmp_path,identity).locked(): pass


def test_lock_and_interrupted_snapshot_fail_explicitly(tmp_path):
    journal=NodeJournal(tmp_path,dict(B=4999))
    with journal.locked():
        with pytest.raises(RuntimeError,match='locked'):
            with NodeJournal(tmp_path,dict(B=4999)).locked(): pass
    assert journal.lock_path.exists()
    journal.path.with_suffix('.tmp').write_text('unfinished write')
    with journal.locked(): assert journal.nodes=={}


def test_process_lock_contention_and_hard_interruption(tmp_path):
    script = '''
import sys
from amr_clonalshare._general_probit_compute import NodeJournal
try:
    with NodeJournal(sys.argv[1], dict(B=4999)).locked():
        print('owned', flush=True)
        sys.stdin.read()
except RuntimeError:
    print('contended', flush=True)
'''
    command=[sys.executable,'-c',script,str(tmp_path)]
    kwargs=dict(text=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
    journal=NodeJournal(tmp_path,dict(B=4999))
    with journal.locked():
        blocked=subprocess.run(command,input='',capture_output=True,text=True,timeout=15)
        assert blocked.returncode==0 and blocked.stdout.strip()=='contended'
    inode=journal.lock_path.stat().st_ino
    owner=subprocess.Popen(command,**kwargs)
    try:
        assert owner.stdout.readline().strip()=='owned'
        with pytest.raises(RuntimeError,match='live process'):
            with journal.locked(): pass
    finally:
        owner.kill()
        owner.communicate(timeout=15)
    with journal.locked():
        assert journal.lock_path.stat().st_ino==inode


def test_unresolved_nodes_never_enter_cache(tmp_path):
    journal=NodeJournal(tmp_path,dict(B=4999))
    engine=SimpleNamespace(input_sizes=[2,2],alternative=None,
                           evaluate_node=lambda *a,**kw:dict(rho=.5,status='unresolved',p_value=None))
    with journal.locked():
        with pytest.raises(RuntimeError,match='not cached'):
            JournalEngine(engine,journal).evaluate_node([1,0],.5,B=4999,seed=713,case_key='fixture')
    assert not journal.path.exists() and journal.nodes=={}


def test_table_preference_expands_or_contracts_without_changing_protocol():
    expanded=memory_plan([2]*2,GeneralComputeOptions(table_cache_mb=1))
    assert expanded['feasible']
    # J=441 exceeds 1 MiB at T=2091; overall memory remains sufficient.
    expanded=memory_plan([20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35],GeneralComputeOptions(table_cache_mb=1))
    assert expanded['feasible'] and expanded['retained_table_cap_bytes']>expanded['preferred_table_bytes']
    contracted=memory_plan([2,2],GeneralComputeOptions(memory_budget_mb=150,table_cache_mb=140))
    assert contracted['feasible'] and contracted['retained_table_cap_bytes']<contracted['preferred_table_bytes']
    assert expanded['reference_replicates']==contracted['reference_replicates']==4999


def test_selected_alpha_preserves_refinement_nodes():
    # Deliberately crosses both thresholds and refinement-band endpoints.
    for function in (lambda r:r/10,lambda r:.1 if r<.43 else 0.,lambda r:.03 if r<.7 else .2):
        engine=SimpleNamespace(input_sizes=np.array([2,2]),alternative=SimpleNamespace(identity=lambda:{}))
        engine.evaluate_node=lambda counts,rho,**kw:dict(rho=rho,status='ok',p_value=function(rho))
        old=kernel.infer_profile_bootstrap([1,0],[2,2],engine=engine,alphas=(.05,.04))
        selected=kernel.infer_profile_bootstrap([1,0],[2,2],engine=engine,alphas=(.04,))
        assert selected['nodes']==old['nodes']
        assert selected['sets'][0]==old['sets'][1]


def test_general_config_default_serialization_and_grouped_filtering():
    cfg=from_dict(dict(dataset=dict(name='fixture',phenotype='p.csv',metadata='m.csv',lineage_column='ST')))
    assert core._config_record(cfg)['population_probit']=={'enabled':False, 'interval_method':'fixed_cutoff'}
    with pytest.raises(ConfigError): PopulationProbitConfig(True,'finite_reference').validate()
    with pytest.raises(ConfigError): PopulationProbitConfig(True,'general',1,2).validate()
    panel=pd.DataFrame({'A':[0.,0.,1.,np.nan], 'B':[0.,1.,0.,1.]})
    rows=core._population_profiles(panel,['x','x',None,'y'],'general',seed=713)
    assert rows['A']['n']==2 and rows['A']['n_groups']==1
    assert rows['A']['n_dropped_non_finite']==1 and rows['A']['n_dropped_untyped']==1
    assert rows['A']['method']=='general_population_probit'


@pytest.mark.parametrize('accepted', [False, True])
def test_general_only_report_uses_object_and_real_counts(accepted, monkeypatch):
    protocol = dict(api._protocol(), accepted=accepted,
                    validation_status='validation_accepted' if accepted else 'validation_not_accepted')
    monkeypatch.setattr(api, "_protocol", lambda: protocol)
    rows={'A':api.general_probit_icc([0,1],[1,1]).as_dict(),
          'B':api.general_probit_icc([0,0],[2,2]).as_dict()}
    record=dict(config=dict(dataset=dict(name='fixture'),attribution=dict(enabled=False)),n_isolates=4,
                metadata_diagnostics=dict(lineage_column='ST',population_probit_icc=rows),seed=713)
    report=render_report(record,_summary(record))
    assert 'With one repeated group, this method reports a one-sided upper limit' in report
    assert '0.04 threshold' in report
    assert ('validated the coverage target' in report) == accepted
    assert ('is not validated' in report) == (not accepted)
    assert 'Fixed LR cutoff' not in report and '2 of 2 population ICC profiles computed' in report
    assert 'Antimicrobials read' in report and 'Full range retained' in report


def test_disabled_general_never_initializes_and_enabled_does_not_change_classical_rng(share_cfg,monkeypatch):
    _,cfg=share_cfg
    original=core.run(cfg,seed=53)
    monkeypatch.setattr(api,'general_probit_icc',lambda *a,**k: pytest.fail('disabled method ran'))
    disabled=core.run(replace(cfg,population_probit=PopulationProbitConfig(False,'general')),seed=53)
    assert disabled['metadata_diagnostics']==original['metadata_diagnostics']
    # Forced resource refusal exercises the enabled path without expensive fits.
    monkeypatch.undo()
    enabled=core.run(replace(cfg,population_probit=PopulationProbitConfig(True,'general',2,1)),seed=53)
    rows=enabled['metadata_diagnostics'].pop('population_probit_icc')
    assert rows and all(r['status']=='resource_incomplete' for r in rows.values())
    assert enabled['metadata_diagnostics']==original['metadata_diagnostics']


def _without_seconds(value):
    if isinstance(value,dict): return {k:_without_seconds(v) for k,v in value.items() if k!='seconds'}
    if isinstance(value,list): return [_without_seconds(v) for v in value]
    return value


def test_actual_selected_reference_and_resume_equivalence(tmp_path):
    # One tiny fixed development fixture, full selected settings; no campaign seeds.
    options=GeneralComputeOptions(cache_dir=str(tmp_path))
    session=GeneralComputeSession(options)
    first=api.general_probit_icc([1,0],[2,2],seed=713,case_key='rc3_fixed_fixture',compute=session)
    assert first.complete and first.internal_alpha==.04
    native=first.computation['native']
    second=api.general_probit_icc([1,0],[2,2],seed=713,case_key='rc3_fixed_fixture',compute=session)
    assert _without_seconds(native)==_without_seconds(second.computation['native'])
    assert second.computation['cache_reused_nodes']==len(native['nodes'])
    # Reuse the actual computed ranks to independently invert both alpha lists.
    by_rho={n['rho']:n for n in native['nodes']}
    engine=SimpleNamespace(input_sizes=np.array([2,2]),alternative=SimpleNamespace(identity=lambda:native['alternative_table_identity']))
    engine.evaluate_node=lambda counts,rho,**kw:by_rho[rho]
    dual=kernel.infer_profile_bootstrap([1,0],[2,2],B=4999,engine=engine,alphas=(.05,.04))
    assert dual['nodes']==native['nodes'] and dual['sets'][1]==native['sets'][0]


def test_actual_upper_uses_exact_selected_division_grid():
    from amr_clonalshare._general_probit_upper import mc_upper_bound
    result=api.general_probit_icc([1,0],[3,1],seed=713,case_key='rc3_upper_grid_fixture')
    raw=result.computation['native']
    selected=[i/100 for i in range(101)]
    assert raw['grid']==selected
    assert [float(x).hex() for x in raw['grid']]==[float(x).hex() for x in selected]
    reference=mc_upper_bound([1,0],[3,1],seed=result.provenance['method_seed'],B=9999,alpha=.05,grid=selected)
    assert _without_seconds(raw)==_without_seconds(reference)
    assert result.ci_high==reference['upper']


@pytest.mark.parametrize('accepted', [False, True])
@pytest.mark.parametrize('counts,sizes', [([], []), ([0, 1], [1, 1]), ([0, 0], [2, 2])])
def test_empty_and_structural_status_follow_bound_protocol(accepted, counts, sizes, monkeypatch):
    protocol = dict(api._protocol(), accepted=accepted,
                    validation_status='validation_accepted' if accepted else 'validation_not_accepted')
    monkeypatch.setattr(api, '_protocol', lambda: protocol)
    result = api.general_probit_icc(counts, sizes)
    assert result.validation_status == protocol['validation_status']
    assert result.provenance['protocol'] == protocol
    assert result.complete == bool(sizes)
    assert result.ci_low == (0. if sizes else None)
    assert result.ci_high == (1. if sizes else None)
