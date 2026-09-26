"""The drafted configuration must be runnable, and every guess must be marked."""
from pathlib import Path
import pytest
import yaml

from amr_clonalshare.cli import main
from amr_clonalshare.config import from_dict
from amr_clonalshare.draft import draft_config

LONG = ('isolate_id,lineage,agent,call\n'
        '001,A,amp,R\n001,A,tet,S\n002,B,amp,S\n002,B,tet,R\n')


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding='utf-8')
    return path


def test_one_long_table_yields_a_configuration_that_loads(tmp_path):
    path = _write(tmp_path, 'both.csv', LONG)
    raw = yaml.safe_load(draft_config([path]))['dataset']
    assert raw['metadata'] == raw['phenotype'] == 'both.csv'
    assert raw['strain_id_column'] == 'isolate_id'
    assert raw['lineage_column'] == 'lineage'
    assert raw['phenotype_antibiotic_column'] == 'agent'
    assert raw['phenotype_call_column'] == 'call'
    assert raw['phenotype_kind'] == 'clinical_sir'
    from_dict({'dataset': raw}, config_path=tmp_path / 'c.yaml').validate()


def test_two_tables_are_given_their_roles(tmp_path):
    calls = _write(tmp_path, 'calls.csv', 'isolate,drug,result\n001,amp,R\n')
    meta = _write(tmp_path, 'meta.csv', 'isolate,st\n001,A\n')
    raw = yaml.safe_load(draft_config([calls, meta]))['dataset']
    assert raw['phenotype'] == 'calls.csv' and raw['metadata'] == 'meta.csv'
    assert raw['lineage_column'] == 'st'


@pytest.mark.parametrize('values,kind', [
    (['R', 'S', 'I'], 'clinical_sir'), (['WT', 'NWT'], 'wt_nwt'), (['0', '1'], 'binary')])
def test_the_recorded_vocabulary_is_read_from_the_values(tmp_path, values, kind):
    rows = ''.join(f'00{i},A,amp,{v}\n' for i, v in enumerate(values))
    path = _write(tmp_path, 'both.csv', 'isolate_id,lineage,agent,call\n' + rows)
    assert yaml.safe_load(draft_config([path]))['dataset']['phenotype_kind'] == kind


def test_a_column_that_cannot_be_guessed_is_marked_not_invented(tmp_path):
    path = _write(tmp_path, 'odd.csv', 'a,b,c\n1,2,3\n')
    text = draft_config([path])
    assert 'REPLACE_ME' in text
    assert 'Every value below is a guess' in text


def test_the_flag_prints_and_refuses_a_missing_table(tmp_path, capsys):
    path = _write(tmp_path, 'both.csv', LONG)
    assert main(['--init', str(path)]) == 0
    assert 'phenotype_kind: clinical_sir' in capsys.readouterr().out
    assert main(['--init', str(tmp_path / 'absent.csv')]) == 2
    assert 'init error' in capsys.readouterr().err


def test_a_run_still_requires_a_configuration(capsys):
    assert main([]) == 2
    assert '--config is required' in capsys.readouterr().err


WIDE = ('isolate_id,lineage,amp,tet\n'
        '001,A,R,S\n002,B,S,R\n')


def test_a_wide_sheet_is_recognised_by_the_values_of_its_columns(tmp_path):
    path = _write(tmp_path, 'wide.csv', WIDE)
    raw = yaml.safe_load(draft_config([path]))['dataset']
    assert raw['phenotype_agent_columns'] == ['amp', 'tet']
    assert 'phenotype_antibiotic_column' not in raw
    assert raw['phenotype_kind'] == 'clinical_sir'
    from_dict({'dataset': raw}, config_path=tmp_path / 'c.yaml').validate()


def test_a_draft_saved_elsewhere_finds_tables_in_other_folders(tmp_path, monkeypatch):
    from amr_clonalshare.config import load_config
    (tmp_path / 'data').mkdir()
    (tmp_path / 'meta').mkdir()
    (tmp_path / 'cfg').mkdir()
    (tmp_path / 'data' / 'calls.csv').write_text(
        'Isolate #,antibiotic,call\n1,tet,R\n2,tet,S\n')
    (tmp_path / 'meta' / 'meta.csv').write_text('Isolate #,ST\n1,1\n2,2\n')
    monkeypatch.chdir(tmp_path)
    text = draft_config([Path('data/calls.csv'), Path('meta/meta.csv')])
    (tmp_path / 'cfg' / 'draft.yaml').write_text(text)
    monkeypatch.chdir(tmp_path / 'cfg')
    cfg = load_config('draft.yaml')
    assert cfg.dataset.strain_id_column == 'Isolate #'
    assert (cfg.data_root / cfg.dataset.metadata).resolve() == (tmp_path / 'meta' / 'meta.csv').resolve()
