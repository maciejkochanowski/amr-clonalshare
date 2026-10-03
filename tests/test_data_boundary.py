"""Regression checks for explicit phenotype semantics and lossless input inventory."""
import numpy as np
import pandas as pd
import pytest

from amr_clonalshare.config import Config, ConfigError, DatasetConfig, from_dict
from amr_clonalshare.io import load_dataset
from amr_clonalshare.phenotype import to_non_susceptible


def calls(values):
    return pd.DataFrame({'Strain_ID': [f'{i:03}' for i in range(len(values))],
                         'antibiotic': 'a', 'resistant_phenotype': values})


def config(tmp_path, call_text=None, mic_text=None, metadata='Strain_ID,ST\n001,A\n002,B\n'):
    (tmp_path / 'meta.csv').write_text(metadata, encoding='utf-8')
    options = dict(name='test', data_dir=str(tmp_path), metadata='meta.csv', lineage_column='ST')
    for name, content in [('phenotype', call_text), ('mic', mic_text)]:
        if content is not None:
            (tmp_path / f'{name}.csv').write_text(content, encoding='utf-8')
            options[name] = f'{name}.csv'
    return options


def test_clinical_abbreviations_r_only_and_missing_unknown_separate():
    p = to_non_susceptible(calls(['S', 'I', 'R', pd.NA, '', 'mystery']), kind='clinical_sir')
    assert p.a.iloc[:3].tolist() == [0, 0, 1]
    assert len(p) == 6
    q = p.attrs['phenotype_reading']
    assert q['n_missing_calls'] == 2
    assert q['n_unrecognized_calls'] == 1
    assert q['positive_definition'] == 'R (resistant)'


@pytest.mark.parametrize('kind,values', [('wt_nwt', ['WT', 'NWT']), ('binary', [0, 1])])
def test_explicit_nonclinical_types(kind, values):
    p = to_non_susceptible(calls(values), kind=kind)
    assert p.a.tolist() == [0, 1]
    assert p.attrs['phenotype_reading']['phenotype_kind'] == kind


def test_duplicate_policies_and_empty_agent_registry():
    d = calls(['S', 'R', ''])
    d.loc[1, 'Strain_ID'] = '000'
    d.loc[2, 'antibiotic'] = 'empty'
    with pytest.raises(ValueError, match='conflict'):
        to_non_susceptible(d, kind='clinical_sir')
    p = to_non_susceptible(d, kind='clinical_sir', duplicate_policy='drop_conflicts')
    assert list(p.columns) == ['a', 'empty']
    assert list(p.index) == ['000', '002']
    assert p.isna().all().all()
    assert p.attrs['phenotype_reading']['n_conflicting_keys'] == 1
    d['resistant_phenotype'] = ['Susceptible', 'Resistant', '']
    positive = to_non_susceptible(d, duplicate_policy='positive_wins')
    assert positive.loc['000', 'a'] == 1


def test_union_raw_universe_preserves_ids_and_metadata_inventory(tmp_path):
    options = config(tmp_path, 'Strain_ID,antibiotic,resistant_phenotype\n001,a,\nNA,empty,\n',
                     'Strain_ID,antibiotic,measurement\n002,a,1\n')
    ds = load_dataset(Config(DatasetConfig(**options, phenotype_kind='clinical_sir')))
    assert list(ds.strain_ids) == ['001', 'NA', '002']
    assert ds.panel.isna().all().all()
    assert ds.input_qc['phenotype_reading']['n_missing_calls'] == 2
    assert ds.input_qc['metadata_join']['n_metadata_only'] == 0


@pytest.mark.parametrize('table', ['phenotype', 'mic', 'metadata'])
def test_conflicting_duplicate_tables_error(tmp_path, table):
    options = config(tmp_path, 'Strain_ID,antibiotic,resistant_phenotype\n001,a,Susceptible\n',
                     'Strain_ID,antibiotic,measurement\n001,a,1\n')
    path = tmp_path / ('meta.csv' if table == 'metadata' else f'{table}.csv')
    extra = {'phenotype': '001,a,Resistant\n', 'mic': '001,a,2\n', 'metadata': '001,Z\n'}[table]
    with path.open('a') as f:
        f.write(extra)
    with pytest.raises(ConfigError, match='conflicting record'):
        load_dataset(Config(DatasetConfig(**options)))


def test_one_file_can_serve_as_metadata_and_calls(tmp_path):
    """A long table repeats the isolate; only a column the analysis reads can conflict."""
    both = ('Strain_ID,ST,antibiotic,resistant_phenotype\n'
            '001,A,a,R\n001,A,b,S\n002,B,a,S\n002,B,b,R\n')
    (tmp_path / 'both.csv').write_text(both)
    options = dict(name='one-file', data_dir=str(tmp_path), metadata='both.csv',
                   lineage_column='ST', phenotype='both.csv',
                   phenotype_id_column='Strain_ID', phenotype_kind='clinical_sir')
    ds = load_dataset(Config(DatasetConfig(**options)))
    assert list(ds.strain_ids) == ['001', '002']
    assert ds.input_qc['metadata_join']['n_identical_duplicates_collapsed'] == 2
    assert ds.input_qc['metadata_join']['n_conflicting_keys'] == 0
    assert ds.input_qc['lineage']['n_groups'] == 2


def test_a_repeated_isolate_still_conflicts_on_a_column_the_analysis_reads(tmp_path):
    both = ('Strain_ID,ST,antibiotic,resistant_phenotype\n'
            '001,A,a,R\n001,Z,b,S\n002,B,a,S\n')
    (tmp_path / 'both.csv').write_text(both)
    options = dict(name='one-file', data_dir=str(tmp_path), metadata='both.csv',
                   lineage_column='ST', phenotype='both.csv',
                   phenotype_id_column='Strain_ID', phenotype_kind='clinical_sir')
    with pytest.raises(ConfigError, match="conflicting records.*'ST'"):
        load_dataset(Config(DatasetConfig(**options)))


def test_empty_valid_cohort_and_empty_trait_prevalence(tmp_path):
    options = config(tmp_path, 'Strain_ID,antibiotic,resistant_phenotype\n')
    ds = load_dataset(Config(DatasetConfig(**options)))
    assert ds.n == 0
    assert ds.input_qc['metadata_join']['n_metadata_only'] == 2
    from amr_clonalshare.qc import trait_adequacy
    assert trait_adequacy(pd.DataFrame({'empty': [np.nan]}))['empty']['prevalence'] is None


def test_mic_wells_reject_off_panel_and_preserve_provenance(tmp_path):
    options = config(tmp_path, mic_text='Strain_ID,antibiotic,measurement\n001,a,3\n')
    with pytest.raises(ConfigError, match='outside.*wells|off.panel'):
        load_dataset(Config(DatasetConfig(**options, mic_wells={'a': [1, 2, 4]})))
    (tmp_path / 'mic.csv').write_text('Strain_ID,antibiotic,measurement\n001,a,<=2\n')
    ds = load_dataset(Config(DatasetConfig(**options, mic_wells={'a': [1, 2, 4]}, mic_units={'a': 'mg/L'})))
    assert ds.mic_join['mic_units'] == {'a': 'mg/L'}


def test_config_reads_new_fields_and_rejects_invalid_wells():
    raw = {'dataset': {'name': 'a', 'phenotype': 'p', 'metadata': 'm', 'lineage_column': 'ST',
                       'phenotype_kind': 'clinical_sir', 'ast_standard': 'supplied', 'ast_version': '2026',
                       'phenotype_source': 'lab export', 'mic_wells': {'a': [1, 2]}, 'mic_units': {'a': 'mg/L'}}}
    cfg = from_dict(raw).validate(check_files_exist=False)
    assert cfg.dataset.ast_standard == 'supplied'
    raw['dataset']['mic_wells'] = {'a': [2, 1]}
    with pytest.raises(ConfigError, match='mic_wells'):
        from_dict(raw).validate(check_files_exist=False)


@pytest.mark.parametrize('table', ['phenotype', 'mic', 'metadata'])
def test_identical_duplicates_collapse_and_drop_conflicts_keep_raw_cohort(tmp_path, table):
    options = config(tmp_path, 'Strain_ID,antibiotic,resistant_phenotype\n001,a,S\n',
                     'Strain_ID,antibiotic,measurement\n001,a,1\n', metadata='Strain_ID,ST\n001,A\n')
    path = tmp_path / ('meta.csv' if table == 'metadata' else f'{table}.csv')
    original = path.read_text()
    path.write_text(original + original.splitlines()[-1] + '\n')
    ds = load_dataset(Config(DatasetConfig(**options, phenotype_kind='clinical_sir')))
    report = ds.input_qc[{'phenotype': 'phenotype_reading', 'mic': 'mic_join', 'metadata': 'metadata_join'}[table]]
    assert report['n_identical_duplicates_collapsed'] == 1
    extra = {'phenotype': '001,a,R\n', 'mic': '001,a,2\n', 'metadata': '001,Z\n'}[table]
    path.write_text(original + extra)
    ds = load_dataset(Config(DatasetConfig(**options, phenotype_kind='clinical_sir', duplicate_policy='drop_conflicts')))
    assert ds.n == 1
    report = ds.input_qc[{'phenotype': 'phenotype_reading', 'mic': 'mic_join', 'metadata': 'metadata_join'}[table]]
    assert report['n_conflict_rows_dropped'] == 2
    if table == 'phenotype':
        assert pd.isna(ds.panel.loc['001', 'a'])
    elif table == 'metadata':
        assert ds.input_qc['lineage']['n_untyped'] == 1
    else:
        assert ds.mic.empty
        assert ds.mic_join['antimicrobials'] == ['a']


def test_semantic_duplicate_equivalence_and_clinical_i_conflict():
    d = calls(['S', 'susceptible'])
    d['Strain_ID'] = '001'
    p = to_non_susceptible(d, kind='clinical_sir')
    assert p.loc['001', 'a'] == 0
    assert p.attrs['phenotype_reading']['n_identical_duplicates_collapsed'] == 1
    d['resistant_phenotype'] = ['I', 'S']
    with pytest.raises(ValueError, match='conflicting records'):
        to_non_susceptible(d, kind='clinical_sir')


def test_nonclinical_vocabulary_never_guesses_clinical_and_records_source():
    p = to_non_susceptible(calls(['WT', 'NWT', 'R']), kind='wt_nwt',
        source='export', ast_standard='supplied scheme', ast_version='v1')
    assert pd.isna(p.a.iloc[2])
    q = p.attrs['phenotype_reading']
    assert q['positive_definition'] == 'NWT (non-wild-type)'
    assert q['n_unrecognized_calls'] == 1
    assert (q['source'], q['ast_standard'], q['ast_version']) == ('export', 'supplied scheme', 'v1')


def test_mic_blank_unknown_invalid_and_operator_conflict(tmp_path):
    options = config(tmp_path, mic_text='Strain_ID,antibiotic,measurement,op\n001,a,,\n002,a,unknown,\n003,a,0,\n004,a,inf,\n')
    ds = load_dataset(Config(DatasetConfig(**options, mic_operator_column='op')))
    assert ds.mic_join['n_missing_values'] == 1
    assert ds.mic_join['n_unparseable_values'] == 3
    assert ds.mic.measurement.isna().all()
    (tmp_path / 'mic.csv').write_text('Strain_ID,antibiotic,measurement,op\n001,a,<2,>\n')
    with pytest.raises(ConfigError, match='conflicting recorded and inline'):
        load_dataset(Config(DatasetConfig(**options, mic_operator_column='op')))


@pytest.mark.parametrize('marker', ['NULL', 'null', '#N/A', '#N/A N/A', '#NA', '-1.#IND',
                                   '-1.#QNAN', '-NaN', '-nan', '1.#IND', '1.#QNAN',
                                   '<NA>', 'N/A', 'NA', 'NaN', 'None', 'n/a', 'nan', ''])
def test_metadata_analysis_missing_markers_preserve_exact_ids(tmp_path, marker):
    options = config(tmp_path, 'Strain_ID,antibiotic,resistant_phenotype\nNULL,a,Resistant\nNA,a,Susceptible\n001,a,Resistant\n',
        metadata=f'Strain_ID,ST,period,batch,untouched\nNULL,{marker},{marker},{marker},NULL\nNA,{marker},p1,b1,NA\n001,A,p2,b2,{marker}\n')
    ds = load_dataset(Config(DatasetConfig(**options, contrast_column='period', contrast_levels=['p1', 'p2'], batch_column='batch')))
    assert list(ds.metadata.index) == ['NULL', 'NA', '001']
    assert ds.metadata.loc[['NULL', 'NA'], 'ST'].isna().all()
    assert pd.isna(ds.metadata.loc['NULL', 'period'])
    assert pd.isna(ds.metadata.loc['NULL', 'batch'])
    assert ds.metadata.loc['001', 'untouched'] == marker
    assert ds.input_qc['lineage']['n_typed'] == 1
    assert ds.input_qc['lineage']['n_untyped'] == 2
    assert ds.input_qc['lineage']['group_sizes'] == {'A': 1}


@pytest.mark.parametrize('intermediate', ['non_susceptible', 'susceptible', 'drop'])
def test_undeclared_vocabulary_reads_resistance_words_only(intermediate):
    values = ['R', 'S', 'I', 'Resistant', 'Susceptible', 'Intermediate',
              'susceptible, increased exposure', 'susceptible increased exposure',
              'Nonsusceptible', 'Non-susceptible', 'other']
    word_map = {'resistant': 1., 'susceptible': 0.,
               'intermediate': {'non_susceptible': 1., 'susceptible': 0., 'drop': np.nan}[intermediate],
               'nonsusceptible': 1., 'non-susceptible': 1.}
    expected = pd.Series(values).str.lower().map(word_map).to_numpy()
    with pytest.warns(RuntimeWarning, match='undeclared'):
        p = to_non_susceptible(calls(values), kind='undeclared', intermediate=intermediate, duplicate_policy='positive_wins')
    np.testing.assert_allclose(p.a.to_numpy(), expected, equal_nan=True)
    assert p.attrs['phenotype_reading']['n_unrecognized_calls'] == 6


def test_a_wide_call_table_reads_as_the_long_one_it_melts_to(tmp_path):
    wide = 'Strain_ID,ST,a,b\n001,A,R,S\n002,B,S,R\n'
    long = ('Strain_ID,antibiotic,resistant_phenotype\n'
            '001,a,R\n001,b,S\n002,a,S\n002,b,R\n')
    (tmp_path / 'wide.csv').write_text(wide)
    (tmp_path / 'meta.csv').write_text('Strain_ID,ST\n001,A\n002,B\n')
    (tmp_path / 'long.csv').write_text(long)
    common = dict(name='shape', data_dir=str(tmp_path), metadata='meta.csv',
                  lineage_column='ST', phenotype_kind='clinical_sir')
    from_wide = load_dataset(Config(DatasetConfig(
        phenotype='wide.csv', phenotype_agent_columns=('a', 'b'), **common)))
    from_long = load_dataset(Config(DatasetConfig(phenotype='long.csv', **common)))
    pd.testing.assert_frame_equal(from_wide.panel, from_long.panel)
    assert (from_wide.input_qc['phenotype_reading']['rows_read']
            == from_long.input_qc['phenotype_reading']['rows_read'] == 4)


def test_a_wide_table_refuses_a_column_it_does_not_have(tmp_path):
    (tmp_path / 'wide.csv').write_text('Strain_ID,ST,a\n001,A,R\n')
    (tmp_path / 'meta.csv').write_text('Strain_ID,ST\n001,A\n')
    options = dict(name='shape', data_dir=str(tmp_path), metadata='meta.csv',
                   lineage_column='ST', phenotype='wide.csv',
                   phenotype_agent_columns=('a', 'absent'))
    with pytest.raises(ConfigError, match="no column"):
        load_dataset(Config(DatasetConfig(**options)))


def test_mic_signs_as_exports_write_them_and_combination_values(tmp_path):
    text = ('Strain_ID,antibiotic,measurement\n001,a,≤0.25\n002,a,≥4\n'
            '001,b,0.5/9.5\n002,b,>=2/38\n')
    ds = load_dataset(Config(DatasetConfig(**config(tmp_path, mic_text=text))))
    m = ds.mic.set_index(['Strain_ID', 'antibiotic'])
    assert m.loc[('001', 'a'), 'measurement'] == .25
    assert m.loc[('002', 'a'), 'measurement'] == 4
    assert m.loc[('001', 'b'), 'measurement'] == .5
    assert m.loc[('002', 'b'), 'measurement'] == 2
    assert ds.mic_join['n_unparseable_values'] == 0


def test_a_censored_value_that_does_not_parse_is_refused(tmp_path):
    text = 'Strain_ID,antibiotic,measurement\n001,a,>see note\n002,a,1\n'
    with pytest.raises(ConfigError, match='censored value'):
        load_dataset(Config(DatasetConfig(**config(tmp_path, mic_text=text))))
