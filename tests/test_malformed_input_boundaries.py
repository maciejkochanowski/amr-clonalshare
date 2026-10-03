"""Public-boundary checks: malformed data must not become evidence."""
import numpy as np
import pytest
from amr_clonalshare.stats import benjamini_hochberg, permutation_pvalue
from amr_clonalshare.attribution import clonal_share
from amr_clonalshare.evalues import e_process


@pytest.mark.parametrize('p', [[-.01, .5], [1.01], [[.1, .2]]])
def test_step_up_rejects_invalid_probability_vectors(p):
    with pytest.raises(ValueError, match='p-value|one-dimensional'):
        benjamini_hochberg(p)


@pytest.mark.parametrize('q', [0., 1., -1., np.nan, np.inf])
def test_step_up_rejects_invalid_error_rates_even_for_empty_family(q):
    with pytest.raises(ValueError, match='q'):
        benjamini_hochberg([], q=q)


@pytest.mark.parametrize('null,observed', [([0., np.nan], 1.), ([0., 1.], np.nan), ([np.inf], 1.)])
def test_permutation_nonfinite_values_do_not_become_minimum_p(null, observed):
    with pytest.raises(ValueError, match='finite'):
        permutation_pvalue(null, observed)


@pytest.mark.parametrize('setting,value', [('folds', 2.5), ('folds', True), ('repeats', 0), ('repeats', -1), ('n_boot', -1), ('n_perm', -1)])
def test_attribution_rejects_invalid_execution_budgets(setting, value):
    kw = dict(folds=2, repeats=2, n_boot=0, n_perm=3, seed=1)
    kw[setting] = value
    with pytest.raises(ValueError, match=setting):
        clonal_share(np.tile([0., 1.], 12), np.repeat(np.arange(4), 6), **kw)


@pytest.mark.parametrize('setting,value', [('folds', 1), ('folds', 2.5), ('repeats', 0)])
def test_split_evidence_rejects_invalid_budgets(setting, value):
    kw = dict(folds=2, repeats=2, seed=1)
    kw[setting] = value
    with pytest.raises(ValueError, match=setting):
        e_process(np.tile([0., 1.], 12), np.repeat(np.arange(4), 6), **kw)


def test_unreadable_inputs_exit_with_the_input_code(tmp_path, capsys):
    from amr_clonalshare import cli
    from amr_clonalshare.config import ConfigError
    from amr_clonalshare.io import _read_table
    assert cli.main(['--config', str(tmp_path)]) == 2
    latin = tmp_path / 'latin.yaml'
    latin.write_bytes('name: caf\xe9\n'.encode('latin-1'))
    assert cli.main(['--config', str(latin)]) == 2
    (tmp_path / 'old.xlsx').write_bytes(b'\xd0\xcf\x11\xe0 not a zip workbook')
    with pytest.raises(ConfigError, match='cannot be read as an .xlsx workbook'):
        _read_table(tmp_path / 'old.xlsx', 'phenotype table', dtype=str)
    (tmp_path / 'semi.csv').write_text('id;agent;call\n1;a;R\n')
    semi = _read_table(tmp_path / 'semi.csv', 'phenotype table', dtype=str)
    assert list(semi.columns) == ['id', 'agent', 'call'] and semi.iloc[0].tolist() == ['1', 'a', 'R']
