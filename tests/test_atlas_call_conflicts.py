"""Public AST calls must not depend on the order of duplicate results."""
import importlib
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))


@pytest.mark.parametrize('module', ['atlas_cross_species', 'vet_atlas', 'resolution_atlas', 'verify_vet_claims'])
def test_conflicting_calls_are_missing_in_both_orders(module):
    parse = importlib.import_module(module).parse_ast
    for field in ('drug=S,drug=R', 'drug=R,drug=S', 'drug=S,drug=I'):
        assert np.isnan(parse(field)['drug'])


@pytest.mark.parametrize('module', ['atlas_cross_species', 'vet_atlas', 'resolution_atlas', 'verify_vet_claims'])
def test_agreeing_duplicates_and_undetermined_calls(module):
    parse = importlib.import_module(module).parse_ast
    assert parse('drug=I,drug=R,drug=ND,other=S,other=S') == {'drug':1.0, 'other':0.0}
    assert parse('drug=ND') == {}


def test_low_prevalence_is_not_reported_as_no_variance():
    from agent_screen import screen_agent
    y = np.zeros(200); y[0] = 1
    result = screen_agent(y, ['a'] * 100 + ['b'] * 100)
    assert not result.analyse
    assert 'prevalence' in result.reason
    assert result.reason != 'no variance'
