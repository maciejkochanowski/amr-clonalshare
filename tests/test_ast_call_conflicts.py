"""Public AST calls must not depend on the order of duplicate results."""
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))

from ast_calls import parse_ast


def test_conflicting_calls_are_missing_in_both_orders():
    for field in ('drug=S,drug=R', 'drug=R,drug=S', 'drug=S,drug=I'):
        assert np.isnan(parse_ast(field)['drug'])


def test_agreeing_duplicates_and_undetermined_calls():
    assert parse_ast('drug=I,drug=R,drug=ND,other=S,other=S') == {'drug': 1.0, 'other': 0.0}
    assert parse_ast('drug=ND') == {}
