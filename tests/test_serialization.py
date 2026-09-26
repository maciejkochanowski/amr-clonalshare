"""Hostile report labels and numpy values remain portable text."""
import json
from pathlib import Path

import numpy as np

from amr_clonalshare.cli import _summary
from amr_clonalshare.report import render_report


def test_markdown_escapes_raw_html_labels_in_tables_and_identity():
    label = '<img src=x onerror="alert(1)">'
    record = {"n_isolates": 60, "seed": 1, "metadata_diagnostics": {
        "lineage_column": label, "clonal_share": {label: {
            "kappa_adj": .4, "ci_low": .2, "ci_high": .6,
            "n_groups": 3, "support": 1., "estimable": True}}}}
    report = render_report(record, _summary(record))
    assert label not in report
    assert "&lt;img" in report
    assert "| &lt;img" in report
    assert "> **" in report  # authored Markdown callouts retain their structure


def test_input_check_markdown_escapes_trait_labels():
    import pandas as pd
    from amr_clonalshare.qc import input_qc, render_markdown
    label = "<img src=x onerror=alert(1)>"
    qc = input_qc(pd.DataFrame({label: [0, 0, 0, 0]}))
    report = render_markdown(qc)
    assert label not in report
    assert "&lt;img" in report


def test_every_conversion_rule_of_the_strict_writer_is_exercised():
    """The module exists to convert numpy and pandas objects rather than to
    stringify them, so each rule is checked on the type it is for."""
    from datetime import datetime

    import pandas as pd
    import pytest

    from amr_clonalshare.jsonio import dumps, to_jsonable

    assert to_jsonable(np.bool_(True)) is True
    assert to_jsonable(np.int64(7)) == 7 and isinstance(to_jsonable(np.int64(7)), int)
    assert to_jsonable(np.float64(1.5)) == 1.5
    assert to_jsonable(np.float64("nan")) is None
    assert to_jsonable(np.float64("inf")) is None
    assert to_jsonable(np.array([[1, 2], [3, 4]])) == [[1, 2], [3, 4]]
    assert to_jsonable(pd.Timestamp("2026-09-21T10:00:00")) == \
        datetime(2026, 9, 21, 10).isoformat()
    assert to_jsonable(pd.Series({"a": 1, "b": 2})) == {"a": 1, "b": 2}
    assert to_jsonable(pd.DataFrame({"a": [1, 2]})) == [{"a": 1}, {"a": 2}]
    assert to_jsonable({1: (2, 3)}) == {"1": [2, 3]}
    assert to_jsonable([{"x": np.float32(0.5)}]) == [{"x": 0.5}]
    assert to_jsonable(Path("/tmp/out.json")) == "/tmp/out.json"
    assert json.loads(dumps({"nan": float("nan")}))["nan"] is None
    with pytest.raises(TypeError) as refused:
        to_jsonable(object())
    assert "rather than stringifying it" in str(refused.value)
