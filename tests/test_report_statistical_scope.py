"""The user-visible report must retain the statistical scope of each result."""
from amr_clonalshare.cli import _summary
from amr_clonalshare.report import render_report
from amr_clonalshare.report_html import render_html_report
from amr_clonalshare.report_model import _interpretation


def test_latent_report_labels_the_saved_interval_as_exploratory():
    record={"n_isolates":200,"seed":1,"metadata_diagnostics":{
        "lineage_column":"group","clonal_share":{"trait":{
            "kappa_adj":.4,"ci_low":.2,"ci_high":.6,"n_groups":10,
            "support":1.,"estimable":True,"superpopulation_low":.1,
            "superpopulation_high":.8,"latent_share":.62,
            "latent_low":.3,"latent_high":.82}}}}
    for output in (render_report(record,_summary(record)),
                   render_html_report(record,_summary(record))):
        assert "exploratory" in output.lower()
        assert "probit" in output.lower()
        assert "estimated prevalence" in output.lower()
        assert "0.620" in output and "0.300 to 0.820" in output
        assert "so the interval keeps its coverage" not in output


def test_lineage_association_does_not_select_an_intervention_or_mechanism():
    rows=[dict(name="large",pt=.8,lo=.3,hi=.9,estimable=True),
          dict(name="small",pt=.2,lo=.1,hi=.4,estimable=True),
          dict(name="unresolved",pt=.01,lo=-.1,hi=.2,estimable=True)]
    output=" ".join(_interpretation(rows,{"large","small"})).lower()
    assert "association" in output
    assert "does not identify" in output
    assert "movement control" not in output
    assert "dosing" not in output
    assert "more plausibly a change in selection" not in output
