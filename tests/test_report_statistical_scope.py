"""The user-visible report must retain the statistical scope of each result."""
from amr_clonalshare.cli import _summary
from amr_clonalshare.report import render_report
from amr_clonalshare.report_html import render_html_report
from amr_clonalshare.report_model import _interpretation


def test_the_bounds_on_the_unobserved_ordering_are_reported_as_bounds():
    """A call bounds the share of the ordering it was cut from; the report
    gives the bounds as bounds, model-free, with the lower limit beside them,
    and does not restate the share on a probit scale in their place."""
    record={"n_isolates":200,"seed":1,"metadata_diagnostics":{
        "lineage_column":"group","clonal_share":{"trait":{
            "kappa_adj":.4,"observed_low":.25,"observed_high":.55,
            "n_groups":10,"support":1.,"estimable":True,
            "latent_order_lower":.1,"latent_order_upper":.7,"latent_order_upper_bound":.7,
            "latent_order_upper_exact":True,"latent_order_lower_limit":.05}}}}
    for output in (render_report(record,_summary(record)),
                   render_html_report(record,_summary(record))):
        assert "0.100 to 0.700; ≥ 0.050" in output
        assert "no model for the latent" in output.lower()
        assert "lower confidence limit" in output.lower()
        assert "0.250 to 0.550" in output
    # an upper end that is certified, not attained, is marked as a bound
    record["metadata_diagnostics"]["clonal_share"]["trait"].update(
        latent_order_upper=.68, latent_order_upper_exact=False)
    for output in (render_report(record,_summary(record)),
                   render_html_report(record,_summary(record))):
        assert "0.100 to ≤ 0.700; ≥ 0.050" in output


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
