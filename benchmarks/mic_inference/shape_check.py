"""Size and power of the check of the Gaussian MIC model.

Every dataset is analysed as the calibrated interval analyses it: the
unrestricted fit, then the parametric bootstrap at the fitted rho, which also
calibrates the two statistics of ``_mic_shape`` (199 simulated datasets). One
row per dataset records both p-values and the decision.

``--source null``       interval-reading cells of ``calibrate_null.py``: the
                        Gaussian core designs (size) and the t4 and
                        contaminated-normal residuals (power);
``--source extension``  cells of ``calibrate_extension.py``: heavy censoring,
                        laboratory effects (size), two-component residuals
                        within lineages (power);
``--source plasmode``   datasets simulated from the Gaussian model fitted to one
                        S. suis agent, on its recorded design (size on the
                        design of the example).

The datasets are those of the coverage campaigns (same generator, root seed and
replicate numbers), so a row joins its coverage row on (cell, replicate) or
(agent, replicate), and coverage can be read by the check's decision. The
WT/non-WT mixture designs, where the check is meant to reject, are run by
``calibrate_mixture.py``, which calls :func:`check_dataset`.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from amr_clonalshare._mic_panels import _validated_problem
from amr_clonalshare.mic_null_bootstrap import _unrestricted_fit, calibrate_null

BOOTSTRAP = 199


def check_dataset(a, b, labels, panel_edges, panel, covariate, seed, n_boot=BOOTSTRAP):
    """The shape check of one dataset, as the interval computes it at the estimate."""
    row = dict(n=int(len(a)))
    started = time.perf_counter()
    try:
        problem, panels, exact = _validated_problem(a, b, labels, panel_edges, panel, covariate)
        fit = _unrestricted_fit(problem, 24)
        test = calibrate_null(problem, panels, exact, rho=fit.rho, n_boot=n_boot, seed=seed,
                              full_fit=fit, shape_check=True)
        row.update(rho_hat=fit.rho, redrawn=test.bootstrap_redrawn, failed=test.bootstrap_failed)
        check = test.shape_check
        if check is None:
            row.update(checked=False)
        else:
            row.update(checked=True, **{k: check[k] for k in (
                "marginal_deviance", "conditional_pearson", "p_marginal", "p_conditional", "p_value",
                "rejected", "simulated_datasets")})
    except (ValueError, ArithmeticError) as exc:
        row.update(checked=False, refusal=f"{type(exc).__name__}: {exc}")
    row["seconds"] = time.perf_counter() - started
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source", choices=["null", "extension", "plasmode"], required=True)
    parser.add_argument("--cell", type=int, default=0)
    parser.add_argument("--agent")
    parser.add_argument("--inputs", type=Path, help="directory written by prepare_empirical.py")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--replicates", type=int, default=10)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    rows = []
    if args.source == "plasmode":
        from benchmarks.mic_inference.calibrate_panels import plasmode_dataset, plasmode_state
        state = plasmode_state(args.inputs, args.agent)
        for rep in range(args.start, args.start + args.replicates):
            a, b, seed = plasmode_dataset(state, rep)
            rows.append(dict(source="plasmode", agent=args.agent, replicate=rep, rho=state["fit"].rho,
                             **check_dataset(a, b, state["labels"], state["record"]["panel_edges"],
                                             state["panel"], state["country"], seed)))
    else:
        if args.source == "null":
            from benchmarks.mic_inference.calibrate_null import generate, null_design_grid as grid
        else:
            from benchmarks.mic_inference.calibrate_extension import extension_design_grid as grid, generate
            import benchmarks.mic_inference.calibrate_extension as extension
        cell = grid()[args.cell]
        if cell["reading"] == "exact":
            parser.error("the check reads interval readings; choose an interval-reading cell")
        for rep in range(args.start, args.start + args.replicates):
            a, b, labels, edges, seed = generate(cell, rep)
            covariate = extension._CURRENT["covariate"] if args.source == "extension" else None
            rows.append(dict(source=args.source, cell=args.cell, replicate=rep, rho=cell["rho"],
                             domain=cell["domain"], residual=cell["residual"], reading=cell["reading"],
                             **check_dataset(a, b, labels, edges, None, covariate, seed)))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    partial = args.output.with_suffix(".partial")
    partial.write_text(json.dumps(rows, indent=1, allow_nan=True, default=float))
    partial.replace(args.output)


if __name__ == "__main__":
    main()
