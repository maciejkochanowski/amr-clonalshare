"""MIC readings drawn as a wild-type / non-wild-type mixture, the shape of real
MIC data, read by the Gaussian MIC interval and by the population model on the
non-wild-type calls.

Latent log2 MIC of isolate i in lineage g:
    Y = mu_wt + DELTA * Z + v_g + e,   e ~ N(0, S^2),  v_g ~ N(0, SV^2)
    Z = 1{u_g + r > c},  u_g ~ N(0, rho_z),  r ~ N(0, 1 - rho_z),  P(Z = 1) = pi
The lineage sizes are those of the S. suis collection. Readings are taken on
one panel with cut points -5 to 3; an isolate is called non-wild-type when its
reading lies wholly above the cut point nearest mu_wt + DELTA / 2.

Two targets, each computed from the model statement by quadrature:
    rho_np  Var_g(E[Y | g]) / Var(Y), the lineage share of latent MIC variance,
            the quantity a variance-ratio reading of the MIC interval claims;
    rho_z   the liability ICC of non-wild-type membership, the target of the
            population model on the calls.

Each dataset records the calibrated MIC test at rho_np (coverage), the check of
the Gaussian model at the estimate (``shape_check.check_dataset``), and the
fixed-cut-off population interval on the calls with its mixing check. The
population interval uses the critical value of the campaign's own calibration,
passed as ``--critical-value``.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.special import roots_hermitenorm
from scipy.stats import norm

from amr_clonalshare._mic_panels import observe_panel
from amr_clonalshare._population_mixing import mixing_check
from amr_clonalshare._population_probit_numerics import fit_profile
from amr_clonalshare.mic_null_bootstrap import null_bootstrap_test
from benchmarks.mic_inference.shape_check import check_dataset

ROOT_SEED = 20261003
SIZES = (161, 60, 57, 50, 44, 42, 40, 36, 31, 22, 20, 17, 13, 11, 11, 10,
         8, 7, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 1)
EDGES = np.arange(-5., 4.)
S, SV, DELTA = .8, .3, 5.


def mixture_design_grid():
    cells = []
    for mu_wt in (-5.5, -1.0):          # wild-type mode below the panel, or inside it
        for pi in (.3, .6):
            for rho_z in (.3, .7):
                cells.append(dict(cell=len(cells), mu_wt=mu_wt, pi=pi, rho_z=rho_z))
    return cells


def rho_np(cell, nodes=200):
    """Lineage share of latent MIC variance, by Gauss-Hermite quadrature over u_g."""
    x, w = roots_hermitenorm(nodes)
    w = w / w.sum()
    c = norm.ppf(1 - cell["pi"])
    q = norm.cdf((np.sqrt(cell["rho_z"]) * x - c) / np.sqrt(1 - cell["rho_z"]))
    between = DELTA ** 2 * (w @ q ** 2 - (w @ q) ** 2) + SV ** 2
    within = S ** 2 + DELTA ** 2 * (w @ (q * (1 - q)))
    return float(between / (between + within))


def generate(cell, replicate):
    rng = np.random.default_rng(np.random.SeedSequence([ROOT_SEED, cell["cell"], replicate]))
    sizes = np.asarray(SIZES)
    labels = np.repeat(np.arange(len(sizes)), sizes)
    c = norm.ppf(1 - cell["pi"])
    u = rng.normal(0, np.sqrt(cell["rho_z"]), len(sizes))[labels]
    z = (u + rng.normal(0, np.sqrt(1 - cell["rho_z"]), len(labels)) > c).astype(int)
    y = cell["mu_wt"] + DELTA * z + rng.normal(0, SV, len(sizes))[labels] + rng.normal(0, S, len(labels))
    a, b = observe_panel(y, EDGES)
    seed = int(np.random.SeedSequence([ROOT_SEED, cell["cell"], replicate, 991]).generate_state(1)[0])
    return a, b, labels, z, seed


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cell", type=int, required=True)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--replicates", type=int, default=10)
    parser.add_argument("--critical-value", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    cell = mixture_design_grid()[args.cell]
    target = rho_np(cell)
    sizes = np.asarray(SIZES)
    threshold = EDGES[np.argmin(np.abs(EDGES - (cell["mu_wt"] + DELTA / 2)))]
    rows = []
    for rep in range(args.start, args.start + args.replicates):
        a, b, labels, z, seed = generate(cell, rep)
        row = dict(cell=args.cell, replicate=rep, rho_np=target, rho_z=cell["rho_z"], mu_wt=cell["mu_wt"],
                   pi=cell["pi"])
        started = time.perf_counter()
        try:
            test = null_bootstrap_test(a, b, labels, rho=target, panel_edges=EDGES, n_boot=199, seed=seed)
            row.update(g_reportable=bool(test.reportable), g_covered=bool(test.accepted),
                       g_rho_hat=test.unrestricted_fit.rho, g_redrawn=test.bootstrap_redrawn)
        except (ValueError, ArithmeticError) as exc:
            row.update(g_reportable=False, g_covered=False, g_refusal=f"{type(exc).__name__}: {exc}")
        row["g_seconds"] = time.perf_counter() - started
        row.update({f"shape_{k}": v for k, v in
                    check_dataset(a, b, labels, EDGES, None, None, seed).items()})
        calls = (a >= threshold).astype(int)
        counts = np.bincount(labels, weights=calls, minlength=len(sizes)).astype(int)
        row["misclassified"] = float(np.mean(calls != z))
        started = time.perf_counter()
        try:
            native = fit_profile(counts, sizes, critical_value=args.critical_value, truth_rho=cell["rho_z"])
            row.update(p_status=native["status"], p_rho_hat=native["rho_hat"], p_low=native["ci_low"],
                       p_high=native["ci_high"],
                       p_covered=bool(native["ci_low"] <= cell["rho_z"] <= native["ci_high"]))
            mix = mixing_check(counts, sizes, native, seed=seed)
            if mix is not None:
                row.update(p_mixing_p=mix["p_value"], p_mixing_rejected=bool(mix["rejected"]))
        except (RuntimeError, ValueError, FloatingPointError, OverflowError, ArithmeticError) as exc:
            row.update(p_status="numerical_failure", p_covered=False, p_refusal=f"{type(exc).__name__}: {exc}")
        row["p_seconds"] = time.perf_counter() - started
        rows.append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    partial = args.output.with_suffix(".partial")
    partial.write_text(json.dumps(rows, indent=1, allow_nan=True, default=float))
    partial.replace(args.output)


if __name__ == "__main__":
    main()
