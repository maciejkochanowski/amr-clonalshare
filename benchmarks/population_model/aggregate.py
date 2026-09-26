#!/usr/bin/env python3
"""Summarise the population-model study and write the package's evidence.

    python -m benchmarks.population_model.aggregate calibrate OUT/lr_calibration OUT/calibration.json
    python -m benchmarks.population_model.aggregate cover OUT/<method>_<family> OUT/<name>.csv
    python -m benchmarks.population_model.aggregate benefit OUT/benefit OUT/benefit.csv
    python -m benchmarks.population_model.aggregate mixing OUT/mixing_<family> OUT/<name>.csv
    python -m benchmarks.population_model.aggregate evidence OUT src/amr_clonalshare

The rules are those of PROTOCOL.md: the fixed cut-off is the largest, over the
calibration cells, of the 95th percentile of the likelihood-ratio statistic at
the generating rho, a failed fit counting as an infinite statistic; a cell
passes when the upper Wilson 95% limit of its coverage reaches 0.95.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from benchmarks.population_model.design import SEEDS, cells

Z = 1.959963984540054


def wilson(covered: int, n: int):
    if n == 0:
        return float("nan"), float("nan")
    p = covered / n
    centre = (p + Z * Z / (2 * n)) / (1 + Z * Z / n)
    half = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / (1 + Z * Z / n)
    return centre - half, centre + half


def read(directory: Path) -> pd.DataFrame:
    frames = [pd.read_csv(p) for p in sorted(directory.glob("cell_*.csv"))]
    if not frames:
        raise SystemExit(f"no results in {directory}")
    data = pd.concat(frames, ignore_index=True)
    duplicated = data.duplicated(["cell", "replicate", "method"])
    if duplicated.any():
        raise SystemExit(f"{int(duplicated.sum())} duplicated replicates in {directory}")
    return data


#: The general method and the mixing check are run on 500 datasets a cell (PROTOCOL.md).
GENERAL_REPLICATES = 500


def _expected(data, family):
    design = {c["cell"]: c for c in cells() if c["family"] == family}
    general = bool(data.method.isin(["general", "mixing"]).all())
    for cell, frame in data.groupby("cell"):
        wanted = GENERAL_REPLICATES if general else design[cell]["replicates"]
        if set(frame.replicate) != set(range(wanted)):
            raise SystemExit(f"cell {cell}: replicates missing or extra")
    if set(data.cell) != set(design):
        raise SystemExit(f"cells missing: {sorted(set(design) - set(data.cell))}")
    return design


def calibrate(directory: Path, output: Path) -> float:
    data = read(directory)
    design = _expected(data, "calibration")
    rows = []
    for cell, frame in data.groupby("cell"):
        lr = frame.lr_at_truth.where(frame.status != "numerical_failure", np.inf).to_numpy(float)
        rank = math.ceil(0.95 * len(lr))
        rows.append(dict(cell=int(cell), pattern=design[cell]["pattern"], groups=design[cell]["groups"],
                         prevalence=design[cell]["prevalence"], rho=design[cell]["rho"], n=len(lr),
                         failures=int((frame.status == "numerical_failure").sum()),
                         unidentified=int(frame.status.isin(["constant_outcome_unidentified",
                                                             "insufficient_repeated_groups"]).sum()),
                         quantile_95=float(np.sort(lr)[rank - 1])))
    worst = max(rows, key=lambda r: r["quantile_95"])
    payload = dict(rule="largest 95th percentile of the likelihood-ratio statistic at the generating rho "
                        "over the calibration cells; a failed fit counts as infinite",
                   seed=SEEDS["calibration"], replicates=int(len(data)), cells=rows,
                   critical_value=worst["quantile_95"], attained_in_cell=worst["cell"])
    output.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: payload[k] for k in ("critical_value", "attained_in_cell", "replicates")}))
    return payload["critical_value"]


def cover(directory: Path, output: Path) -> pd.DataFrame:
    data = read(directory)
    family = data.family.iloc[0]
    design = _expected(data, family)
    rows = []
    for cell, frame in data.groupby("cell"):
        failed = frame.status == "numerical_failure"
        covered = int(frame.covered.fillna(0).astype(int).sum())
        n = len(frame)
        low, high = wilson(covered, n)
        width = (frame.ci_high - frame.ci_low)[~failed]
        c = design[cell]
        rows.append(dict(cell=int(cell), family=family, method=frame.method.iloc[0], pattern=c["pattern"],
                         groups=c["groups"], isolates=int(sum(c["sizes"])), law=c["law"],
                         prevalence=c["prevalence"], rho=c["rho"], n=n, covered=covered,
                         coverage=covered / n, wilson_low=low, wilson_high=high,
                         passes=bool(high >= 0.95), failures=int(failed.sum()),
                         unidentified=int(frame.status.isin(["constant_outcome_unidentified",
                                                             "insufficient_repeated_groups"]).sum()),
                         full_range=int(((frame.ci_low <= 0) & (frame.ci_high >= 1)).sum()),
                         median_width=float(width.median()) if len(width) else float("nan"),
                         mean_seconds=float(frame.seconds.mean())))
    table = pd.DataFrame(rows)
    table.to_csv(output, index=False)
    print(table.groupby(["family", "law"]).agg(cells=("cell", "size"), passing=("passes", "sum"),
                                               min_coverage=("coverage", "min")).to_string())
    return table


def benefit(directory: Path, output: Path) -> pd.DataFrame:
    data = read(directory)
    _expected(data, "benefit")
    ok = data.status == "ok"
    table = data[ok].groupby(["cell", "rho", "prevalence"]).agg(
        n=("replicate", "size"), rho_hat_mean=("rho_hat", "mean"), rho_hat_sd=("rho_hat", "std"),
        share_mean=("share_kappa_adj", "mean"), share_sd=("share_kappa_adj", "std"),
        latent_restatement_mean=("share_latent", "mean")).reset_index()
    table["unidentified_or_failed"] = [int((~ok & data.cell.eq(c)).sum()) for c in table.cell]
    table.to_csv(output, index=False)
    print(table.to_string())
    return table


def mixing(directory: Path, output: Path) -> pd.DataFrame:
    """Rejection rate of the mixing check, and coverage of the fixed-cut-off
    interval among the datasets the check passes and among those it rejects."""
    data = read(directory)
    family = data.family.iloc[0]
    design = _expected(data, family)
    rows = []
    for cell, frame in data.groupby("cell"):
        c = design[cell]
        checked = frame.mixing_rejected.notna()
        rejected = frame.mixing_rejected.fillna(0).astype(int).astype(bool) & checked
        covered = frame.covered.fillna(0).astype(int)
        row = dict(cell=int(cell), family=family, pattern=c["pattern"], groups=c["groups"], law=c["law"],
                   prevalence=c["prevalence"], rho=c["rho"], n=len(frame), checked=int(checked.sum()),
                   rejected=int(rejected.sum()),
                   rejection_rate=float(rejected.sum() / checked.sum()) if checked.any() else float("nan"),
                   coverage=float(covered.mean()))
        for name, mask in (("passed", checked & ~rejected), ("rejected", rejected)):
            k, n = int(covered[mask].sum()), int(mask.sum())
            low, high = wilson(k, n)
            row.update({f"n_{name}": n, f"coverage_{name}": k / n if n else float("nan"),
                        f"wilson_high_{name}": high})
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(output, index=False)
    print(table.groupby(["family", "law"]).agg(cells=("cell", "size"),
                                               mean_rejection=("rejection_rate", "mean"),
                                               max_rejection=("rejection_rate", "max"),
                                               min_coverage_passed=("coverage_passed", "min")).to_string())
    return table


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def evidence(results: Path, package: Path) -> None:
    """Write the two packaged records the population model reads at run time."""
    numerical = _sha(package / "_population_probit_numerics.py")
    calibration = json.loads((results / "calibration.json").read_text())
    table = pd.read_csv(results / "validation_lr.csv")
    gaussian = table[table.law.eq("gaussian")]
    designs = [dict(cell_id=int(r.cell), design=r.pattern, prevalence=float(r.prevalence),
                    rho_latent=float(r.rho), n_groups=int(r.groups),
                    sizes=[int(s) for s in cells()[int(r.cell)]["sizes"]],
                    coverage_all=float(r.coverage), n_all=int(r.n),
                    n_identified=int(r.n - r.unidentified - r.failures),
                    n_failures=int(r.failures), median_width=float(r.median_width))
               for r in gaussian.itertuples()]
    record = dict(
        method="grouped_binomial_probit_profile_fixed_cutoff",
        target="gaussian_population_liability_icc",
        nominal_level=0.95,
        critical_value=float(calibration["critical_value"]),
        calibration_id=f"population-model-calibration-{SEEDS['calibration']}",
        calibration_sha256=_sha(results / "calibration.json"),
        validation_id=f"population-model-validation-{SEEDS['validation']}",
        validation_status=("accepted_finite_grid_model_dependent" if bool(gaussian.passes.all())
                           else "not_accepted_finite_grid"),
        validation_summary_sha256=_sha(results / "validation_lr.csv"),
        numerical_source_sha256=numerical,
        calibration_seed=SEEDS["calibration"], validation_seed=SEEDS["validation"],
        tested_group_count_range=[int(gaussian.groups.min()), int(gaussian.groups.max())],
        tested_max_group_size=int(max(max(cells()[int(c)]["sizes"]) for c in gaussian.cell)),
        tested_prevalence_range=[float(gaussian.prevalence.min()), float(gaussian.prevalence.max())],
        scope=("Finite-grid Gaussian random-intercept model calibration and independent "
               "validation on fresh draws; not a universal 95% coverage guarantee."),
        assumptions=("Independent Gaussian group effects, conditional binomial sampling, group "
                     "sizes fixed independently of effects, and ignorable selection for the "
                     "modeled population."),
        calibration_replicates=int(calibration["replicates"]),
        validation_replicates=int(gaussian.n.sum()),
        validation_numerical_failures=int(gaussian.failures.sum()),
        validation_all_draw_coverage_range=[round(float(gaussian.coverage.min()), 3),
                                            round(float(gaussian.coverage.max()), 3)],
        validated_designs=designs)
    target = package / "population_probit_validation.json"
    target.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print("wrote", target)

    general = pd.read_csv(results / "validation_general.csv")
    gaussian = general[general.law.eq("gaussian")]
    passed = bool(gaussian.passes.all())
    scope = ("Finite-panel empirical validation of the coverage target on fresh draws from the "
             "Gaussian model. For R>=2, the bootstrap uses an internal 0.04 threshold against a 95% "
             "hull-coverage target; R1 uses a separately justified 95% upper limit. Fitted-nuisance "
             "bootstrap and continuous-parameter interpolation remain approximate, with no exact "
             "continuum guarantee.")
    target = package / "general_probit_protocol.json"
    settings = json.loads(target.read_text())
    sources = ("_general_probit/histogram_table/histogram_table.py",
               "_general_probit/profile_bootstrap/bootstrap.py", "_general_probit_upper.py")
    protocol = {key: settings[key] for key in (
        "B", "a_grid_count", "alternative_rho_count", "components", "dispatch", "internal_alpha",
        "nominal_level", "primary_object", "refinement_levels")}
    protocol.update(
        protocol_id="general_population_probit_1.0.0",
        sources={name: _sha(package / name) for name in sources},
        numerical_source_sha256=numerical,
        accepted=passed,
        validation_status="validation_accepted" if passed else "validation_not_accepted",
        coverage_scope=scope,
        validation_evidence=dict(
            validation_summary_sha256=_sha(results / "validation_general.csv"),
            validation_seed=SEEDS["validation"], validation_replicates=int(gaussian.n.sum()),
            validation_cells=int(len(gaussian)),
            coverage_range=[round(float(gaussian.coverage.min()), 3),
                            round(float(gaussian.coverage.max()), 3)],
            scope=scope))
    target.write_text(json.dumps(protocol, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("wrote", target)


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 3 or argv[0] not in ("calibrate", "cover", "benefit", "mixing", "evidence"):
        print(__doc__, file=sys.stderr)
        return 2
    step, source, target = argv[0], Path(argv[1]), Path(argv[2])
    {"calibrate": calibrate, "cover": cover, "benefit": benefit, "mixing": mixing, "evidence": evidence}[step](source, target)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
