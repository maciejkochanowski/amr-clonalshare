"""Optional grouped-binomial profile for a Gaussian population liability ICC.

This population target is distinct from the finite-collection membership
share. The fixed cutoff is tied to a recorded finite-grid calibration and
validation campaign; it is not a universal confidence-coverage guarantee.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache
import hashlib
from importlib.resources import files
import json
import time
from typing import Any

import numpy as np

from . import _population_probit_numerics as _numerics
from . import _population_mixing as _mixing

__all__ = ["population_probit_icc", "PopulationProbitResult"]


@dataclass(frozen=True)
class PopulationProbitResult:
    """A population-model result, its numerical status and calibration domain.

    ``identified`` follows the declared repeated-group/constant-outcome
    reporting policy; it is not a proof of model validity. ``informative``
    records whether a computed interval is narrower than the entire [0,1]
    range. Numerical failure has missing point and interval endpoints.
    """
    rho_hat: float | None
    ci_low: float | None
    ci_high: float | None
    prevalence_hat: float | None
    n: int
    n_groups: int
    n_repeated_groups: int
    repeat_isolate_support: float
    identified: bool
    informative: bool
    status: str
    nll: float | None
    n_eval: int | None
    seconds: float
    profile_topology_check: str
    critical_value: float
    nominal_level: float
    calibration_id: str
    calibration_sha256: str
    validation_id: str
    validation_status: str
    validation_summary_sha256: str | None
    numerical_source_sha256: str
    packaged_evidence_sha256: str
    method: str
    target: str
    domain_diagnostics: dict[str, Any]
    failure_reason: str | None = None
    #: check of the Gaussian law of lineage effects (``mixing_check``); None
    #: where no interval was computed
    mixing_check: dict | None = None

    def as_dict(self) -> dict:
        """Return an independent JSON-compatible record, without imputing values."""
        return asdict(self)


@lru_cache(maxsize=1)
def _evidence() -> dict:
    data = files("amr_clonalshare").joinpath("population_probit_validation.json").read_bytes()
    record = json.loads(data)
    record["packaged_evidence_sha256"] = hashlib.sha256(data).hexdigest()
    return record


def _domain(sizes, prevalence) -> dict:
    evidence = _evidence()
    sizes = np.asarray(sizes, dtype=int)
    ng = len(sizes)
    warnings = []
    lo, hi = evidence["tested_group_count_range"]
    if ng and not lo <= ng <= hi:
        warnings.append("group count outside the tested finite-grid range")
    maximum = int(sizes.max()) if ng else None
    if maximum is not None and maximum > evidence["tested_max_group_size"]:
        warnings.append("group size exceeds the tested maximum")
    plo, phi = evidence["tested_prevalence_range"]
    if prevalence is not None and not plo <= prevalence <= phi:
        warnings.append("fitted prevalence outside the tested finite-grid range")
    return {
        "finite_grid_only": True,
        "tested_group_count_range": list(evidence["tested_group_count_range"]),
        "tested_max_group_size": evidence["tested_max_group_size"],
        "tested_prevalence_range": list(evidence["tested_prevalence_range"]),
        "observed_group_count": ng, "observed_max_group_size": maximum,
        "observed_min_group_size": int(sizes.min()) if ng else None,
        "extrapolation_warnings": warnings,
        "range_inclusion_does_not_establish_validation": True,
        "minimum_repeated_groups_policy": _numerics.MIN_REPEATED_GROUPS,
        "implementation_max_group_size": _numerics.MAX_GROUP_SIZE,
        "scope": evidence["scope"], "assumptions": evidence["assumptions"],
    }


def _result(native: dict, sizes, reason=None, check=None) -> PopulationProbitResult:
    evidence = _evidence()
    names = ("rho_hat", "ci_low", "ci_high", "prevalence_hat", "n", "n_groups",
             "n_repeated_groups", "repeat_isolate_support", "identified", "informative",
             "status", "nll", "n_eval", "seconds")
    fields = {name: native[name] for name in names}
    fields.update({name: evidence[name] for name in (
        "critical_value", "nominal_level", "calibration_id", "calibration_sha256", "validation_id",
        "validation_status", "validation_summary_sha256", "numerical_source_sha256",
        "packaged_evidence_sha256", "method", "target")})
    fields["profile_topology_check"] = native.get("profile_topology_check", "not_applicable")
    fields["domain_diagnostics"] = _domain(sizes, native["prevalence_hat"])
    fields["failure_reason"] = reason
    fields["mixing_check"] = check
    return PopulationProbitResult(**fields)


def _unavailable_result(n: int, n_groups: int, n_repeated_groups: int,
                        support: float, status: str, reason: str, *, sizes=(), seconds=0.) -> PopulationProbitResult:
    native = dict(n=n, n_groups=n_groups, n_repeated_groups=n_repeated_groups,
                  repeat_isolate_support=support, status=status, rho_hat=None,
                  ci_low=None, ci_high=None, prevalence_hat=None, identified=False,
                  informative=False, nll=None, n_eval=None, seconds=seconds,
                  profile_topology_check="not_completed")
    return _result(native, sizes, reason)


def population_probit_icc(counts, sizes) -> PopulationProbitResult:
    """Fit grouped successes/tested counts under a Gaussian random-intercept model.

    Counts and sizes must be equal nonempty one-dimensional integer vectors,
    with 0 <= counts <= sizes and 1 <= sizes <= the implementation size cap.
    No observation is imputed. Constant outcomes and fewer than two repeated
    groups retain an explicit uninformative [0,1] set and a missing point.
    Invalid inputs raise ValueError; detected numerical failures return a
    missing point/interval with their failure reason. The jointly fitted
    marginal prevalence need not equal the raw prevalence for unequal groups.
    ``mixing_check`` compares the Gaussian law of lineage effects with an
    unrestricted one; where it rejects, the interval's coverage has not been
    established for the data.
    """
    start = time.perf_counter()
    try:
        counts, sizes = _numerics.validate_counts(counts, sizes)
    except (TypeError, OverflowError) as error:
        raise ValueError("Counts and sizes must be representable integer vectors") from error
    critical_value = _evidence()["critical_value"]
    if not np.isfinite(critical_value) or critical_value <= 0:
        raise RuntimeError("Packaged calibration cutoff must be positive and finite")
    try:
        native = _numerics.fit_profile(counts, sizes, critical_value=critical_value)
    except (RuntimeError, ValueError, FloatingPointError, OverflowError) as error:
        n = int(sizes.sum())
        repeated = sizes >= 2
        return _unavailable_result(n, len(sizes), int(repeated.sum()),
                                   float(sizes[repeated].sum()/n), "numerical_failure",
                                   str(error), sizes=sizes, seconds=time.perf_counter()-start)
    # The check's simulated datasets are drawn from a seed fixed by the data,
    # so the result stays a function of the counts and sizes alone.
    seed = int.from_bytes(hashlib.sha256(np.concatenate([counts, sizes]).astype('<i8').tobytes())
                          .digest()[:8], 'little')
    try:
        check = _mixing.mixing_check(counts, sizes, native, seed=seed)
    except (RuntimeError, ValueError, FloatingPointError, OverflowError):
        check = None
    return _result(native, sizes, check=check)
