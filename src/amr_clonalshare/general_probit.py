"""General population-liability confidence objects with bound empirical validation.

The selected fitted-nuisance bootstrap targets 95% coverage with an internal
0.04 threshold. Neither fitted-nuisance nor continuum coverage is exact.
"""
from __future__ import annotations

from contextlib import nullcontext
from dataclasses import asdict, dataclass
import hashlib
from importlib.resources import files
import json
from pathlib import Path
import time
from typing import Any

import numpy as np

from . import _population_probit_numerics as numerical
from ._general_probit_compute import (
    GeneralComputeOptions, GeneralComputeSession, JournalEngine, NodeJournal,
    digest, kernel, memory_plan, runtime_identity,
)
from ._general_probit_upper import mc_upper_bound

__all__ = ["general_probit_icc", "GeneralProbitResult", "GeneralComputeOptions", "GeneralComputeSession"]


def _protocol() -> dict:
    root = files("amr_clonalshare")
    data = root.joinpath("general_probit_protocol.json").read_bytes()
    record = json.loads(data)
    for path, expected in record["sources"].items():
        if hashlib.sha256(root.joinpath(*path.replace("\\", "/").split("/")).read_bytes()).hexdigest() != expected:
            raise RuntimeError("a source file of the general method differs from the recorded one")
    if hashlib.sha256(Path(numerical.__file__).read_bytes()).hexdigest() != record["numerical_source_sha256"]:
        raise RuntimeError("the numerical source differs from the recorded one")
    record["packaged_protocol_sha256"] = hashlib.sha256(data).hexdigest()
    record["adapter_source_sha256"] = {
        name: hashlib.sha256(root.joinpath(name).read_bytes()).hexdigest()
        for name in ("general_probit.py", "_general_probit_compute.py")}
    return record


@dataclass(frozen=True)
class GeneralProbitResult:
    """Explicit confidence object, computation status and empirical scope.

    An incomplete calculation exposes only provisional components in the
    computation record. Its public interval endpoints remain missing.
    Structural full sets and boundary-filled upper limits are completed objects.
    """
    ci_low: float | None
    ci_high: float | None
    components: list[list[float]]
    empty: bool
    confidence_kind: str
    status: str
    complete: bool
    informative: bool
    n: int
    n_groups: int
    n_repeated_groups: int
    repeat_isolate_support: float
    seconds: float
    failure_reason: str | None
    computation: dict[str, Any]
    provenance: dict[str, Any]
    nominal_level: float = .95
    method: str = "general_population_probit"
    target: str = "Gaussian population liability ICC"
    validation_status: str = "validation_pending"
    rho_hat: None = None
    prevalence_hat: None = None
    internal_alpha: float | None = None
    raw_acceptance_empty: bool = False
    boundary_fill_applied: bool = False
    monte_carlo_rejects_at_zero: bool | None = None

    def as_dict(self) -> dict:
        result = asdict(self)
        result["n_components"] = len(self.components)
        result["component_length"] = sum(hi - lo for lo, hi in self.components) if self.complete else None
        result["hull_width"] = (self.ci_high - self.ci_low if self.ci_low is not None and self.ci_high is not None else
                                0. if self.empty and self.complete else None)
        return result


def _unavailable_general_result(n, n_groups, n_repeated_groups, support, status, reason, *, sizes=(), seconds=0.):
    record = _protocol()
    return GeneralProbitResult(None, None, [], False, "unavailable", status, False, False,
                               int(n), int(n_groups), int(n_repeated_groups), float(support),
                               seconds, reason, {}, dict(protocol=record, sizes=list(map(int, sizes))),
                               validation_status=record["validation_status"])


def general_probit_icc(counts, sizes, *, seed: int = 42, case_key: str = "query",
                       compute: GeneralComputeOptions | GeneralComputeSession | None = None) -> GeneralProbitResult:
    """Compute the general method for the complete retained input geometry.

    No exact-bank registration is required. Zero repeated groups retain [0,1].
    With one repeated group, this method reports a one-sided upper limit.
    Otherwise it reports the bootstrap hull, with components retained
    as diagnostics. Fixed method settings cannot be tuned through this API.
    Its coverage target was validated empirically on a finite panel of designs;
    fitted-nuisance and continuum coverage remain approximate.
    """
    started = time.perf_counter()
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    if not isinstance(case_key, str) or not case_key:
        raise ValueError("case_key must be a nonempty, outcome-independent analysis identifier")
    if compute is not None and not isinstance(compute, (GeneralComputeOptions, GeneralComputeSession)):
        raise ValueError("compute must be GeneralComputeOptions or GeneralComputeSession")
    session = compute if isinstance(compute, GeneralComputeSession) else GeneralComputeSession(compute)
    protocol = _protocol()
    if np.asarray(counts).shape == (0,) and np.asarray(sizes).shape == (0,):
        return _unavailable_general_result(0, 0, 0, 0., "no_retained_calls", "No retained grouped counts")
    k, m = numerical.validate_counts(counts, sizes)
    n, groups = int(m.sum()), len(m)
    repeated = m >= 2
    r = int(repeated.sum())
    support = float(m[repeated].sum() / n)
    identity = dict(counts=k.tolist(), sizes=m.tolist(), B=4999, seed=int(seed), case_key=case_key,
                    protocol_sha256=protocol["packaged_protocol_sha256"], runtime=runtime_identity(),
                    adapter_source_sha256=protocol["adapter_source_sha256"],
                    stream_namespace="general_population_probit_v1")
    # Independent deterministic namespace; never touches the classical generator.
    method_seed = int(digest(["general_population_probit_v1", int(seed), case_key])[:16], 16)
    identity["method_seed"] = method_seed
    provenance = dict(protocol=protocol, query_sha256=digest(identity), runtime=identity["runtime"],
                      seed=int(seed), method_seed=method_seed, case_key=case_key,
                      complete_sizes=m.tolist(), complete_counts_sha256=digest(k.tolist()),
                      assumptions="Independent Gaussian group effects, conditional binomial counts, fixed or noninformative group sizes and ignorable selection")
    computation: dict[str, Any] = dict(cache_reused_nodes=0, dispatch_repeated_groups=r)
    components: list[list[float]] = []
    kind, status, complete, reason = "model_confidence_set", "ok", True, None
    alpha = .04 if r >= 2 else None
    raw_empty = boundary = False
    zero_rejected = None
    try:
        if r == 0:
            components = [[0., 1.]]
            kind, status = "structural_full_range", "structurally_unidentified"
            computation["scope"] = "All groups are singletons: the grouped model cannot identify ICC"
        elif r == 1:
            raw = mc_upper_bound(k, m, seed=method_seed, B=9999, alpha=.05,
                                 grid=[i / 100 for i in range(101)])
            components = [[0., raw["upper"]]]
            kind, status = "one_sided_upper_limit", raw["status"]
            raw_empty = not any(raw["grid_acceptance"])
            boundary = raw_empty
            zero_rejected = raw["monte_carlo_rejects_at_zero"]
            computation["native"] = raw
            computation["scope"] = "With one repeated group, this method reports a one-sided upper limit"
        else:
            constant = int(k.sum()) in (0, n)
            plan = memory_plan(m, session.options)
            computation["memory"] = plan
            if not constant and not plan["feasible"]:
                raise MemoryError(f"Planned table and reference arrays need at least {plan['minimum_total_bytes']} bytes; set memory_budget_mb to at least {plan['suggested_memory_budget_mb']} and use a worker with additional process overhead available")
            engine = None if constant else session.engine(m)
            journal = NodeJournal(session.options.cache_dir, identity) if session.options.cache_dir and not constant else None
            with journal.locked() if journal else nullcontext():
                delegated = JournalEngine(engine, journal) if journal else engine
                raw = kernel.infer_profile_bootstrap(k, m, B=4999, seed=method_seed,
                                                    case_key=case_key, alphas=(.04,), engine=delegated)
                if journal:
                    computation["cache_reused_nodes"] = journal.reused
                    computation["cache_query_sha256"] = journal.key
            computation["native"] = raw
            status, complete = raw["status"], bool(raw["usable"])
            selected = raw["sets"][0]
            if complete:
                components = selected["components"]
            else:
                computation["provisional_components"] = selected["components"]
                kind = "incomplete"
                reason = "One or more null calculations are unresolved; no completed confidence object is reported"
    except (MemoryError, RuntimeError, FloatingPointError, OSError) as error:
        complete = False
        status = "resource_incomplete" if isinstance(error, MemoryError) else "computation_incomplete"
        kind, reason = "incomplete", str(error)
        components = []
    low = components[0][0] if components else None
    high = components[-1][1] if components else None
    empty = complete and not components
    if empty:
        status = "empty_confidence_set"
    informative = complete and not empty and low is not None and high is not None and high - low < 1.
    return GeneralProbitResult(low, high, components, empty, kind, status, complete, informative,
                               n, groups, r, support, time.perf_counter() - started, reason,
                               computation, provenance, internal_alpha=alpha,
                               raw_acceptance_empty=raw_empty, boundary_fill_applied=boundary,
                               monte_carlo_rejects_at_zero=zero_rejected,
                               validation_status=protocol["validation_status"])
