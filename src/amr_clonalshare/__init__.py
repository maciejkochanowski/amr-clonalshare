"""amr-clonalshare — multi-layer trait clustering for bacterial genotype data.

The package does three things, in this order of importance:

1. **Diagnoses the failure modes** that make multi-layer trait clusterings look
   more meaningful than they are: fusion collapsing onto a single layer
   (:mod:`~amr_clonalshare.influence`), the trait-absent stratum masquerading
   as an archetype and tied affinities making the answer depend on input row
   order (:mod:`~amr_clonalshare.stats`,
   :mod:`~amr_clonalshare.fusion`), co-inherited loci inflating feature
   counts, and clonal population structure
   (:mod:`~amr_clonalshare.lineage`).
2. **Provides valid post-clustering inference** for binary panels by multi-split
   feature splitting with exchangeable p-value merging
   (:mod:`~amr_clonalshare.inference`), plus an exact, closed-form
   phylogenetic convergence test for small cohorts
   (:mod:`~amr_clonalshare.archephy`).
3. **Reads a susceptibility panel at its recorded resolution** and says how
   much of it is the lineage: a clonal share and a partition attribution
   (:mod:`~amr_clonalshare.attribution`), the same share from a call, a
   recorded concentration or a censored reading
   (:mod:`~amr_clonalshare.censored`), evidence that may be re-inspected at
   any time (:mod:`~amr_clonalshare.evalues`), and a prevalence difference
   split into lineage mix and within-lineage rate
   (:mod:`~amr_clonalshare.clonality`).
4. **Clusters**: per-layer distances, similarity-network fusion, consensus
   clustering, k selection that can return "no clusters", and effect-size
   profiles (:mod:`~amr_clonalshare.core`), benchmarked against the
   baselines any such pipeline has to beat
   (:mod:`~amr_clonalshare.baselines`).
"""
from __future__ import annotations

from importlib import import_module

__version__ = "1.0.0"

__all__ = [
    "__version__", "Config", "ConfigError", "load_config", "run", "validate",
    "inference", "influence", "lineage", "baselines", "archephy", "tva",
    "clonality", "attribution", "censored", "evalues", "realised",
]

_PUBLIC_MODULES = {"inference", "influence", "lineage", "baselines", "archephy",
                   "tva", "clonality", "attribution", "censored", "evalues",
                   "realised"}
_CONFIG_EXPORTS = {"Config", "ConfigError", "load_config"}
_CORE_EXPORTS = {"run", "validate"}


def __getattr__(name: str):
    """Load the numerical stack only when a public object is first requested.

    Python imports the package before executing ``python -m
    amr_clonalshare.cli``. Eager imports here therefore initialized BLAS
    before the CLI could apply ``--threads``. PEP 562 lazy attributes preserve
    the public API while making the execution limit effective.
    """
    if name in _PUBLIC_MODULES:
        value = import_module(f".{name}", __name__)
    elif name in _CONFIG_EXPORTS:
        value = getattr(import_module(".config", __name__), name)
    elif name in _CORE_EXPORTS:
        value = getattr(import_module(".core", __name__), name)
    else:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
