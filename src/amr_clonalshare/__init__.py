"""amr-clonalshare — how much of an antimicrobial resistance phenotype the
lineages of a bacterial collection carry.

From a lineage label and a susceptibility call or a recorded dilution, with
no sequence data read: the lineage-membership share of a collection and the
design-corrected component ratio beside it, with a superpopulation share for
the species (:mod:`~amr_clonalshare.attribution`,
:mod:`~amr_clonalshare.realised`, restated on the liability scale by
:mod:`~amr_clonalshare.latent`); the same share read from an interval-censored
panel (:mod:`~amr_clonalshare.censored`); evidence that may be re-inspected at
any time, including a sequential e-process over intakes
(:mod:`~amr_clonalshare.evalues`); and a prevalence difference split into
lineage mix and within-lineage rate (:mod:`~amr_clonalshare.clonality`). Each
sits behind a gate that refuses where the collection cannot identify the
quantity or the residuals contradict the law an interval assumes.
:func:`~amr_clonalshare.run` reads a configuration that names the two
columns, writes the record, and stops.
"""
from __future__ import annotations

from importlib import import_module

__version__ = "1.0.0"

__all__ = [
    "__version__", "Config", "ConfigError", "load_config", "run",
    "attribution", "realised", "latent", "censored", "clonality", "evalues",
    "phenotype", "stats",
]

_PUBLIC_MODULES = {"attribution", "realised", "latent", "censored", "clonality",
                   "evalues", "phenotype", "stats"}
_CONFIG_EXPORTS = {"Config", "ConfigError", "load_config"}
_CORE_EXPORTS = {"run"}


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
