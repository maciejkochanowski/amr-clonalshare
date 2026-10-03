"""Lineage structure in antimicrobial susceptibility tables: the lineage share
of calls and of the recorded MIC ordering, the lineage share of the
latent MIC ordering established by censored readings, decomposition of
prevalence differences and comparison of lineage definitions.

Every result has its own target, stated with it. Genome sequences are not
inputs.
"""
from __future__ import annotations

from importlib import import_module

__version__ = "1.0.0"

__all__ = [
    "__version__", "Config", "ConfigError", "load_config", "run",
    "attribution", "censored", "clonality", "evalues", "latent_order",
    "mic_order", "phenotype", "stats",
]

_PUBLIC_MODULES = {"attribution", "censored", "clonality", "evalues", "latent_order",
                   "mic_order", "phenotype", "stats"}
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
