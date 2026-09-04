"""
schemas.py — pandera validation schemas for amr-clonalshare inputs.

Provides the **fail-loud** structural contract between raw CSV inputs (and the
loaded :class:`amr_clonalshare.io.Dataset`) and the downstream multi-layer
trait-clustering pipeline. The pipeline assumes every loaded ``wide_binary``
layer is strictly 0/1 with a strain-indexed table; this module makes that
assumption explicit and surfaces violations early with actionable error
messages, instead of letting them silently corrupt Jaccard / SNF / consensus
computations downstream.

Two schemas
-----------
:func:`make_layer_schema`
    Factory returning a :class:`pandera.pandas.DataFrameSchema` for one binary
    layer. Per-layer policy is supplied via :class:`LayerPolicy` (e.g. AMR
    binary 0/1, virulence binary, capsule one-hot — capsule may optionally be
    enforced as "at most one 1 per row" via :attr:`LayerPolicy.one_hot`).

:func:`make_metadata_schema`
    Factory returning a :class:`pandera.pandas.DataFrameSchema` for an optional
    isolate-metadata table, with a **foreign-key check** that every metadata
    row's isolate index must appear in the master layer index (no orphan rows).

Both schemas are deliberately *strict* about NaN propagation: NaN in a binary
layer is the most common silent-corruption mode (a single missing cell can
trigger Jaccard NaN that swallows the entire row from SNF), so it is rejected
up-front.

Integration
-----------
:func:`validate_dataset` is the single entry point used by the loader at the
end of :func:`amr_clonalshare.io.load_dataset` (when ``cfg.validation`` has
``schema_check: true``, the default). It iterates over every loaded
``wide_binary`` layer, applies the layer schema, and (if metadata is supplied)
applies the metadata schema with the layer index as FK target.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Optional, Sequence

import pandas as pd

# pandera>=0.18: prefer the `pandera.pandas` namespace; fall back to the
# legacy top-level import on older 0.17.x for users who pin the older floor.
try:  # pragma: no cover - import branch
    from pandera.pandas import Check, Column, DataFrameSchema, Index
except ImportError:  # pragma: no cover
    from pandera import Check, Column, DataFrameSchema, Index

from pandera.errors import SchemaError, SchemaErrors

__all__ = [
    "LayerPolicy",
    "make_layer_schema",
    "make_metadata_schema",
    "validate_layer",
    "validate_metadata",
    "validate_dataset",
    "SchemaError",
    "SchemaErrors",
]


# --------------------------------------------------------------------------- #
# Per-layer policy
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class LayerPolicy:
    """Per-layer validation policy.

    Attributes
    ----------
    name
        Role name (e.g. ``"amr"``, ``"vir"``, ``"cap"``).
    binary
        If True (default), every cell must be exactly 0 or 1.
    one_hot
        If True, additionally require that every row has **at most one** 1 —
        the classical one-hot encoding used by capsule (K-locus/O-type)
        categorical layers.
    min_columns
        Minimum number of feature columns the layer must have. Default 1.
    """
    name: str
    binary: bool = True
    one_hot: bool = False
    min_columns: int = 1


# --------------------------------------------------------------------------- #
# Schema factories
# --------------------------------------------------------------------------- #
def _binary_cell_check() -> Check:
    """Check.isin([0, 1]) — every cell must be exactly 0 or 1."""
    return Check.isin([0, 1], error="cell must be 0 or 1 (binary layer)")


def make_layer_schema(policy: LayerPolicy, strain_index_name: str = "isolate"
                      ) -> DataFrameSchema:
    """Build a pandera schema for one wide_binary layer.

    The schema enforces:
      * integer dtype + values in {0, 1} on every column (if ``policy.binary``)
      * **no NaN** (NaN tolerance disabled at the column level)
      * unique strain index — duplicate isolates are a silent-corruption mode
        the Jaccard step would not catch; reject up-front
      * one-hot row constraint if ``policy.one_hot``
    """
    df_checks = []
    if policy.min_columns > 0:
        df_checks.append(
            Check(
                lambda df: df.shape[1] >= policy.min_columns,
                error=(f"layer {policy.name!r} requires at least "
                       f"{policy.min_columns} columns"),
            )
        )
    if policy.one_hot:
        df_checks.append(
            Check(
                lambda df: (df.sum(axis=1) <= 1).all(),
                error=(f"layer {policy.name!r} is declared one-hot but at least "
                       "one row has more than one 1"),
            )
        )

    return DataFrameSchema(
        columns={},  # columns are validated by the regex Column below
        checks=df_checks,
        index=Index(str, name=strain_index_name, unique=True, nullable=False),
        strict=False,
        coerce=True,
    ).add_columns(
        {
            r".*": Column(  # apply binary check to ALL columns via regex
                int,
                checks=[_binary_cell_check()] if policy.binary else None,
                nullable=False,
                coerce=True,
                regex=True,
            )
        }
    )


def make_metadata_schema(
    fk_index: Optional[pd.Index] = None,
    strain_index_name: str = "isolate",
    required_columns: Sequence[str] = (),
) -> DataFrameSchema:
    """Build a metadata schema with foreign-key (FK) check vs ``fk_index``.

    Parameters
    ----------
    fk_index
        Master strain index (typically the aligned layer index). Every
        metadata row's index must appear in this set. ``None`` disables the
        FK check (structure-only validation).
    strain_index_name
        Expected index name on the metadata frame.
    required_columns
        Metadata column names that must be present (the column dtypes are not
        constrained — metadata is heterogeneous).
    """
    df_checks = []
    if fk_index is not None:
        fk_set = set(fk_index)

        def _fk(df: pd.DataFrame) -> bool:
            return set(df.index).issubset(fk_set)

        df_checks.append(
            Check(
                _fk,
                error=("metadata.index contains isolates not present in the "
                       "layer index (FK violation)"),
            )
        )

    columns = {
        c: Column(object, nullable=True, required=True) for c in required_columns
    }

    return DataFrameSchema(
        columns=columns,
        checks=df_checks,
        index=Index(str, name=strain_index_name, unique=True, nullable=False),
        strict=False,
        coerce=True,
    )


# --------------------------------------------------------------------------- #
# Validation entry points
# --------------------------------------------------------------------------- #
def validate_layer(df: pd.DataFrame, policy: LayerPolicy,
                   strain_index_name: str = "isolate") -> pd.DataFrame:
    """Validate one binary layer in-place; returns the (possibly coerced) frame.

    Raises :class:`pandera.errors.SchemaError` on any violation.
    """
    schema = make_layer_schema(policy, strain_index_name=strain_index_name)
    # Make sure the index has a name (pandera requires it for Index checks).
    if df.index.name != strain_index_name:
        df = df.copy()
        df.index = df.index.rename(strain_index_name)
    return schema.validate(df, lazy=False)


def validate_metadata(df: pd.DataFrame,
                      fk_index: Optional[pd.Index] = None,
                      strain_index_name: str = "isolate",
                      required_columns: Sequence[str] = ()) -> pd.DataFrame:
    """Validate the metadata frame; returns it. Raises SchemaError on failure."""
    schema = make_metadata_schema(
        fk_index=fk_index,
        strain_index_name=strain_index_name,
        required_columns=required_columns,
    )
    if df.index.name != strain_index_name:
        df = df.copy()
        df.index = df.index.rename(strain_index_name)
    return schema.validate(df, lazy=False)


def validate_dataset(dataset, *, policies: Optional[Mapping[str, LayerPolicy]] = None,
                     metadata: Optional[pd.DataFrame] = None,
                     metadata_required_columns: Sequence[str] = (),
                     ) -> None:
    """Validate every wide_binary layer in a loaded :class:`io.Dataset`.

    Parameters
    ----------
    dataset
        :class:`amr_clonalshare.io.Dataset` returned by ``load_dataset``.
    policies
        Optional per-role :class:`LayerPolicy`. Roles not listed receive a
        default binary policy. Use this to mark e.g. ``cap`` as one-hot.
    metadata
        Optional metadata frame; if given, FK-validated against
        ``dataset.strain_ids``.
    metadata_required_columns
        Column names required to be present in ``metadata``.

    Raises
    ------
    pandera.errors.SchemaError
        On the first schema violation (lazy=False).
    """
    strain_name = dataset.cfg.dataset.strain_id_column
    policies = dict(policies or {})
    for role in dataset.wide_binary_roles:
        policy = policies.get(role, LayerPolicy(name=role))
        validate_layer(dataset.binary(role), policy,
                       strain_index_name=strain_name)
    if metadata is not None:
        validate_metadata(metadata, fk_index=dataset.strain_ids,
                          strain_index_name=strain_name,
                          required_columns=metadata_required_columns)


# --------------------------------------------------------------------------- #
# Numpy convenience — used by the Hypothesis property tests
# --------------------------------------------------------------------------- #
