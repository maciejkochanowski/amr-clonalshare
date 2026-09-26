"""One definition of the random stream an estimator draws from.

Drawing every antimicrobial from a single generator in column order made each
reported number depend on where its agent stood in the table: another row
order, or one more drug in the file, and the numbers of the drugs already
there moved. Each agent now draws from a fresh generator with the same
prespecified stream, so the same folds and the same permutations are applied
to every antimicrobial, and a per-agent estimate is a property of that agent's
calls, of the seed and of nothing else. Common random numbers across agents
also make the drugs comparable: what differs between two of them is their
calls, not their draws.
"""
from __future__ import annotations

import numpy as np

__all__ = ["estimator_rng"]


def estimator_rng(entropy: int, stream: int) -> np.random.Generator:
    """A fresh generator for one estimator, identical for every agent."""
    return np.random.default_rng(
        np.random.SeedSequence(entropy=int(entropy), spawn_key=(int(stream),)))
