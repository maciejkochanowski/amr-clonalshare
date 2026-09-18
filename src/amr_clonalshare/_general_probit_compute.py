"""Computation controls around the unchanged selected statistical kernel."""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
from typing import Any

import numpy as np
import scipy

from ._general_probit.profile_bootstrap_research import bootstrap as kernel


def digest(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def runtime_identity() -> dict:
    return dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                platform=sys.platform, machine=platform.machine(), byteorder=sys.byteorder,
                host=platform.node(), processor=platform.processor(),
                thread_environment={name: os.environ.get(name) for name in
                                    ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")})


@dataclass(frozen=True)
class GeneralComputeOptions:
    """Planning allowance, table cap and optional local resumable node journal.

    The memory allowance is an array/workspace estimate, not an operating-system
    RSS limit. Python, allocator and numerical-library overhead remain external.
    The preferred table allowance expands to the required minimum, or contracts
    to fit the total allowance. Process-owned journal locks release on exit.
    """
    memory_budget_mb: int = 512
    table_cache_mb: int = 128
    cache_dir: str | None = None

    def __post_init__(self):
        for name in ("memory_budget_mb", "table_cache_mb"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{name} must be a positive integer")
        if self.table_cache_mb >= self.memory_budget_mb:
            raise ValueError("table_cache_mb must be smaller than memory_budget_mb")
        if self.cache_dir is not None and (not isinstance(self.cache_dir, str) or not self.cache_dir.strip()):
            raise ValueError("cache_dir must be a nonempty path string or None")


def memory_plan(sizes, options: GeneralComputeOptions, *, B=4999) -> dict:
    """Account for retained tables, dense scoring, count encoding and workspace.

    Count/CSR allowance uses a worst-case distinct bin per group, including
    duplicate intermediates. Workspace is deliberately reserved in addition to
    retained arrays; it is an estimate rather than a certified process bound.
    """
    m = np.asarray(sizes, dtype=np.int64)
    bins = sum(int(size) + 1 for size in set(m.tolist()))
    groups = len(m)
    parameter_count = len(kernel.A_GRID) * len(kernel.ALT_RHOS)
    alternative = 17 * bins * parameter_count
    one_null = 17 * bins * len(kernel.A_GRID)
    scores = 8 * B * parameter_count
    count_encoding = 128 * B * groups + 64 * B
    maximum = int(m.max(initial=1))
    conditional_order = max(256, 16 * math.ceil(32 * math.sqrt(maximum + 1) / 16))
    normal_order = max(96, 4 * maximum)
    # Frozen numerical LRUs: at most 128 conditional bases, 64 conditional
    # rules and 32 normal rules. Bases can also hold observed unique counts.
    selected_counts = min(maximum + 1, max(128, groups))
    numerical_cache = (128 * selected_counts * conditional_order * 8
                       + 64 * conditional_order * 16 + 32 * normal_order * 16)
    workspace = 32 * 2**20 + 6 * 8 * 128 * max(normal_order, conditional_order)
    non_table = scores + count_encoding + numerical_cache + workspace
    minimum_table = alternative + one_null
    available_table = options.memory_budget_mb * 2**20 - non_table
    preferred_table = options.table_cache_mb * 2**20
    table_budget = max(minimum_table, min(preferred_table, available_table))
    total = table_budget + non_table
    return dict(geometry_bins=bins, n_groups=groups, reference_replicates=B,
                statistic_parameters=parameter_count, alternative_table_bytes=alternative,
                one_null_table_bytes=one_null, retained_table_cap_bytes=table_budget,
                preferred_table_bytes=preferred_table, minimum_table_bytes=minimum_table,
                minimum_total_bytes=minimum_table + non_table,
                suggested_memory_budget_mb=math.ceil((minimum_table + non_table) / 2**20),
                reference_score_matrix_bytes=scores, count_and_csr_allowance_bytes=count_encoding,
                numerical_cache_reserve_bytes=numerical_cache,
                workspace_reserve_bytes=workspace, planned_array_bytes=total,
                allowance_bytes=options.memory_budget_mb * 2**20,
                feasible=(alternative + one_null <= table_budget and total <= options.memory_budget_mb * 2**20),
                scope="Array planning allowance including worst-case numerical caches for this geometry; excludes prior larger-geometry global numerical caches, Python, allocator, BLAS and process overhead; not an RSS cap")


class GeneralComputeSession:
    """Reuse one exact ordered geometry; never retain several geometry engines.

    Use one session per worker. It is intentionally not thread safe. A change
    of geometry releases the previous engine before constructing another.
    """
    def __init__(self, options: GeneralComputeOptions | None = None):
        self.options = options or GeneralComputeOptions()
        self._engine: Any = None
        self._geometry: tuple = ()

    def engine(self, sizes):
        geometry = tuple(map(int, sizes))
        if geometry != self._geometry:
            plan = memory_plan(sizes, self.options)
            if not plan["feasible"]:
                raise MemoryError("Geometry exceeds the overall array planning allowance")
            self._engine = None
            self._geometry = ()
            self._engine = kernel.BootstrapEngine(geometry, memory_budget_bytes=plan["retained_table_cap_bytes"])
            self._geometry = geometry
        return self._engine

    def clear(self):
        self._engine = None
        self._geometry = ()


class NodeJournal:
    """Single atomic JSON snapshot per query, with checksummed completed nodes.

    Query identities include outcomes, complete sizes, protocol, RNG and runtime.
    Temporary writes never replace a valid snapshot until flushed successfully.
    Corrupt existing records fail explicitly; unresolved nodes are not cached.
    """
    def __init__(self, directory, identity):
        self.identity = identity
        self.key = digest(identity)
        self.directory = Path(directory)
        self.path = self.directory / (self.key + ".json")
        self.lock_path = self.directory / (self.key + ".lock")
        self.nodes: dict[str, dict] = {}
        self.reused = 0

    @contextmanager
    def locked(self):
        self.directory.mkdir(parents=True, exist_ok=True)
        # Keep this inode: replacing/deleting an advisory-lock file can create
        # two independent locks for the same query. The OS releases ownership
        # even after a hard process interruption.
        lock = self.lock_path.open("a+b")
        acquired = False
        try:
            lock.seek(0, os.SEEK_END)
            if lock.tell() == 0:
                lock.write(b"\0")
                lock.flush()
            lock.seek(0)
            try:
                if os.name == "nt":
                    import msvcrt
                    msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
                else:
                    import fcntl
                    fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except OSError as error:
                raise RuntimeError("Query cache is locked by another live process; retry after it finishes") from error
            acquired = True
            if self.path.exists():
                try:
                    stored = json.loads(self.path.read_text(encoding="utf-8"))
                    payload = stored["payload"]
                    if stored["sha256"] != digest(payload) or payload["identity"] != self.identity:
                        raise ValueError("identity or checksum mismatch")
                    self.nodes = payload["nodes"]
                    if not isinstance(self.nodes, dict):
                        raise ValueError("invalid node mapping")
                    for key, node in self.nodes.items():
                        if not isinstance(node, dict):
                            raise ValueError("invalid node record")
                        if key != float(node["rho"]).hex() or node["status"] not in ("ok", "constant_outcome", "impossible_null"):
                            raise ValueError("invalid completed node")
                        if node["B"] != self.identity["B"] or not 0 <= node["exceedances"] <= node["B"]:
                            raise ValueError("invalid reference count")
                        if node["p_value"] != (1 + node["exceedances"]) / (node["B"] + 1):
                            raise ValueError("rank mismatch")
                except (KeyError, ValueError, TypeError) as error:
                    raise RuntimeError("Query cache failed integrity verification") from error
            yield self
        finally:
            if acquired:
                lock.seek(0)
                if os.name == "nt":
                    import msvcrt
                    msvcrt.locking(lock.fileno(), msvcrt.LK_UNLCK, 1)
                else:
                    import fcntl
                    fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
            lock.close()

    def put(self, node):
        if node.get("status") not in ("ok", "constant_outcome", "impossible_null"):
            raise RuntimeError("Only supported completed nodes may enter the query cache")
        self.nodes[float(node["rho"]).hex()] = node
        payload = dict(identity=self.identity, nodes=self.nodes)
        temporary = self.path.with_suffix(".tmp")
        with temporary.open("w", encoding="utf-8", newline="\n") as stream:
            json.dump(dict(payload=payload, sha256=digest(payload)), stream, sort_keys=True, allow_nan=False)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, self.path)


class JournalEngine:
    """Import-only numerical delegation; cache hits return the original node."""
    def __init__(self, engine, journal: NodeJournal):
        self.engine = engine
        self.journal = journal
        self.input_sizes = engine.input_sizes
        self.alternative = engine.alternative

    def evaluate_node(self, counts, rho, *, B, seed, case_key):
        key = float(rho).hex()
        if key in self.journal.nodes:
            self.journal.reused += 1
            return self.journal.nodes[key]
        result = self.engine.evaluate_node(counts, rho, B=B, seed=seed, case_key=case_key)
        if result.get("status") in ("ok", "constant_outcome", "impossible_null"):
            self.journal.put(result)
        else:
            raise RuntimeError("Unexpected unresolved native node; it was not cached")
        return result
