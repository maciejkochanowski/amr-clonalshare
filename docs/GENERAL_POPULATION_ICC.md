# General population-model confidence limits

The general method computes limits for the actual retained group sizes, rather
than reading a fixed cut-off calibrated on a grid of designs. It is used only
when it is enabled explicitly; otherwise the population model uses the fixed
cut-off. In the validation study of `benchmarks/population_model/PROTOCOL.md`
its coverage ranged from 93.4% to 99.6% over 144 Gaussian designs with 500
fresh datasets each, and every design met the prespecified rule (upper Wilson
95% limit at least 0.95).

```yaml
population_probit:
  enabled: true
  interval_method: general
  memory_budget_mb: 512
  table_cache_mb: 128
  cache_dir: ./profile-cache   # relative to the configuration file
```

The target is a population liability ICC under independent Gaussian group
effects and conditional binomial counts. Group sizes must be fixed or
noninformative, and selection must be ignorable. The method does not change the
collection-level attribution target, its admission gates or its random streams.

With two or more repeated groups, the fitted-nuisance bootstrap uses
4,999 reference draws and an internal threshold of 0.04 for a **95% confidence
target**. This is an approximation, not an exact guarantee over fitted nuisance
parameters or the continuous parameter range. Its hull is the primary result;
all retained components are diagnostics. A finite validation panel does not
establish an exact fitted-nuisance or continuous-parameter guarantee.

With exactly one repeated group, the general method reports a standalone one-sided 95% upper confidence limit [0,U]. The construction uses 9,999 finite-Monte-Carlo reference draws, inclusive upper-tail ranks with the plus-one correction, and the prespecified grid `tuple(i/100 for i in range(101))`; it takes the next grid endpoint after the last accepted node, capped at 1. Its model-based coverage statement averages over both the data and the prespecified algorithm randomness under independent Gaussian group effects, conditional binomial counts, fixed or noninformative group sizes, and ignorable selection. It is not a two-sided interval or a guarantee conditional on a fixed deterministic seed. It must not be intersected with, or selected against, another 95% procedure without a separately justified error allocation.

Geometry and outcome states remain distinct. When at least one group is retained and every retained group is a singleton, rho is structurally unidentified and the method returns the full range [0,1] without a point estimate. With no retained data, it returns no confidence object and public limits are missing. With exactly one repeated group, a homogeneous repeated-group outcome has zero centrality and returns the unrestricted upper range [0,1], irrespective of singleton outcomes. If the one-repeat inversion accepts no raw grid node, the boundary convention reports [0,0] and retains explicit raw-empty and boundary-fill flags; this is not an empty confidence set. With at least two repeated groups, a globally constant outcome returns [0,1], whereas a genuine empty bootstrap confidence set remains empty. A failed or unresolved calculation has missing public limits and is not relabelled as a completed full, upper-only, or empty result.

## Python use and repeat computations

```python
from amr_clonalshare import (
    general_probit_icc, GeneralComputeOptions, GeneralComputeSession,
)

session = GeneralComputeSession(GeneralComputeOptions(cache_dir="./profile-cache"))
result = general_probit_icc([1, 0], [2, 2], seed=713,
                            case_key="example-trait", compute=session)
record = result.as_dict()
session.clear()
```

Choose the seed and a stable analysis identifier before examining the result.
Repeating the same query reuses completed null calculations when input counts,
sizes, seed, protocol, source and runtime identities match. A changed query has
a different cache identity. Failed nodes are never saved as completed nodes.
Cache files use atomic snapshots and checksums. A process-owned operating-system lock prevents two
workers writing simultaneously and releases automatically when its owner exits,
including after an abrupt termination. Persistent `.lock` files are normal; do
not remove or replace them while workers may run. In a configuration a relative `cache_dir` is read against the
configuration file; in Python it is read against the process working
directory, so an absolute path is preferable for scheduled jobs.

One session retains one exact ordered geometry engine. Use one session per
worker and keep its queries serial. Distinct workers can process distinct
queries, each with its own memory allowance and the same predetermined case
identifiers; worker order does not determine the random streams. No cluster
submission or remote data transfer is performed by this API.

## Resource planning

The memory allowance includes retained likelihood tables, the complete reference
score matrix, count/CSR encoding, a numerical-cache reserve and temporary array
workspace. At the fixed settings, the 4,999 by 2,091 score matrix alone uses
83,623,272 bytes. Many groups or many distinct large sizes can need substantially
more memory and table-build time. A failed preflight reports an incomplete
calculation; it does not substitute another statistical method.

Increase the overall memory allowance when necessary. The preferred table
allowance must be smaller than the overall allowance; it expands to the minimum
needed for a geometry or contracts to leave room for the reference arrays. A
resource refusal records the required bytes and suggested overall setting. This is an array planning estimate, **not a
hard operating-system memory limit**. Python, allocator, numerical-library and
process overhead, plus numerical caches left by other computations in the same
process, remain outside the estimate. Use fresh worker processes and
explicit scheduler memory limits for large jobs. Memory limits never change
the reference count, threshold or statistical grids.
