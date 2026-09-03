# Reproducing the reported numbers

Three files carry three different claims and they are deliberately not merged.

`pyproject.toml` states the versions the package supports. It uses ranges,
because a user installing the tool should not be forced onto one stack.

`uv.lock` and `requirements-lock.txt` state the versions that produced the
numbers in the manuscript. `requirements-lock.txt` is an exact pin list read
from the installed distribution metadata of the environment that ran them, for
readers who do not use the uv toolchain.

`Dockerfile` builds that environment from the pin list on a pinned Python base
image. A difference between a rerun inside the image and a number in the paper
is therefore a difference in the work rather than in the stack.

## The two stacks the results were computed on

Every result was produced on one of two environments and each artefact records
which, in its own provenance block.

| | interpreter | numpy | scipy | pandas |
|---|---|---|---|---|
| workstation | 3.12.13 | 2.5.2 | 1.18.1 | 3.0.5 |
| cluster | 3.11.16 | 2.4.6 | 1.17.1 | 3.0.5 |

Running the same analysis on two stacks is not redundancy. A failure that
appears on one and not on the other is an artefact of the machine, and
reporting such an artefact as a defect of the work is the most damaging error
available in an audit.

## Seeds

Every stochastic component is seeded and the seed is written into the
artefact. The estimator seed is 42 in the veterinary atlas and 20260902 in the
estimator benchmark. The permutation control carries its own seed, 20260902,
so that the real arm and the falsifying arm cannot share a draw.

## Building and running

    docker build -t amr-clonalshare:1.0.0 .
    docker run --rm amr-clonalshare:1.0.0

    # or, without a container
    python -m venv .venv && . .venv/bin/activate
    pip install -r requirements-lock.txt
    pip install --no-deps .
