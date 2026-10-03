#!/usr/bin/env python3
"""Run the commands of one array task of a queue file on the cores of a node.

    python benchmarks/campaign/run_queue.py QUEUE TASK [--cores N]

Every line of QUEUE is ``task<TAB>cores<TAB>command``. The lines of TASK are
started in file order whenever enough cores are free, each in its own shell
with the environment of the batch script and its core count in QUEUE_CORES,
so that one task keeps the node busy until its last command ends however
unequal the commands are. The task fails
if any of its commands fails; the command, its exit status and its time are
printed as each ends.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("queue", type=Path)
    parser.add_argument("task", type=int)
    parser.add_argument("--cores", type=int, default=len(os.sched_getaffinity(0)))
    args = parser.parse_args()
    items = []
    for line in args.queue.read_text(encoding="utf-8").splitlines():
        task, cores, command = line.split("\t", 2)
        if int(task) == args.task:
            items.append((int(cores), command))
    if not items:
        print(f"no commands for task {args.task} in {args.queue}", file=sys.stderr)
        return 2
    if max(cores for cores, _ in items) > args.cores:
        print("a command needs more cores than the task has", file=sys.stderr)
        return 2
    free, running, failed = args.cores, {}, 0
    pending = list(items)
    while pending or running:
        # start what fits, in file order; a wide command waits for its cores
        while pending and pending[0][0] <= free:
            cores, command = pending.pop(0)
            process = subprocess.Popen(["bash", "-c", command],
                                       env={**os.environ, "QUEUE_CORES": str(cores)})
            running[process] = (cores, command, time.time())
            free -= cores
        if not running:
            break
        finished = [p for p in running if p.poll() is not None]
        if not finished:
            time.sleep(1)
            continue
        for process in finished:
            cores, command, started = running.pop(process)
            free += cores
            failed += process.returncode != 0
            print(f"exit={process.returncode} cores={cores} seconds={time.time() - started:.0f} {command}",
                  flush=True)
    print(f"task {args.task}: {len(items)} commands, {failed} failed", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
