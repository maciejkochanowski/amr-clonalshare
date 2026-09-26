"""Record Slurm CPU-hours of the release campaigns and the final test count.

Usage on the cluster: python account_jobs.py OUTPUT.json PYTEST_LOG NAME=JOBID [NAME=JOBID ...]
CPU-hours are elapsed seconds times allocated CPUs, summed over array tasks.
"""
import json
import re
import subprocess
import sys
from pathlib import Path

output, log, *pairs = sys.argv[1:]
campaigns = {}
for pair in pairs:
    name, job = pair.split('=')
    text = subprocess.run(['sacct', '-X', '-n', '-P', '-j', job, '--format=JobID,State,ElapsedRaw,AllocCPUS'],
                          check=True, capture_output=True, text=True).stdout
    rows = [line.split('|') for line in text.splitlines() if line]
    states = {}
    for _, state, _, _ in rows:
        states[state] = states.get(state, 0) + 1
    seconds = sum(int(e)*int(c) for _, _, e, c in rows)
    cpus = sorted({int(c) for _, _, _, c in rows})
    campaigns[name] = dict(job_id=job, tasks=len(rows), states=states, cpu_hours=seconds/3600,
                           cpus_per_task=cpus)
passed = re.search(r'(\d+) passed', Path(log).read_text())
result = dict(campaigns=campaigns, cpu_hours_total=sum(c['cpu_hours'] for c in campaigns.values()),
              tests_passed=int(passed.group(1)) if passed else None,
              empirical_cpus_per_task=campaigns.get('mic_empirical', {}).get('cpus_per_task', [None])[0])
Path(output).write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
