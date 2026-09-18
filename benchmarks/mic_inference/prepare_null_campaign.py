"""Generate the prespecified confirmation task table without reading outcomes."""
from pathlib import Path
import json
import numpy as np
from benchmarks.mic_inference.calibrate_null import null_design_grid


def main():
    tasks = []
    counts = {}
    for cell in null_design_grid():
        if cell['reading'] == 'exact' or cell['domain'] == 'nonidentifiable_control':
            continue
        replicates = {'core': 1000, 'distribution_stress': 200, 'off_grid_nuisance': 500}[cell['domain']]
        counts[cell['domain']] = counts.get(cell['domain'], 0) + replicates
        for start in range(0, replicates, 100):
            tasks.append((cell['cell'], start, 25, 199, 500))
    rng = np.random.default_rng(20260916074)
    rng.shuffle(tasks)
    output = Path('benchmarks/mic_inference/tasks_null_confirmation.tsv')
    if output.exists():
        raise FileExistsError(output)
    output.write_text(''.join(' '.join(map(str, row))+'\n' for row in tasks))
    print(json.dumps(dict(tasks=len(tasks), datasets=counts, total=sum(counts.values())), indent=2))


if __name__ == '__main__':
    main()
