"""Write the design list of the extension campaign, which the figures read.

    python -m benchmarks.mic_inference.write_extension_grid <output.json>
"""
import json
import sys
from pathlib import Path

from benchmarks.mic_inference.calibrate_extension import extension_design_grid


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 1:
        print(__doc__, file=sys.stderr)
        return 2
    Path(argv[0]).write_text(json.dumps(extension_design_grid(), indent=1) + '\n', encoding='utf-8')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
