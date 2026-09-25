"""Explicit MIC interval inference with strict JSON and source provenance.

Input columns lo and hi are on the log2 concentration scale; equal bounds mean
genuinely exact data, not numeric imputation of a dilution interval. Untyped
labels are rejected.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import tempfile
import numpy as np
import scipy
from .mic_inference import exact_gaussian_interval
from .mic_null_bootstrap import null_bootstrap_interval, null_bootstrap_test


def _strict(value):
    if isinstance(value, dict):
        return {str(k): _strict(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_strict(v) for v in value]
    if isinstance(value, (float, np.floating)) and not math.isfinite(value):
        return None
    return value


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, nargs='?')
    parser.add_argument('--input', dest='input_option', type=Path,
                        help='the input table, named as the other commands name theirs; '
                             'give either this or the positional argument')
    parser.add_argument('--method', choices=['exact', 'bootstrap', 'null-test'], required=True)
    parser.add_argument('--panel-edges', help='JSON array of known log2 cut points, including negative values; '
                        'with --panel-column, a JSON object mapping each panel name to its cut points')
    parser.add_argument('--panel-column', help='input column naming the panel each record was read on')
    parser.add_argument('--covariate-column', action='append', help='input column naming the level of one '
                        'categorical covariate (a laboratory, a country) whose fixed effects enter the model; '
                        'repeat for several covariates')
    parser.add_argument('--bootstrap', type=int, default=199)
    parser.add_argument('--seed', type=int)
    parser.add_argument('--rho', type=float, help='candidate population rho; null-test only')
    parser.add_argument('--alpha', type=float, default=.05)
    parser.add_argument('--tolerance', type=float, default=.002)
    parser.add_argument('--workers', type=int, default=1,
                        help='processes for the bootstrap fits (0: every CPU the process may use); results do not depend on it')
    parser.add_argument('--output', type=Path, required=True)
    try:
        args = parser.parse_args(argv)
        if (args.input is None) == (args.input_option is None):
            raise ValueError('name the input table once, either as the positional '
                             'argument or with --input')
        args.input = args.input if args.input is not None else args.input_option
        delattr(args, 'input_option')
        if args.input.resolve() == args.output.resolve() or args.output.exists():
            raise ValueError('output already exists or equals input; overwriting is not allowed')
        # The output location is checked now, not after a bootstrap of hours.
        args.output.parent.mkdir(parents=True, exist_ok=True)
        tempfile.TemporaryFile(dir=args.output.parent).close()
        with args.input.open(newline='', encoding='utf-8-sig') as handle:
            reader = csv.DictReader(handle)
            if not {'lo', 'hi', 'lineage'} <= set(reader.fieldnames or []):
                raise ValueError('input must contain lo, hi and lineage columns')
            rows = list(reader)

        def column(name):
            values = []
            for line, row in enumerate(rows, start=2):
                try:
                    values.append(float(row[name]))
                except (TypeError, ValueError):
                    raise ValueError(f'line {line}: {name} {row[name]!r} is not a number; '
                                     f'write -inf or inf for an open end') from None
            return np.array(values)
        lo = column('lo')
        hi = column('hi')
        lineage = np.array([row['lineage'] for row in rows], dtype=object)
        if args.method == 'exact':
            if not np.array_equal(lo, hi):
                raise ValueError('exact inference requires genuinely exact observations with lo == hi')
            if args.covariate_column is not None:
                raise ValueError('the exact pivot has no covariate; use --method bootstrap')
            if args.rho is not None or args.panel_edges is not None:
                raise ValueError('--rho and --panel-edges do not apply to the exact pivot')
            answer = exact_gaussian_interval(lo, lineage, alpha=args.alpha)
            reportable = True
        else:
            if args.panel_edges is None or args.seed is None:
                missing = [name for name, value in (('--panel-edges', args.panel_edges),
                                                    ('--seed', args.seed)) if value is None]
                raise ValueError('bootstrap inference requires ' + ' and '.join(missing))
            edges = json.loads(args.panel_edges)
            panel = None
            if args.panel_column is not None:
                if args.panel_column not in (reader.fieldnames or []):
                    raise ValueError(f'input has no column {args.panel_column!r}')
                panel = np.array([row[args.panel_column] for row in rows], dtype=object)
            covariate = None
            if args.covariate_column:
                for column in args.covariate_column:
                    if column not in (reader.fieldnames or []):
                        raise ValueError(f'input has no column {column!r}')
                covariate = [np.array([row[column] for row in rows], dtype=object) for column in args.covariate_column]
            common = dict(panel_edges=edges, seed=args.seed, panel=panel, covariate=covariate,
                          n_boot=args.bootstrap, alpha=args.alpha, workers=args.workers)
            if args.method == 'null-test':
                if args.rho is None:
                    raise ValueError('null-test requires --rho')
                answer = null_bootstrap_test(lo, hi, lineage, rho=args.rho, **common)
            else:
                if args.rho is not None:
                    raise ValueError('--rho applies to null-test only; the bootstrap method inverts the test')
                answer = null_bootstrap_interval(lo, hi, lineage,
                                                tolerance=args.tolerance, **common)
            reportable = answer.reportable
        package = Path(__file__).parent
        dependencies = ['mic_inference.py', '_mic_likelihood.py', '_mic_quadrature.py', '_mic_adaptive.py',
                        '_mic_piecewise.py', 'mic_null_bootstrap.py', '_mic_panels.py', '_normal_numerics.py',
                        '_mic_boundary.py', 'attribution.py', 'mic_inference_cli.py']
        manifest = dict(schema='amr_clonalshare.mic_inference/1',
            input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
            source_sha256={name: hashlib.sha256((package/name).read_bytes()).hexdigest()
                           for name in dependencies},
            python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
            arguments={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
            nonfinite_encoding='null; inspect status and failed-bootstrap count',
            reportable=reportable, result=answer.as_dict())
        text = json.dumps(_strict(manifest), indent=2, allow_nan=False) + '\n'
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(mode='w', dir=args.output.parent,
                                             prefix='.mic-', suffix='.tmp', delete=False) as handle:
                temporary = Path(handle.name)
                handle.write(text)
                handle.flush()
                os.fsync(handle.fileno())
            mask = os.umask(0)
            os.umask(mask)
            temporary.chmod(0o666 & ~mask)
            # Linking is atomic and fails when the output path already exists.
            # A file system without hard links (FAT, many network shares)
            # gets an exclusive create instead.
            try:
                os.link(temporary, args.output)
            except FileExistsError:
                raise
            except OSError:
                with open(args.output, 'x', encoding='utf-8') as out:
                    out.write(text)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
        return 0 if reportable else 3
    except (ValueError, OSError, ArithmeticError, KeyError, TypeError) as exc:
        print(f'MIC inference: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
