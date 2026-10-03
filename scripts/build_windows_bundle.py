#!/usr/bin/env python3
"""Build the Windows folder that runs without any installation.

    python scripts/build_windows_bundle.py [--out dist/windows] [--wheel PATH]

The folder holds the embeddable distribution of Python from python.org, the
runtime part of the pinned stack of requirements-lock.txt as Windows wheels
(the test runner and its helpers stay out), the package itself, the example
collections of the repository (examples/, listed on the Examples page of the
form, ready to run), the manual as HTML pages (manual/, served by the form
under /manual), a two-page quick start (QUICK_START.pdf), and two files to
start from: AMR-ClonalShare.bat, which starts the local form
(amr_clonalshare.gui) in the browser, and amr-clonalshare.cmd, the command
line. A reader unpacks the zip anywhere and starts it; nothing is installed,
no administrator rights are needed, and the folder is removed by deleting it.

The script runs on any operating system: pip downloads the Windows wheels
for the pinned interpreter, and python.org serves the embeddable zip, whose
SHA-256 is pinned below. It writes <out>/AMR-ClonalShare-<version>-windows/,
the zip of that folder and the zip's SHA-256 beside it.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import re
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
#: The interpreter of the bundle: the newest 3.13 with a Windows build at
#: the time of the release. The wheels of the pinned stack exist for it.
PYTHON_VERSION = "3.13.15"
PYTHON_SHA256 = "d1f04d990aee1253d8569e8e5104e30fa9f5fa830899f14843448872d936a2cf"
PYTHON_URL = "https://www.python.org/ftp/python/{v}/python-{v}-embed-amd64.zip"

#: The packages of requirements-lock.txt the folder needs to run; the test
#: runner and its helpers are left out of the interpreter's path.
RUNTIME = ("numpy", "pandas", "scipy", "pyyaml", "python-dateutil", "six", "threadpoolctl",
           "openpyxl", "et-xmlfile")
#: The test runner and the property-testing library the test fixtures import,
#: kept in a folder of their own that the interpreter does not see;
#: scripts/check_windows_bundle.py puts it on the path for the check.
CHECK_HELPERS = ("pytest", "pluggy", "iniconfig", "packaging", "pygments", "hypothesis", "sortedcontainers")

LAUNCHER = """@echo off
setlocal
title AMR-ClonalShare {version}
echo AMR-ClonalShare {version}
echo Starting the local form; the browser opens in a moment.
echo Keep this window open while you work. Close it to stop.
echo.
set "AMR_CLONALSHARE_EXAMPLES=%~dp0examples"
"%~dp0python\\python.exe" -X utf8 -m amr_clonalshare.gui --root "%USERPROFILE%\\AMR-ClonalShare" --examples "%~dp0examples" --manual "%~dp0manual"
if errorlevel 1 (
  echo.
  echo AMR-ClonalShare did not start. The message above says why.
  pause
)
"""

COMMAND = """@echo off
set "AMR_CLONALSHARE_EXAMPLES=%~dp0examples"
"%~dp0python\\python.exe" -X utf8 -m amr_clonalshare.cli %*
"""

README = """AMR-ClonalShare {version} for Windows
=====================================

Quick start, three steps:

  1. Double-click AMR-ClonalShare.bat. A window opens and stays open while
     you work; the browser opens on the Examples page.
  2. Click "Run the analysis" under the fictional teaching collection. The
     page of the run refreshes itself; after a few seconds it shows the
     results table and the report, and lists the files of the run.
  3. For your own data, go to "New run": choose or paste the tables; one
     sheet with the isolate, its lineage and one column per antimicrobial is
     enough. The next page shows the first rows of every table beside the
     column choices; check them and run. "Check the input only" is the right
     first run on a new table.

QUICK_START.pdf shows these steps with the pages of the form; the manual is
in the folder too (the "Manual" link at the foot of every page of the form).

Windows may warn before the first start that the folder is not from a known
publisher (SmartScreen): choose "More info", then "Run anyway". The folder
is not signed; its SHA-256 is beside the zip.

The Examples page also runs the Escherichia coli, Streptococcus suis and
Salmonella collections (seconds to a few minutes each) and any configuration
file already on this computer. Every run
writes its configuration, its record and its reports into a folder of its own
under AMR-ClonalShare in your user folder (C:\\Users\\<you>\\AMR-ClonalShare), and
the page of the run links to them. Close the window to stop.

Nothing is installed: the folder carries its own Python {python} and the
packages of the release, and it is removed by deleting it. Nothing leaves the
computer; the form answers this computer only.

The command line is here too. In a command prompt opened in this folder:

    amr-clonalshare.cmd example ecoli --results-dir <folder>
    amr-clonalshare.cmd run --config examples\\ssuis\\config.yaml --results-dir <folder>
    amr-clonalshare.cmd init <metadata.csv> <mic.csv> > analysis.yaml
    amr-clonalshare.cmd compare --help

The same package is on PyPI (pip install amr-clonalshare) and its manual is
at https://maciejkochanowski.github.io/amr-clonalshare/. The licence of the
package is in LICENSE; the licence of Python is in python\\LICENSE.txt.
"""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def version() -> str:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    found = re.search(r'^version\s*=\s*"([^"]+)"', text, re.M)
    assert found, "pyproject.toml names no version"
    return found.group(1)


def download_python(cache: Path) -> Path:
    cache.mkdir(parents=True, exist_ok=True)
    target = cache / f"python-{PYTHON_VERSION}-embed-amd64.zip"
    if not target.exists():
        url = PYTHON_URL.format(v=PYTHON_VERSION)
        print(f"downloading {url}")
        with urllib.request.urlopen(url, timeout=120) as response:
            target.write_bytes(response.read())
    digest = sha256(target)
    if digest != PYTHON_SHA256:
        target.unlink()
        raise SystemExit(f"the embeddable Python does not match its pinned SHA-256: {digest}")
    return target


def build_wheel(out: Path) -> Path:
    print("building the wheel")
    subprocess.run([sys.executable, "-m", "build", "--wheel", "--outdir", str(out)], cwd=ROOT, check=True)
    wheels = sorted(out.glob("amr_clonalshare-*.whl"))
    assert wheels, "no wheel was built"
    return wheels[-1]


def build_manual(cache: Path) -> "Path | None":
    """The manual as HTML pages, built with mkdocs when it is installed;
    without it the folder ships no manual and the form says where it is."""
    if shutil.which("mkdocs") is None:
        print("mkdocs is not installed: the folder ships no manual")
        return None
    target = cache / "manual-site"
    print("building the manual")
    subprocess.run(["mkdocs", "build", "--quiet", "--site-dir", str(target)], cwd=ROOT, check=True)
    return target


def build_quick_start(cache: Path) -> "Path | None":
    """QUICK_START.pdf from docs/quick_start.md, through pandoc and
    LibreOffice when both are installed."""
    source = ROOT / "docs" / "quick_start.md"
    if not source.is_file() or shutil.which("pandoc") is None or shutil.which("soffice") is None:
        print("pandoc or LibreOffice is not installed: the folder ships no quick start")
        return None
    work = cache / "quick-start"
    work.mkdir(parents=True, exist_ok=True)
    docx = work / "QUICK_START.docx"
    print("building the quick start")
    command = ["pandoc", str(source), "--resource-path", str(source.parent), "-o", str(docx)]
    reference = _plain_reference(work)
    if reference is not None:
        command += ["--reference-doc", str(reference)]
    subprocess.run(command, check=True)
    subprocess.run(["soffice", "--headless", "--convert-to", "pdf", "--outdir", str(work), str(docx)],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    pdf = work / "QUICK_START.pdf"
    return pdf if pdf.is_file() else None


def _plain_reference(work: Path) -> "Path | None":
    """pandoc's default Word template with black headings in the text font,
    when python-docx is installed to make the change."""
    try:
        import docx
        from docx.shared import RGBColor
    except ImportError:
        return None
    reference = work / "reference.docx"
    subprocess.run(["pandoc", "-o", str(reference), "--print-default-data-file", "reference.docx"], check=True)
    document = docx.Document(str(reference))
    for style in document.styles:
        if style.name.startswith(("Heading", "Title")):
            style.font.color.rgb = RGBColor(0, 0, 0)
            style.font.name = "Calibri"
    from docx.shared import Cm
    for section in document.sections:
        section.top_margin = section.bottom_margin = Cm(1.8)
        section.left_margin = section.right_margin = Cm(2.0)
    document.save(str(reference))
    return reference


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=ROOT / "dist" / "windows")
    parser.add_argument("--wheel", type=Path, default=None, help="the wheel of the package (default: built now)")
    parser.add_argument("--cache", type=Path, default=ROOT / "dist" / "cache", help="where downloads are kept")
    args = parser.parse_args(argv)
    ver = version()
    out: Path = args.out
    out.mkdir(parents=True, exist_ok=True)
    name = f"AMR-ClonalShare-{ver}-windows"
    folder = out / name
    if folder.exists():
        shutil.rmtree(folder)
    python_dir = folder / "python"
    python_dir.mkdir(parents=True)

    with zipfile.ZipFile(download_python(args.cache)) as archive:
        archive.extractall(python_dir)
    # The embeddable interpreter ignores site-packages until its ._pth file
    # says otherwise; the packages go under Lib\\site-packages.
    pth = next(python_dir.glob("python3*._pth"))
    stdlib = next(python_dir.glob("python3*.zip")).name
    pth.write_text(f"{stdlib}\n.\nLib\\site-packages\nimport site\n", encoding="ascii")

    wheel = args.wheel or build_wheel(out)
    site = python_dir / "Lib" / "site-packages"
    site.mkdir(parents=True)
    minor = ".".join(PYTHON_VERSION.split(".")[:2])
    abi = "cp" + minor.replace(".", "")
    print(f"installing the pinned stack for Windows, Python {minor}")
    lock = {line.split("==")[0].strip().lower(): line.strip()
            for line in (ROOT / "requirements-lock.txt").read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.startswith("#")}
    missing = [name for name in RUNTIME + CHECK_HELPERS if name not in lock]
    assert not missing, f"requirements-lock.txt does not pin {missing}"
    # pip would compile the modules with the interpreter that runs this
    # script; the folder's own Python cannot read those files and writes its
    # own on first use.
    for target, names in ((site, RUNTIME), (python_dir / "Lib" / "check-helpers", CHECK_HELPERS)):
        subprocess.run([sys.executable, "-m", "pip", "install", "--quiet", "--no-deps", "--no-compile",
                        "--only-binary=:all:", "--platform", "win_amd64", "--python-version", minor,
                        "--implementation", "cp", "--abi", abi, "--abi", "none", "--target", str(target),
                        *[lock[n] for n in names]], check=True)
        shutil.rmtree(target / "bin", ignore_errors=True)
    subprocess.run([sys.executable, "-m", "pip", "install", "--quiet", "--no-deps", "--no-compile", "--target",
                    str(site), str(wheel)], check=True)
    # pip writes launcher scripts for a Unix layout under bin/; the folder
    # starts its programmes through python.exe -m, so they only mislead.
    shutil.rmtree(site / "bin", ignore_errors=True)
    # The test suites of NumPy, SciPy and pandas are a third of their size
    # and nothing of the folder runs them.
    for tests in sorted(p for p in site.rglob("*") if p.is_dir() and p.name in ("tests", "test")):
        if tests.exists() and site in tests.parents and "amr_clonalshare" not in tests.parts:
            shutil.rmtree(tests, ignore_errors=True)

    (folder / "AMR-ClonalShare.bat").write_text(LAUNCHER.format(version=ver), encoding="ascii", newline="\r\n")
    (folder / "amr-clonalshare.cmd").write_text(COMMAND, encoding="ascii", newline="\r\n")
    (folder / "README.txt").write_text(README.format(version=ver, python=PYTHON_VERSION), encoding="utf-8",
                                       newline="\r\n")
    shutil.copyfile(ROOT / "LICENSE", folder / "LICENSE")
    shutil.copytree(ROOT / "examples", folder / "examples",
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    manual = build_manual(args.cache)
    if manual is not None:
        shutil.copytree(manual, folder / "manual", ignore=shutil.ignore_patterns("*.map"))
    quick_start = build_quick_start(args.cache)
    if quick_start is not None:
        shutil.copyfile(quick_start, folder / "QUICK_START.pdf")

    archive_path = out / f"{name}.zip"
    if archive_path.exists():
        archive_path.unlink()
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(folder.rglob("*")):
            if path.is_file():
                archive.write(path, f"{name}/{path.relative_to(folder).as_posix()}")
    archive_path.write_bytes(buffer.getvalue())
    digest = sha256(archive_path)
    (out / f"{name}.zip.sha256").write_text(f"{digest}  {archive_path.name}\n", encoding="ascii")
    size = archive_path.stat().st_size / 2**20
    print(f"{archive_path} ({size:.1f} MiB), sha256 {digest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
