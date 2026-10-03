"""Write the web fonts of the reports and the local form into src/amr_clonalshare/fonts.

    python scripts/build_fonts.py

The text is set in Latin Modern, the typeface of LaTeX, and the charts in TeX
Gyre Heros, a Helvetica-type sans-serif. Both families come with TeX Live under
the GUST Font License, which allows them to be redistributed and embedded; the
licence is written beside the fonts. Each font is cut to the characters a report
uses (Latin with its accented letters, punctuation, Greek, arithmetic and
comparison signs) and stored as WOFF2, so that a report stays one file that
opens without an internet connection. Needs fontTools with Brotli and the
fonts of a TeX Live installation (kpsewhich).
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from fontTools import subset

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "src" / "amr_clonalshare" / "fonts"
FONTS = {
    "lmroman10-regular": "lmroman10-regular.otf",
    "lmroman10-bold": "lmroman10-bold.otf",
    "lmroman10-italic": "lmroman10-italic.otf",
    "lmmono10-regular": "lmmono10-regular.otf",
    "texgyreheros-regular": "texgyreheros-regular.otf",
    "texgyreheros-bold": "texgyreheros-bold.otf",
}
UNICODES = "U+0020-007E,U+00A0-017F,U+0391-03C9,U+2010-2027,U+2030-203A,U+2070-209F,U+2190-2193,U+2212-2213,U+2215,U+221E,U+2248,U+2260-2265,U+00D7"


def locate(name: str) -> Path:
    return Path(subprocess.run(["kpsewhich", name], capture_output=True, text=True, check=True).stdout.strip())


def main() -> None:
    OUT.mkdir(exist_ok=True)
    for stem, source in FONTS.items():
        subset.main([str(locate(source)), f"--unicodes={UNICODES}", "--flavor=woff2", "--layout-features=kern,liga",
                     "--no-hinting", "--desubroutinize", f"--output-file={OUT / (stem + '.woff2')}"])
        print(stem, (OUT / (stem + ".woff2")).stat().st_size, "bytes")
    licence = locate("lmroman10-regular.otf").parents[4] / "doc" / "fonts" / "lm" / "GUST-FONT-LICENSE.TXT"
    shutil.copyfile(licence, OUT / "GUST-FONT-LICENSE.txt")


if __name__ == "__main__":
    main()
