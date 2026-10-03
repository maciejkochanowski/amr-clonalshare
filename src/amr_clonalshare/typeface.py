"""The typefaces of the reports and the local form, embedded so that a page
looks the same on every computer and needs no internet connection.

Text is set in Latin Modern, the typeface of LaTeX, and charts in TeX Gyre
Heros, a Helvetica-type sans-serif; both are under the GUST Font License
(``fonts/GUST-FONT-LICENSE.txt``) and are cut to the characters a report uses
by ``scripts/build_fonts.py``.
"""
from __future__ import annotations

import base64
from functools import lru_cache
from pathlib import Path

__all__ = ["font_faces", "TEXT", "MONO", "CHART"]

TEXT = '"Latin Modern Roman","CMU Serif",Cambria,Georgia,serif'
MONO = '"Latin Modern Mono","CMU Typewriter Text",Consolas,ui-monospace,monospace'
CHART = '"TeX Gyre Heros","Helvetica Neue",Helvetica,Arial,"Liberation Sans",sans-serif'

_FACES = (("Latin Modern Roman", "lmroman10-regular", 400, "normal"),
          ("Latin Modern Roman", "lmroman10-bold", 700, "normal"),
          ("Latin Modern Roman", "lmroman10-italic", 400, "italic"),
          ("Latin Modern Mono", "lmmono10-regular", 400, "normal"),
          ("TeX Gyre Heros", "texgyreheros-regular", 400, "normal"),
          ("TeX Gyre Heros", "texgyreheros-bold", 700, "normal"))


@lru_cache(maxsize=1)
def font_faces() -> str:
    """@font-face rules with every font as a data URI; empty when the fonts
    are not beside the module, and the pages then fall back to the system's."""
    folder = Path(__file__).resolve().parent / "fonts"
    rules = []
    for family, stem, weight, style in _FACES:
        path = folder / f"{stem}.woff2"
        if not path.is_file():
            return ""
        data = base64.b64encode(path.read_bytes()).decode("ascii")
        rules.append(f'@font-face{{font-family:"{family}";src:url(data:font/woff2;base64,{data}) format("woff2");'
                     f'font-weight:{weight};font-style:{style};font-display:block}}')
    return "\n".join(rules)
