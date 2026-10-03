"""One name per antimicrobial across the tables of a run.

A call table and an MIC table of the same collection often name the agents
differently: a laboratory system writes the WHONET code (``CIP``, with the
test method appended as ``CIP_NM`` for an MIC or ``CIP_ND5`` for a disk), a
spreadsheet the abbreviation of a breakpoint table, a report the full name.
The run keeps one name per agent so that the two tables meet, and records
every name it changed. A name it does not know is kept as written.
"""
from __future__ import annotations

import re
from typing import Dict, Iterable, Tuple

__all__ = ["canonical_agent", "unify_agents"]

#: Codes and abbreviations, lower case, to the full name the run uses. The
#: three-letter codes are those of WHONET and of the EUCAST and CLSI tables;
#: a code is only mapped as a whole token, so ``st`` in a column name is not
#: streptomycin.
_SYNONYMS: Dict[str, str] = {
    "amk": "amikacin", "amp": "ampicillin", "amx": "amoxicillin", "amo": "amoxicillin",
    "amc": "amoxicillin-clavulanic acid", "aug": "amoxicillin-clavulanic acid",
    "amoxicillin/clavulanic acid": "amoxicillin-clavulanic acid",
    "amoxicillin-clavulanate": "amoxicillin-clavulanic acid",
    "amoxicillin/clavulanate": "amoxicillin-clavulanic acid",
    "co-amoxiclav": "amoxicillin-clavulanic acid",
    "apr": "apramycin", "azm": "azithromycin", "atm": "aztreonam",
    "caz": "ceftazidime", "ctx": "cefotaxime", "cro": "ceftriaxone", "fep": "cefepime",
    "fox": "cefoxitin", "cef": "cephalothin", "cfq": "cefquinome", "xnl": "ceftiofur", "tio": "ceftiofur",
    "chl": "chloramphenicol", "cip": "ciprofloxacin", "cli": "clindamycin", "col": "colistin", "cst": "colistin",
    "dox": "doxycycline", "enr": "enrofloxacin", "ery": "erythromycin",
    "ffn": "florfenicol", "flr": "florfenicol", "fof": "fosfomycin", "gen": "gentamicin", "gm": "gentamicin",
    "ipm": "imipenem", "imp": "imipenem", "kan": "kanamycin", "lin": "lincomycin", "lnz": "linezolid",
    "mar": "marbofloxacin", "mem": "meropenem", "mer": "meropenem", "nal": "nalidixic acid", "neo": "neomycin",
    "nit": "nitrofurantoin", "oxa": "oxacillin", "pen": "penicillin", "pen g": "penicillin",
    "rif": "rifampicin", "smx": "sulfamethoxazole", "spt": "spectinomycin", "spe": "spectinomycin",
    "str": "streptomycin", "sxt": "trimethoprim-sulfamethoxazole",
    "trimethoprim/sulfamethoxazole": "trimethoprim-sulfamethoxazole",
    "trimethoprim/sulfonamide": "trimethoprim-sulfamethoxazole",
    "co-trimoxazole": "trimethoprim-sulfamethoxazole", "cotrimoxazole": "trimethoprim-sulfamethoxazole",
    "tet": "tetracycline", "tcy": "tetracycline", "tgc": "tigecycline", "tia": "tiamulin",
    "til": "tilmicosin", "tmp": "trimethoprim", "tul": "tulathromycin", "tyl": "tylosin",
    "van": "vancomycin",
}

#: What a laboratory system appends to a code: the WHONET test-method suffix
#: (``_NM`` an MIC, ``_NE`` an Etest, ``_ND10`` a 10 µg disk) and a unit.
_SUFFIX = re.compile(r"[\s_]+(?:n[med]\d*|mic|disk|disc|zone)$|\s*\((?:mg/l|µg/ml|ug/ml)\)$", re.IGNORECASE)


def canonical_agent(name: str) -> str:
    """The full name of an agent written as a code, else the name as given,
    with surrounding spaces removed."""
    text = str(name).strip()
    stripped = _SUFFIX.sub("", text).strip()
    key = re.sub(r"\s+", " ", stripped.lower())
    full = _SYNONYMS.get(key)
    if full is None and stripped != text:
        full = _SYNONYMS.get(text.lower())
    return full if full is not None else text


def unify_agents(names: Iterable[str]) -> Tuple[Dict[str, str], Dict[str, str]]:
    """``(mapping, changed)``: every distinct name to the name the run uses,
    and the subset that differs from what was written."""
    mapping = {str(n): canonical_agent(n) for n in dict.fromkeys(names)}
    changed = {k: v for k, v in mapping.items() if k != v}
    return mapping, changed
