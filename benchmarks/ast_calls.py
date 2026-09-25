"""Reading of contributed NCBI Pathogen Detection fields: AST calls, serovar, year."""
from collections import defaultdict

import numpy as np

SEROVAR_MARK = " serovar "


def parse_ast(field: str) -> dict[str, float]:
    """Return I/R versus S; retain a conflicting agent as a missing value.

    Repeated calls that agree after binary recoding are one observation.
    Opposing calls have no defensible precedence in the contributed field and
    are excluded for that isolate-agent pair, not for other agents. Keeping
    the key preserves the count of agents present. ND is not a measurement.
    Grouping I with R is an analysis convention, not a harmonized
    interpretation of the EUCAST increased-exposure category.
    """
    if not isinstance(field, str):
        return {}
    calls: dict[str, set[float]] = defaultdict(set)
    for part in field.strip('"').split(','):
        if '=' not in part:
            continue
        drug, call = part.rsplit('=', 1)
        drug, call = drug.strip().lower(), call.strip().upper()
        if drug and call in ('S', 'I', 'R'):
            calls[drug].add(float(call != 'S'))
    return {drug: next(iter(values)) if len(values) == 1 else float('nan')
            for drug, values in calls.items()}


def serovar(name):
    """The serovar written after the serovar mark of an isolate name, or ''."""
    if not isinstance(name, str) or SEROVAR_MARK not in name:
        return ""
    return name.split(SEROVAR_MARK, 1)[1].strip()


def year(value):
    """The year of a collection date that starts with four digits, else NaN."""
    if not isinstance(value, str) or len(value) < 4:
        return np.nan
    head = value[:4]
    return float(head) if head.isdigit() else np.nan
