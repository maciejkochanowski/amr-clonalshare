"""Order-independent interpretation of contributed Pathogen Detection calls."""
from collections import defaultdict


def parse_ast(field: str) -> dict[str, float]:
    """Return I/R versus S; retain a conflicting agent as a missing value.

    Repeated calls that agree after binary recoding are one observation.
    Opposing calls have no defensible precedence in the contributed field and
    are excluded for that isolate-agent pair, not for other agents. Keeping
    the key preserves the count of agents present. ND is not a measurement.
    This historical I/R grouping is an analysis convention, not a harmonized
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
