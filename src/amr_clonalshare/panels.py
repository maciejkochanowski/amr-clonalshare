"""Dilution ranges of the harmonised surveillance panels, as presets.

The tested concentrations of an agent decide which readings sit on an end
well, so they must come from the panel sheet. For the panels that the
European Union prescribes for the monitoring of resistance in food-producing
animals, the ranges are public: Commission Implementing Decision (EU)
2020/1729 of 17 November 2020, Annex, Part A, Tables 2 to 5. They are
transcribed here as twofold series from the lowest to the highest
concentration printed in the Decision. A preset is a convenience, not a
certificate: a laboratory's plate may differ from the Decision in a lot or a
product, so the record names the preset and the report asks the reader to
verify it against the plate.
"""
from __future__ import annotations

from typing import Dict, List

__all__ = ["PANEL_PRESETS", "preset_wells", "preset_names"]

_SOURCE = "Commission Implementing Decision (EU) 2020/1729, Annex, Part A"


def _series(low: float, high: float) -> List[float]:
    """Every doubling from ``low`` to ``high`` as exact powers of two: the
    Decision prints 0.015, 0.03, 0.06 and 0.12 for 1/64 to 1/8 mg/L, and a
    reading written either way is matched to its well."""
    import math
    first, last = round(math.log2(low)), round(math.log2(high))
    return [float(2.0 ** k) for k in range(first, last + 1)]


#: Each preset: the agents of the panel with their tested concentrations,
#: the source table and the organisms the panel is prescribed for.
PANEL_PRESETS: Dict[str, dict] = {
    "eu-2020-1729-salmonella-ecoli": {
        "title": "EU harmonised panel for Salmonella spp. and E. coli (first panel)",
        "source": f"{_SOURCE}, Table 2",
        "organisms": "Salmonella spp., Escherichia coli",
        "wells": {
            "amikacin": _series(4, 128), "ampicillin": _series(1, 32), "azithromycin": _series(2, 64),
            "cefotaxime": _series(0.25, 4), "ceftazidime": _series(0.25, 8), "chloramphenicol": _series(8, 64),
            "ciprofloxacin": _series(0.015, 8), "colistin": _series(1, 16), "gentamicin": _series(0.5, 16),
            "meropenem": _series(0.03, 16), "nalidixic acid": _series(4, 64), "sulfamethoxazole": _series(8, 512),
            "tetracycline": _series(2, 32), "tigecycline": _series(0.25, 8), "trimethoprim": _series(0.25, 16),
        },
    },
    "eu-2020-1729-salmonella-ecoli-second": {
        "title": "EU harmonised second panel for Salmonella spp. and E. coli (isolates resistant to cefotaxime, ceftazidime or meropenem)",
        "source": f"{_SOURCE}, Table 5",
        "organisms": "Salmonella spp., Escherichia coli",
        "wells": {
            "cefepime": _series(0.06, 32), "cefotaxime": _series(0.25, 64),
            "cefotaxime-clavulanic acid": _series(0.06, 64), "cefoxitin": _series(0.5, 64),
            "ceftazidime": _series(0.25, 128), "ceftazidime-clavulanic acid": _series(0.125, 128),
            "ertapenem": _series(0.015, 2), "imipenem": _series(0.12, 16), "meropenem": _series(0.03, 16),
            "temocillin": _series(0.5, 128),
        },
    },
    "eu-2020-1729-campylobacter": {
        "title": "EU harmonised panel for Campylobacter jejuni and C. coli",
        "source": f"{_SOURCE}, Table 3",
        "organisms": "Campylobacter jejuni, Campylobacter coli",
        "wells": {
            "chloramphenicol": _series(2, 64), "ciprofloxacin": _series(0.12, 32), "ertapenem": _series(0.125, 4),
            "erythromycin": _series(1, 512), "gentamicin": _series(0.25, 16), "tetracycline": _series(0.5, 64),
        },
    },
    "eu-2020-1729-enterococcus": {
        "title": "EU harmonised panel for Enterococcus faecalis and E. faecium",
        "source": f"{_SOURCE}, Table 4",
        "organisms": "Enterococcus faecalis, Enterococcus faecium",
        "wells": {
            "ampicillin": _series(0.5, 64), "chloramphenicol": _series(4, 128), "ciprofloxacin": _series(0.12, 16),
            "daptomycin": _series(0.25, 32), "erythromycin": _series(1, 128), "gentamicin": _series(8, 1024),
            "linezolid": _series(0.5, 64), "quinupristin-dalfopristin": _series(0.5, 64),
            "teicoplanin": _series(0.5, 64), "tetracycline": _series(1, 128), "tigecycline": _series(0.03, 4),
            "vancomycin": _series(1, 128),
        },
    },
}


def preset_names() -> List[str]:
    return list(PANEL_PRESETS)


def preset_wells(name: str) -> Dict[str, List[float]]:
    """The tested concentrations per agent of a preset, by its name."""
    if name not in PANEL_PRESETS:
        raise KeyError(name)
    return {agent: list(wells) for agent, wells in PANEL_PRESETS[name]["wells"].items()}
