#!/usr/bin/env python3
"""Reproducible literature search behind the novelty claims in this package.

Two claims need a search rather than an assertion. The first is that no
published work splits an observed difference in antimicrobial non-wild-type
prevalence into a lineage-composition component and a within-lineage rate
component with an interval attached - the question a surveillance laboratory
asks when it must choose between transmission control and stewardship. The
second is that the demographic and econometric estimator for exactly that split
(Kitagawa 1955, and its Blinder-Oaxaca and Fairlie descendants) has not been
carried into bacterial population data.

An absence claim is only as good as the search behind it, so the search is a
script: it records the exact query string, the hit count and every returned
record, and writes the lot to JSON beside this file.

    python docs/novelty_search.py --out docs/novelty_search_results.json

Re-running it later will not reproduce the counts, because the corpus grows.
That is the point. The recorded artefact says what the corpus held on the date
in its header; a reader who re-runs it sees what has appeared since.

Europe PMC indexes PubMed, PubMed Central, Agricola and the bioRxiv/medRxiv
preprint servers in one place. Crossref is queried separately because the
methodological half of this claim lives in economics and demography journals
that Europe PMC does not index at all - which is the whole reason the transfer
is worth checking.
"""
from __future__ import annotations

import argparse
import json
import time
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

EUROPEPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CROSSREF = "https://api.crossref.org/works"

_DECOMP = ('("Kitagawa decomposition" OR "Oaxaca-Blinder" OR "Blinder-Oaxaca" '
           'OR "Oaxaca decomposition" OR "Fairlie decomposition" '
           'OR "decomposition of a difference" OR "components of a difference")')
_AMR = ('("antimicrobial resistance" OR "antibiotic resistance" '
        'OR "non-wild-type" OR "non wild type" OR resistome '
        'OR "resistance gene" OR "susceptibility testing")')
_POP = ('(bacterial OR microbial OR pathogen OR lineage OR lineages '
        'OR "sequence type" OR "sequence types" OR MLST OR "clonal complex")')

# Each entry is (id, what the query is meant to catch, query string).
QUERIES = [
    ("Q1", "the decomposition estimator applied to resistance",
     f"{_DECOMP} AND {_AMR}"),
    ("Q2", "the decomposition estimator applied to any bacterial population",
     f"{_DECOMP} AND {_POP}"),
    ("Q3", "the biological question, asked quantitatively",
     '("clonal expansion" OR "clonal spread" OR "clonal dissemination") '
     'AND ("horizontal transfer" OR "horizontal acquisition" OR '
     '"independent acquisition" OR "de novo acquisition") '
     'AND (quantify OR quantified OR decompose OR decomposed OR '
     'attributable OR "relative contribution") '
     f'AND {_AMR}'),
    ("Q4", "composition versus rate framing of a prevalence change",
     '("composition effect" OR "compositional change" OR "structure effect" '
     'OR "within-group component" OR "between-group component" '
     'OR "mix effect" OR "shift-share") '
     'AND (prevalence OR incidence OR trend) '
     f'AND {_AMR}'),
    ("Q5", "lineage-resolved or dereplicated prevalence reporting",
     '("per-lineage" OR "one isolate per" OR dereplicated OR dereplication '
     'OR "lineage-resolved" OR "clone-corrected" OR "clone corrected") '
     'AND (prevalence OR surveillance OR reporting) '
     f'AND {_AMR}'),
    ("Q6", "tree-based phylogenetic signal for a binary resistance trait",
     '("phylogenetic signal" OR "D statistic" OR "Pagel\'s lambda" '
     'OR "Blomberg\'s K" OR homoplasy) '
     'AND ("binary trait" OR "discrete trait" OR presence-absence) '
     f'AND {_AMR}'),
    ("Q7", "the estimator anywhere in the life sciences",
     '"Kitagawa decomposition" OR "Kitagawa-Blinder-Oaxaca"'),
    ("Q9", "the estimator and the subject in one query, the decisive test",
     f'{_DECOMP} AND {_AMR} AND {_POP}'),
    ("Q10", "the two components named as such, on resistance",
     '("lineage composition" OR "clonal composition" OR "strain composition") '
     'AND ("within-lineage" OR "within-clone" OR "within-strain" OR '
     '"within-serotype") '
     f'AND {_AMR}'),
    ("Q8", "surveillance trend decomposition, any pathogen",
     '("decomposition of the change" OR "decomposing the change" '
     'OR "decomposition analysis") AND (surveillance OR monitoring) '
     f'AND {_POP}'),
]

# Crossref scores a bibliographic query by relevance and returns a hit count for
# an implicit OR of every term, so the count is meaningless as evidence: three
# earlier queries here reported 1.2, 6.3 and 1.9 million hits and none of the
# top records was on the subject. What is measurable is whether any *returned*
# record couples the estimator with microbiology, so the returned titles are
# screened and the screen is what gets recorded.
CROSSREF_QUERIES = [
    ("X1", "Kitagawa / Oaxaca-Blinder applied to resistance",
     "Kitagawa Oaxaca Blinder decomposition antimicrobial resistance"),
    ("X2", "decomposition applied to bacterial lineage composition",
     "decomposition composition effect bacterial lineage sequence type resistance"),
    ("X3", "clonal expansion versus acquisition, quantified",
     "quantifying clonal expansion versus horizontal acquisition resistance"),
]

_METHOD_TERMS = ("kitagawa", "oaxaca", "blinder", "fairlie", "shift-share",
                 "decomposition of a difference", "composition effect",
                 "components of a difference")
_SUBJECT_TERMS = ("bacter", "microb", "antimicrobial", "antibiotic",
                  "resistance", "lineage", "sequence type", "mlst",
                  "serotype", "clone", "clonal", "pathogen", "isolate")


# A third claim needs its own search, its own screen and its own terms, so it is
# kept in a separate family rather than folded into the two above: the counts
# reported for the decomposition claim must not move because a query about a
# different subject was added beside them.
#
# The claim is about an estimand, not an estimator. Inverting a noncentral
# distribution to bound its noncentrality parameter is Venables 1975, exact
# effect-size intervals built that way in the fixed-effects analysis of variance
# are Steiger 2004, and treating the group effects a design actually holds as a
# finite population with its own variance is Cornfield and Tukey 1956. None of
# that is claimed. What is searched for is whether a bacterial collection has
# been reported with both shares - the one a new draw of lineages would show and
# the one the lineages in hand do - from a single fit.
_COND = ('("finite population" OR "finite-population" OR "narrow inference" '
         'OR "broad inference" OR superpopulation OR "super-population" '
         'OR "conditional inference" OR realised OR realized)')
_SHARE = ('("intraclass correlation" OR "variance component" '
          'OR "variance components" OR "proportion of variance" '
          'OR "variance explained" OR heritability)')

ESTIMAND_QUERIES = [
    ("E1", "the conditional share of a trait variance in a bacterial population",
     f'{_COND} AND {_SHARE} AND {_POP}'),
    ("E2", "a noncentral F or noncentrality interval used on bacterial data",
     '("noncentral F" OR "non-central F" OR "noncentrality parameter" '
     'OR "noncentral chi-square" OR "non-central chi-square") '
     f'AND {_POP}'),
    ("E3", "both estimands named together, anywhere in the life sciences",
     f'{_COND} AND {_SHARE} AND (estimand OR "target of inference" '
     'OR "two questions" OR "which question")'),
    ("E4", "the small-number-of-groups limit stated for such an interval",
     f'{_SHARE} AND ("few clusters" OR "small number of clusters" '
     'OR "small number of groups" OR "number of groups is small" '
     'OR "wide confidence interval" OR "wide confidence intervals") '
     f'AND {_POP}'),
    ("E5", "a lineage share of a susceptibility trait, with an interval",
     '("clonal" OR lineage OR "sequence type" OR clone) '
     f'AND {_SHARE} AND {_AMR}'),
]

ESTIMAND_CROSSREF = [
    ("XE1", "the conditional share, applied to bacterial lineages",
     "finite population intraclass correlation bacterial lineage "
     "antimicrobial resistance"),
    ("XE2", "the two estimands reported from one fit",
     "realised versus superpopulation variance component conditional "
     "intraclass correlation noncentral F"),
]

_ESTIMAND_METHOD_TERMS = (
    "noncentral f", "non-central f", "noncentrality", "noncentral chi",
    "non-central chi", "finite population", "finite-population",
    "narrow inference", "broad inference", "superpopulation",
    "super-population", "conditional inference", "intraclass correlation",
    "variance component", "proportion of variance", "variance explained",
    "heritability")


def screen_estimand(records) -> dict:
    """How many returned titles couple a conditional-share term with microbiology.

    The subject terms are the same as above; only the method half differs,
    because the claim being tested is a different one.
    """
    both = []
    for record in records:
        title = (record.get("title") or "").lower()
        if (any(t in title for t in _ESTIMAND_METHOD_TERMS)
                and any(t in title for t in _SUBJECT_TERMS)):
            both.append(record)
    return {"n_screened": len(records), "n_method_and_subject": len(both),
            "matches": both}


# A fourth family, added on 2026-09-02 after a review found a record the first
# three could not have returned. Plumb and colleagues (Vaccine 2020;38:4273,
# doi:10.1016/j.vaccine.2020.04.048) split a rise in pneumococcal
# non-susceptibility into a serotype-redistribution part and a within-serotype
# part by direct standardisation. That is this decomposition, in a bacterial
# population, published. The first three families could not find it because
# they searched for the estimator by the names it carries in demography and
# economics, and the epidemiological literature that does the same arithmetic
# calls it standardisation, serotype replacement or clonal replacement. A
# search that can only return the vocabulary it was given is not a search, and
# the absence it produced was an artefact of the query rather than a property
# of the corpus. This family exists so the claim is tested against the words
# the field actually uses.
_EPI_METHOD = ('("direct standardisation" OR "direct standardization" '
               'OR "indirect standardisation" OR "indirect standardization" '
               'OR standardised OR standardized OR "counterfactual" '
               'OR "attributable fraction" OR "expected if" '
               'OR "holding constant" OR "adjusting for the distribution")')
_EPI_STRUCTURE = ('("serotype replacement" OR "serotype redistribution" '
                  'OR "serotype shift" OR "clonal replacement" '
                  'OR "strain replacement" OR "lineage replacement" '
                  'OR "within-serotype" OR "within serotype" '
                  'OR "within-lineage" OR "within-clone" '
                  'OR "serotype distribution")')

EPIDEMIOLOGICAL_QUERIES = [
    ("P1", "standardisation applied to a resistance change, with structure",
     f"{_EPI_METHOD} AND {_EPI_STRUCTURE} AND {_AMR}"),
    ("P2", "replacement versus within-type change, on resistance",
     f"{_EPI_STRUCTURE} AND {_AMR} AND (increase OR decrease OR trend "
     "OR change OR rise)"),
    ("P3", "the pneumococcal case, which is where the precedent was found",
     '("Streptococcus pneumoniae" OR pneumococcal OR pneumococci) '
     f'AND {_EPI_STRUCTURE} AND {_AMR}'),
    ("P4", "the same split named as components, in any pathogen",
     '("contribution of" OR "contributions of" OR "accounted for by") '
     f'AND {_EPI_STRUCTURE} AND {_AMR}'),
    ("P5", "an interval or a test attached to such a split",
     f'{_EPI_STRUCTURE} AND {_AMR} AND ("confidence interval" OR bootstrap '
     'OR "credible interval" OR "standard error")'),
]

_EPI_TERMS = ("standardis", "standardiz", "replacement", "redistribution",
              "within-serotype", "within serotype", "within-lineage",
              "within-clone", "counterfactual", "attributable fraction",
              "serotype shift", "serotype distribution", "contribution")


def screen_epidemiological(records) -> dict:
    """Returned titles that couple the split with resistance in a population.

    The screen for the first family asks for a decomposition term by its
    econometric name. This one asks for the epidemiological name, because that
    is the vocabulary the precedent is written in.
    """
    both = []
    for record in records:
        title = (record.get("title") or "").lower()
        if (any(t in title for t in _EPI_TERMS)
                and any(t in title for t in _SUBJECT_TERMS)):
            both.append(record)
    return {"n_screened": len(records), "n_method_and_subject": len(both),
            "matches": both}


def screen(records) -> dict:
    """How many returned titles couple the estimator with microbiology.

    A hit count answers "how many documents matched some term". This answers
    the question the novelty claim actually rests on, and it is the number to
    quote.
    """
    both = []
    for record in records:
        title = (record.get("title") or "").lower()
        if (any(t in title for t in _METHOD_TERMS)
                and any(t in title for t in _SUBJECT_TERMS)):
            both.append(record)
    return {"n_screened": len(records), "n_method_and_subject": len(both),
            "matches": both}


def _get(url: str, params: dict, retries: int = 3) -> dict:
    query = urllib.parse.urlencode(params)
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(f"{url}?{query}", timeout=60) as fh:
                return json.loads(fh.read().decode("utf-8"))
        except Exception as exc:                       # network, 5xx, timeout
            if attempt == retries - 1:
                return {"error": f"{type(exc).__name__}: {exc}"}
            time.sleep(2 * (attempt + 1))
    return {}


def search_europepmc(query: str, page_size: int = 50) -> dict:
    raw = _get(EUROPEPMC, {"query": query, "format": "json",
                           "pageSize": page_size, "resultType": "core"})
    if "error" in raw:
        return {"hits": None, "records": [], "error": raw["error"]}
    result = raw.get("resultList", {}).get("result", [])
    return {
        "hits": int(raw.get("hitCount", 0)),
        "records": [{"id": r.get("id"), "doi": r.get("doi"),
                     "year": r.get("pubYear"), "journal": r.get("journalTitle"),
                     "title": r.get("title")} for r in result],
    }


def search_crossref(query: str, rows: int = 50) -> dict:
    raw = _get(CROSSREF, {"query.bibliographic": query, "rows": rows})
    if "error" in raw:
        return {"hits": None, "records": [], "error": raw["error"]}
    msg = raw.get("message", {})
    return {
        "hits": int(msg.get("total-results", 0)),
        "records": [{"doi": it.get("DOI"), "year": (it.get("issued") or {})
                     .get("date-parts", [[None]])[0][0],
                     "journal": (it.get("container-title") or [None])[0],
                     "title": (it.get("title") or [None])[0]}
                    for it in msg.get("items", [])],
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).resolve().parent
                    / "novelty_search_results.json")
    args = ap.parse_args()

    out = {"searched": date.today().isoformat(),
           "systems": {"europepmc": EUROPEPMC, "crossref": CROSSREF},
           "europepmc": {}, "crossref": {}}

    for qid, intent, query in QUERIES:
        res = search_europepmc(query)
        res["intent"], res["query"] = intent, query
        res["screen"] = screen(res["records"])
        out["europepmc"][qid] = res
        print(f"{qid} {str(res['hits']):>8} hits, "
              f"{res['screen']['n_method_and_subject']} of "
              f"{res['screen']['n_screened']} returned titles carry both a "
              f"decomposition term and a microbiology term   {intent}")
        time.sleep(1)

    for qid, intent, query in CROSSREF_QUERIES:
        res = search_crossref(query)
        res["intent"], res["query"] = intent, query
        res["screen"] = screen(res["records"])
        res["hit_count_is_not_evidence"] = (
            "Crossref counts an implicit OR of every term; read the screen")
        out["crossref"][qid] = res
        print(f"{qid} {str(res['hits']):>8} hits (not evidence), "
              f"{res['screen']['n_method_and_subject']} of "
              f"{res['screen']['n_screened']} returned titles carry both   "
              f"{intent}")
        time.sleep(1)

    out["europepmc_epidemiological"] = {}
    for qid, intent, query in EPIDEMIOLOGICAL_QUERIES:
        res = search_europepmc(query)
        res["intent"], res["query"] = intent, query
        res["screen"] = screen_epidemiological(res["records"])
        out["europepmc_epidemiological"][qid] = res
        print(f"{qid} {str(res['hits']):>8} hits, "
              f"{res['screen']['n_method_and_subject']} of "
              f"{res['screen']['n_screened']} returned titles carry both a "
              f"standardisation term and a microbiology term   {intent}")
        time.sleep(1)
    epi_total = sum(r["screen"]["n_method_and_subject"]
                    for r in out["europepmc_epidemiological"].values())
    out["epidemiological_summary"] = {
        "n_queries": len(EPIDEMIOLOGICAL_QUERIES),
        "n_returned_titles_coupling_split_and_resistance": epi_total,
        "why_this_family_exists": (
            "the first family searched for the estimator by its econometric "
            "names and could not have returned a paper that does the same "
            "arithmetic under the epidemiological name; one such paper exists "
            "and is cited in NOVELTY_EVIDENCE.md"),
        "reading": ("a non-zero count here is the expected outcome, not a "
                    "failure: the claim that survives is about the packaged "
                    "estimator and its gates, not about the arithmetic"),
    }

    out["europepmc_estimand"], out["crossref_estimand"] = {}, {}
    for qid, intent, query in ESTIMAND_QUERIES:
        res = search_europepmc(query)
        res["intent"], res["query"] = intent, query
        res["screen"] = screen_estimand(res["records"])
        out["europepmc_estimand"][qid] = res
        print(f"{qid} {str(res['hits']):>8} hits, "
              f"{res['screen']['n_method_and_subject']} of "
              f"{res['screen']['n_screened']} returned titles carry both a "
              f"conditional-share term and a microbiology term   {intent}")
        time.sleep(1)

    for qid, intent, query in ESTIMAND_CROSSREF:
        res = search_crossref(query)
        res["intent"], res["query"] = intent, query
        res["screen"] = screen_estimand(res["records"])
        res["hit_count_is_not_evidence"] = (
            "Crossref counts an implicit OR of every term; read the screen")
        out["crossref_estimand"][qid] = res
        print(f"{qid} {str(res['hits']):>8} hits (not evidence), "
              f"{res['screen']['n_method_and_subject']} of "
              f"{res['screen']['n_screened']} returned titles carry both   "
              f"{intent}")
        time.sleep(1)

    estimand_total = sum(r["screen"]["n_method_and_subject"]
                         for r in list(out["europepmc_estimand"].values())
                         + list(out["crossref_estimand"].values()))
    out["estimand_summary"] = {
        "n_queries": len(ESTIMAND_QUERIES) + len(ESTIMAND_CROSSREF),
        "n_returned_titles_coupling_conditional_share_and_microbiology":
            estimand_total,
        "reading": ("a returned title carrying both terms is a record to read "
                    "and cite, not a refutation on its own: the claim is that "
                    "both shares are not reported together from one fit of a "
                    "bacterial cohort, which only the record itself can settle"),
    }

    total = sum(r["screen"]["n_method_and_subject"]
                for r in list(out["europepmc"].values())
                + list(out["crossref"].values()))
    out["summary"] = {
        "n_queries": len(QUERIES) + len(CROSSREF_QUERIES),
        "n_returned_titles_coupling_estimator_and_microbiology": total,
        "reading": ("zero is the claim: no returned record applies a "
                    "Kitagawa-family decomposition to a bacterial population. "
                    "A non-zero count is a record to read and cite, not a "
                    "failure of the search"),
    }
    print(f"\ntitles coupling the estimator with microbiology, "
          f"across every query: {total}")

    args.out.write_text(json.dumps(out, indent=2, ensure_ascii=False))
    print(f"\nwritten: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
