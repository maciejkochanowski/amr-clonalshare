"""EUCAST epidemiological cut-offs for the S. suis panel.

Four of the sixteen agents have a published EUCAST (T)ECOFF. For the other
twelve the cut-off is derived from the observed MIC distribution using the
method EUCAST itself uses to set its cut-offs: iterative fitting of a
log-normal to the wild-type subpopulation, Turnidge, Kahlmeter & Kronvall 2006,
Clin Microbiol Infect 12:418-425 ("ECOFFinder" / normalised resistance
interpretation).

This stays inside the microbiological frame. The cut-off separates the
wild-type MIC distribution from everything else; it is not a clinical
breakpoint, depends on no dosing regimen, and no S/I/R category is produced.

The derivation is validated where it can be: on the four agents that have a
published EUCAST value, the fitted cut-off is compared with it. Agreement there
is the evidence that the twelve derived values are trustworthy.
"""
import csv
import html
import json
import math
import os
import re
import urllib.parse
import urllib.request

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
EUCAST_SPECIES = 508
UA_E = {"User-Agent": "Mozilla/5.0 (research)"}
COVERAGE = 0.99          # EUCAST reports the 99 % cut-off
MIN_WT_FRACTION = 0.15   # a wild-type peak smaller than this is not a peak
MIN_R2 = 0.95


def published_ecoffs():
    q = urllib.parse.urlencode({"search[method]": "mic", "search[antibiotic]": "-1",
                                "search[species]": str(EUCAST_SPECIES),
                                "search[disk_content]": "-1", "search[limit]": "500"})
    s = urllib.request.urlopen(urllib.request.Request(
        "https://mic.eucast.org/search/?" + q, headers=UA_E),
        timeout=120).read().decode("utf-8", "replace")
    s = re.sub(r"<script.*?</script>", "", s, flags=re.S)
    out = {}
    for row in re.findall(r"<tr[^>]*>(.*?)</tr>", s, re.S)[1:]:
        c = [html.unescape(re.sub("<.*?>", " ", x)).strip()
             for x in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", row, re.S)]
        if len(c) < 24 or not c[0]:
            continue
        m = re.fullmatch(r"\(?\s*([0-9.]+)\s*\)?", c[-2].strip())
        out[c[0].lower()] = {"raw": c[-2].strip(),
                             "value": float(m.group(1)) if m else None,
                             "tentative": c[-2].strip().startswith("("),
                             "n_obs": c[-3]}
    return out


def parse_mic(text):
    m = re.fullmatch(r"((?:[<>]=?|=)*)\s*([0-9.]+)", str(text).strip().replace(",", "."))
    if not m:
        return None, None
    signs = re.findall(r"[<>]=?|=", m.group(1))
    return (signs[0] if signs else "="), float(m.group(2))


def load():
    per = {}
    for r in csv.DictReader(open(os.path.join(DATA, "mic_long.csv"))):
        sign, v = parse_mic(r["measurement"])
        if v is None:
            continue
        per.setdefault(r["antibiotic"].lower(), []).append((r["genome_id"], sign, v))
    return per


def snap(values):
    """Collapse MICs onto the standard two-fold dilution series.

    The submitters report the same dilution at different precision -- 0.03 and
    0.031, 0.12 and 0.125, 0.015 and 0.016 all appear in this panel. Left alone
    they become separate bins, which widens every fitted distribution and pushes
    the 99 % cut-off one to two dilutions too high.
    """
    return np.round(np.log2(np.asarray(values, float))).astype(int)


def wild_type_mode(bins, counts):
    """Lowest local maximum carrying real mass: the wild-type peak.

    Both conditions matter. Non-strict comparison makes a flat run of singletons
    at the bottom of the panel look like a peak -- for tetracycline that put the
    "wild-type mode" at 0.06 mg/L, on a single isolate, when the actual peak is
    at 2 mg/L on 47. The mass floor rejects those tails.
    """
    floor = 0.05 * counts.sum()
    for i, b in enumerate(bins):
        if counts[i] < floor:
            continue
        up = counts[i - 1] < counts[i] if i else True
        down = counts[i + 1] < counts[i] if i + 1 < len(bins) else True
        if up and down:
            return b
    return bins[int(np.argmax(counts))]


def ecoffinder(values):
    """Fit the wild-type log-normal and return the 99 % cut-off.

    For each candidate upper endpoint of the wild-type range a three-parameter
    cumulative normal (mean, sd, subpopulation size) is fitted by least squares
    to the cumulative counts up to that endpoint, as in Turnidge et al. 2006.

    Selection is anchored on the wild-type peak. Without that anchor the search
    happily walks the endpoint past the resistant population and returns the fit
    that covers everything -- for doxycycline here, where 84 % of isolates are
    non-wild-type, that gave a cut-off seven dilutions above the published one.
    """
    log2 = snap(values)
    bins = np.array(sorted(set(log2)))
    counts = np.array([(log2 == b).sum() for b in bins], float)
    cum = np.cumsum(counts)
    n_total = counts.sum()
    wt_mode = wild_type_mode(bins, counts)

    best = None
    for k in range(2, len(bins) + 1):
        x, y = bins[:k].astype(float), cum[:k]
        if x[-1] < wt_mode or y[-1] < MIN_WT_FRACTION * n_total:
            continue

        def resid(p, x=x, y=y):
            mu, log_sd, log_n = p
            return norm.cdf(x, mu, math.exp(log_sd)) * math.exp(log_n) - y

        p0 = [float(wt_mode), math.log(0.6), math.log(max(y[-1], 1.0))]
        try:
            sol = least_squares(resid, p0, max_nfev=4000)
        except Exception:
            continue
        mu, sd, n_wt = sol.x[0], math.exp(sol.x[1]), math.exp(sol.x[2])
        # the fitted peak must stay on the wild-type population, and the fitted
        # subpopulation must not be an extrapolation far beyond what was seen
        if abs(mu - wt_mode) > 1.5 or not (0.15 < sd < 2.5):
            continue
        if n_wt > 1.35 * y[-1] or n_wt < 0.70 * y[-1]:
            continue
        ss_res = float((resid(sol.x) ** 2).sum())
        ss_tot = float(((y - y.mean()) ** 2).sum()) or 1.0
        r2 = 1 - ss_res / ss_tot
        if best is None or r2 > best[0]:
            cut = mu + norm.ppf(COVERAGE) * sd
            best = (r2, dict(endpoint=float(2.0 ** x[-1]), mu=mu, sd=sd,
                             n_wt=n_wt, r2=r2, wt_mode=float(2.0 ** wt_mode),
                             raw_cut=float(2.0 ** cut),
                             ecoff=float(2.0 ** math.ceil(cut - 1e-9)),
                             wt_fraction=n_wt / n_total))
    return best[1] if best else None


def main():
    pub = published_ecoffs()
    per = load()
    genomes = sorted({g for v in per.values() for g, _, _ in v})
    rows = []
    print(f"{'agent':<16}{'n':>5}{'derived':>9}{'wt frac':>9}{'sd':>6}{'R2':>7}"
          f"{'published':>11}{'agreement':>12}")
    for a in sorted(per):
        vals = [v for _, _, v in per[a]]
        fit = ecoffinder(vals)
        p = pub.get(a) or pub.get(a.split("/")[0])
        pv = p["value"] if p else None
        agree = ""
        if fit and pv:
            d = math.log2(fit["ecoff"] / pv)
            agree = ("exact" if abs(d) < 1e-9 else
                     f"{d:+.0f} dilution" + ("s" if abs(d) > 1 else ""))
        c_e = format(fit["ecoff"], "g") if fit else "no fit"
        c_f = format(fit["wt_fraction"], ".2f") if fit else "-"
        c_s = format(fit["sd"], ".2f") if fit else "-"
        c_r = format(fit["r2"], ".3f") if fit else "-"
        c_p = p["raw"] if p else "-"
        print(f"{a:<16}{len(vals):>5}{c_e:>9}{c_f:>9}{c_s:>6}{c_r:>7}"
              f"{c_p:>11}{agree:>12}")
        rows.append(dict(agent=a, n=len(vals), fit=fit,
                         published=p["raw"] if p else None, published_value=pv,
                         published_tentative=p["tentative"] if p else None,
                         agreement=agree))

    val = [r for r in rows if r["published_value"] and r["fit"]]
    ok = sum(1 for r in val
             if abs(math.log2(r["fit"]["ecoff"] / r["published_value"])) <= 1)
    print(f"\nvalidation: {ok}/{len(val)} published cut-offs reproduced "
          f"within one doubling dilution")

    json.dump(rows, open(os.path.join(DATA, "ecoff_derived.json"), "w"),
              indent=1, default=str)

    # ---- wild-type / non-wild-type matrix -------------------------------
    def call(sign, v, cut):
        if sign in ("<=", "="):
            return 0 if v <= cut else (None if sign == "<=" else 1)
        if sign == "<":
            return 0 if v <= cut else None
        if sign == ">":
            return 1 if v >= cut else None
        if sign == ">=":
            return 1 if v > cut else None
        return None

    keep, cell, stats = [], {}, {}
    for r in rows:
        if not r["fit"]:
            continue
        a = r["agent"]
        # a published value always wins over a derived one
        cut = r["published_value"] or r["fit"]["ecoff"]
        col = {}
        for g, sign, v in per[a]:
            c = call(sign, v, cut)
            if c is not None:
                col[g] = c
        n1 = sum(col.values())
        stats[a] = dict(cut=cut, source="published" if r["published_value"] else "derived",
                        n=len(col), nwt=n1, wt=len(col) - n1)
        if min(n1, len(col) - n1) >= 20:
            keep.append(a)
            cell[a] = col

    out = os.path.join(DATA, "wildtype_matrix.csv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["genome_id"] + keep)
        for g in genomes:
            w.writerow([g] + [cell[a].get(g, "") for a in keep])

    print(f"\n{'agent':<16}{'cut-off':>9}{'source':>11}{'non-WT':>8}{'WT':>6}  usable")
    for a in sorted(stats):
        s = stats[a]
        print(f"{a:<16}{s['cut']:>9g}{s['source']:>11}{s['nwt']:>8}{s['wt']:>6}"
              f"  {'yes' if a in keep else 'no (minority class < 20)'}")
    print(f"\n{len(keep)} agents x {len(genomes)} isolates -> {out}")


main()
