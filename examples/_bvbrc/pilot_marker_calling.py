"""Pilot: are BLAST-based marker calls plausible where annotation-based ones are not?

apxIVA is the internal control. It is present in essentially every
A. pleuropneumoniae isolate -- it is the species-specific diagnostic target -- so
any method that calls it at 16 % is measuring annotation coverage, not biology.
"""
import json
import os
import subprocess
import sys
import tempfile
import urllib.request

UA = {"Accept": "application/json", "User-Agent": "amr-clonalshare/1.0.0 (research)"}


def api(path, q, accept="application/json", timeout=180):
    url = f"https://www.bv-brc.org/api/{path}/?{q}&http_accept={accept}"
    h = dict(UA)
    h["Accept"] = accept
    raw = urllib.request.urlopen(
        urllib.request.Request(url, headers=h), timeout=timeout).read()
    return json.loads(raw) if accept == "application/json" else raw.decode()


def reference_protein(taxon, keyword):
    """Longest annotated exemplar of a marker, used as the BLAST query."""
    r = api("genome_feature",
            f'and(eq(taxon_id,{taxon}),keyword(%22{keyword.replace(" ", "+")}%22),'
            f'eq(feature_type,CDS))&select(feature_id,aa_sequence_md5,product,'
            f'aa_length)&limit(60)')
    r = [x for x in r if x.get("aa_length") and x.get("aa_sequence_md5")]
    if not r:
        return None, None
    best = max(r, key=lambda x: x["aa_length"])
    seq = api("protein_feature" if False else "feature_sequence",
              f'eq(md5,{best["aa_sequence_md5"]})&limit(1)')
    if not seq:
        return None, None
    return best, ">ref\n" + seq[0]["sequence"] + "\n"


def proteome(genome_id):
    return api("genome_feature",
               f"and(eq(genome_id,{genome_id}),eq(feature_type,CDS))&limit(12000)",
               accept="application/protein+fasta")


def call_markers(faa_refs, genome_ids, min_ident=80.0, min_cov=0.70):
    """Presence/absence by blastp of each reference against each proteome."""
    with tempfile.TemporaryDirectory() as td:
        qpath = os.path.join(td, "q.faa")
        lens = {}
        with open(qpath, "w") as fh:
            for name, (rec, faa) in faa_refs.items():
                seq = "".join(faa.split("\n")[1:]).strip()
                lens[name] = len(seq)
                fh.write(f">{name}\n{seq}\n")
        out = {}
        for gid in genome_ids:
            try:
                p = proteome(gid)
            except Exception:
                continue
            spath = os.path.join(td, "s.faa")
            with open(spath, "w") as fh:
                fh.write(p)
            subprocess.run(["makeblastdb", "-in", spath, "-dbtype", "prot",
                            "-out", os.path.join(td, "db")],
                           capture_output=True, check=True)
            res = subprocess.run(
                ["blastp", "-query", qpath, "-db", os.path.join(td, "db"),
                 "-outfmt", "6 qseqid pident length", "-evalue", "1e-20",
                 "-max_target_seqs", "5", "-num_threads", "2"],
                capture_output=True, text=True, check=True).stdout
            hit = {k: 0 for k in faa_refs}
            for line in res.strip().splitlines():
                q, pid, ln = line.split("\t")
                if float(pid) >= min_ident and int(ln) >= min_cov * lens[q]:
                    hit[q] = 1
            out[gid] = hit
    return out


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "app"
    if which == "app":
        taxon, n = 715, 30
        kws = {"apxIVA": "ApxIVA", "apxIIA": "ApxIIA", "apxIA": "ApxIA",
               "tbpA": "transferrin-binding protein A", "omlA": "OmlA"}
    else:
        taxon, n = 1307, 30
        kws = {"mrp": "muramidase-released", "epf": "extracellular factor",
               "sly": "suilysin", "ofs": "serum opacity factor"}

    refs = {}
    for name, kw in kws.items():
        rec, faa = reference_protein(taxon, kw)
        if faa:
            refs[name] = (rec, faa)
            print(f"ref {name:8s} {rec['aa_length']:4d} aa  {rec['product'][:46]}",
                  flush=True)
    gl = api("genome", f"and(eq(taxon_id,{taxon}),eq(genome_quality,Good))"
                       f"&select(genome_id)&limit({n})")
    gids = [g["genome_id"] for g in gl]
    print(f"\nblastp against {len(gids)} proteomes ...", flush=True)
    calls = call_markers(refs, gids)
    print(f"got calls for {len(calls)} genomes\n")
    for name in refs:
        k = sum(v[name] for v in calls.values())
        print(f"  {name:8s} {k:3d}/{len(calls)}  = {100*k/max(len(calls),1):5.1f} %")


main()
