"""Strand/allele check between AADR (Affymetrix Human Origins) and the genome (Illumina GSA).

Usage: aadr_strand.py <ho.bim> <self.bim> <out_keep_ids> <out_json>
Keeps chr:pos sites present in both, not strand-ambiguous in AADR, whose called SELF alleles fit AADR's
two alleles on the forward strand. Reports how many fit only after complementing (those are dropped).
"""

import json
import sys

import pandas as pd

COMP = str.maketrans("ACGT", "TGCA")
AMBIG = {frozenset("AT"), frozenset("CG")}


def main(ho_fp: str, self_fp: str, out_ids: str, out_json: str) -> None:
    cols = ["chr", "id", "cm", "pos", "a1", "a2"]
    ho = pd.read_csv(ho_fp, sep="\t", header=None, names=cols, dtype=str)
    me = pd.read_csv(self_fp, sep="\t", header=None, names=cols, dtype=str)
    m = ho.merge(me, on="id", suffixes=("_ho", "_me"))
    counts = {"ho_sites": len(ho), "self_sites": len(me), "shared_positions": len(m)}
    keep, fwd, rc, bad, ambig = [], 0, 0, 0, 0
    for r in m.itertuples(index=False):
        pair = {r.a1_ho, r.a2_ho} - {"0", "."}
        if frozenset(pair) in AMBIG:
            ambig += 1
            continue
        called = {r.a1_me, r.a2_me} - {"0", "."}
        if called <= pair:
            fwd += 1
            keep.append(r.id)
        elif {a.translate(COMP) for a in called} <= pair:
            rc += 1
        else:
            bad += 1
    counts.update(
        strand_ambiguous_dropped=ambig,
        forward_consistent=fwd,
        only_if_complemented=rc,
        inconsistent=bad,
    )
    pd.Series(keep).to_csv(out_ids, index=False, header=False)
    with open(out_json, "w") as fh:
        json.dump(counts, fh, indent=1)
    print(counts)


if __name__ == "__main__":
    main(*sys.argv[1:5])
