"""Report sections: AADR ancestry, pharmacogenomics (PharmCAT), HLA, polygenic scores, traits.

Each builder returns a stub when its results are missing, so report.py can run at any point.
"""

import html
import json
from pathlib import Path

import numpy as np
import pandas as pd

from report_ancestry import nearest, scatter


def fmt_af(v) -> str:
    from report import fmt_af as f

    return f(v)


ROOT = Path.cwd()  # the run directory runs/<sample>/ (see workflow/Snakefile)
RES = ROOT / "results"
esc = html.escape


def stub(title: str, what: str) -> str:
    return f'<section><h2>{esc(title)}</h2><p class="note">Not run yet: {esc(what)}</p></section>'


# ---------- AADR ----------


def aadr_regions() -> dict[str, list[str]]:
    from run_config import get
    return get("ancestry.regions")


def zoom_populations() -> list[str]:
    from run_config import get
    return get("ancestry.zoom_populations")


def load_aadr(label):
    s = pd.read_csv(RES / "aadr" / f"proj_{label}.sscore", sep="\t").rename(columns={"#FID": "FID", "#IID": "IID"})
    pcs = [c for c in s.columns if c.endswith("_AVG")]
    pops = pd.read_csv(RES / "aadr" / "id_pop.tsv", sep=r"\s+", names=["IID", "pop"])
    s = s.merge(pops, on="IID", how="left")
    s.loc[s["IID"] == "SELF", "pop"] = "SELF"
    return s, pcs


TITLE6 = "6 · Ancestry against present-day West Eurasian populations"


def aadr_section(num) -> str:
    from run_config import get
    M = RES / "membership"
    need = ("aadr_weur_pooled.json", "aadr_zoom.json", "aadr_weur.json", "aadr_likelihood.json")
    if not (RES / "aadr" / "proj_weur.sscore").exists() or not all((M / f).exists() for f in need):
        return stub(TITLE6, "scripts/aadr.sh, ancestry_membership.py and ancestry_likelihood.py").replace("<section>", '<section id="aadr">', 1)
    strand = json.loads((RES / "aadr" / "strand.json").read_text())
    w, wp = load_aadr("weur")
    b, bp = load_aadr("zoom")
    b = b[b["pop"].isin(zoom_populations()) | (b["IID"] == "SELF")]
    counts = w[w["IID"] != "SELF"]["pop"].value_counts()
    pools = get("ancestry.pooled_groups")
    region_of = {p: r for r, ps in aadr_regions().items() for p in ps}
    region_of.update({g: region_of.get(ps[0], g) for g, ps in pools.items()})
    n_snps = int(w.loc[w["IID"] == "SELF", "ALLELE_CT"].iloc[0]) // 2
    out_fp = RES / "aadr" / "zoom_outliers.tsv"
    n_out = sum(1 for l in out_fp.read_text().splitlines() if l.strip()) if out_fp.exists() else 0
    pop_short = lambda s: s.replace("_", " ")
    pooled = json.loads((M / "aadr_weur_pooled.json").read_text())
    zoom = json.loads((M / "aadr_zoom.json").read_text())
    per_pop = json.loads((M / "aadr_weur.json").read_text())
    lik = json.loads((M / "aadr_likelihood.json").read_text())
    top_region = next(iter(pooled["self_region_posterior"]))
    lda_rows = "".join(
        f"<tr><td>{esc(pop_short(r['pop']))}</td><td>{esc(region_of.get(r['pop'], ''))}</td><td>{r['n']}</td>"
        f"<td>{r['posterior']:.2f}</td><td>{r['loo_accuracy'] * 100:.0f}%</td></tr>" for r in pooled["self_top"][:6])
    lik_rows = "".join(
        f"<tr><td>{esc(pop_short(r['pop']))}</td><td>{r['n_total']}</td><td>{r['rank1_fraction'] * 100:.0f}%</td>"
        f"<td>{r['mean_delta_logL']:.0f}</td></tr>" for r in lik["ranking"][:6])
    # Leave-one-out accuracy for the populations closest to this person (before pooling).
    near = [r["pop"] for r in per_pop["self_top"][:8]]
    loo_rows = "".join(f"<tr><td>{esc(pop_short(g))}</td><td>{counts.get(g, 0)}</td><td>{per_pop['loo_accuracy'].get(g, 0) * 100:.0f}%</td></tr>"
                       for g in near)
    weak = [pop_short(g) for g in near if per_pop["loo_accuracy"].get(g, 1) < 0.6]
    top_loo = sorted(per_pop["loo_accuracy"].items(), key=lambda kv: -kv[1])[:3]
    top2 = [pop_short(r["pop"]) for r in pooled["self_top"][:2]]
    pool_txt = "; ".join(f"{esc(g)} ({', '.join(pop_short(p) for p in ps)})" for g, ps in pools.items()) or "none"
    weak_txt = (f"Among the groups closest to you, {', '.join(esc(x) for x in weak)} are recognised less than 60% of the time, so "
                f"labels finer than their region are not supported." if weak else
                "The groups closest to you are well separated, so the population-level label is supported.")
    return f"""
<section id="aadr">
  <h2>{TITLE6}</h2>
  <p class="plainline"><b>In plain words:</b> {len(counts)} present-day European, Anatolian and Caucasus populations, many more than
  1000 Genomes has. Two independent methods ask which populations your DNA resembles most. A third check asks how well the data can
  tell each population apart at all. <b>Result: your closest match is {esc(top2[0])}, then {esc(top2[1])} ({esc(top_region)}).</b></p>
  <p>Reference: Allen Ancient DNA Resource (AADR) v66, present-day individuals genotyped on the Affymetrix Human Origins array
  ({num(int(counts.sum()))} people; shotgun-sequenced genomes left out so everyone was measured the same way). Only positions typed on
  both arrays were used, each checked for strand: {num(strand["forward_consistent"])} kept, {num(strand["only_if_complemented"])} strand-flipped,
  {num(strand["strand_ambiguous_dropped"])} strand-ambiguous. {num(n_snps)} SNPs after LD pruning. A population missing from this reference
  shows up as its nearest available neighbour, so read a population name as "most similar of these", not as an origin.</p>

  <h3>Method 1: membership probability on the genetic map</h3>
  <p class="note">You are placed on the PCA map <i>together with</i> the reference people, not projected onto it afterwards, because projection
  pulls a single sample toward the centre. Only the {pooled["n_pcs_used"]} components that rise above random noise are used
  (Marchenko–Pastur limit). A linear discriminant analysis then gives the probability that you belong to each group, if you belong to
  one of them; such probabilities tend to be overconfident, so read them as a ranking. Pooled groups: {pool_txt}. The last column is the
  check: how often that group's own members are assigned back to it when left out.</p>
  <div class="tablewrap"><table class="num"><thead><tr><th>Group</th><th>Region</th><th>n</th><th>Your membership probability</th>
  <th>Own members recognised</th></tr></thead><tbody>{lda_rows}</tbody></table></div>
  <p class="note">By region: {", ".join(f"{esc(k)} {v:.2f}" for k, v in list(pooled["self_region_posterior"].items())[:3])}.
  Fit test for your best group: p = {pooled["fit_p_best"]:.2f} (a small p would mean "none of these").</p>

  <h3>Method 2: how likely your genotypes are under each group's allele frequencies</h3>
  <p class="note">No map involved. Every group is cut down to {lik["subsample_n"]} random people, {lik["replicates"]} times, so larger groups don't win
  simply by having better-estimated frequencies. "Ranked first" is how often a group came out top.</p>
  <div class="tablewrap"><table class="num"><thead><tr><th>Group</th><th>n</th><th>Ranked first</th><th>Mean log-likelihood gap to best</th>
  </tr></thead><tbody>{lik_rows}</tbody></table></div>

  <h3>How far can labels go? Leave-one-out accuracy</h3>
  <p class="note">For the populations closest to you (before pooling): the share of each group's own members the method recognises when they
  are left out. The most distinct groups in the reference score {", ".join(f"{esc(pop_short(g))} {v * 100:.0f}%" for g, v in top_loo)}.
  {weak_txt}</p>
  <div class="tablewrap"><table class="num"><thead><tr><th>Group</th><th>n</th><th>Own members recognised</th></tr></thead>
  <tbody>{loo_rows}</tbody></table></div>

  <h3>The map (for orientation only)</h3>
  <p class="note">West Eurasian PCA with you projected onto it. Distances on this picture are distorted (see above), so read it as a picture,
  not a measurement.</p>
  {scatter(w, wp, "pop", pop_short)}
  <details><summary>Regional zoom: {esc(get("ancestry.zoom_name"))}</summary>
  <p class="note">PCA recomputed on these groups only. {n_out} outlying reference individuals (more than 6 SD on PC1–4) removed; axes show the
  central 98%. {zoom["n_pcs_used"]} component(s) rise above noise here.</p>
  {scatter(b, bp, "pop", pop_short, robust=True)}</details>
</section>"""


# ---------- HLA ----------


def hla_section() -> str:
    fp = RES / "hla" / "hla_calls.tsv"
    if not fp.exists():
        return stub("8 · HLA types", "scripts/hla.R").replace(
            "<section>", '<section id="hla">', 1
        )
    h = pd.read_csv(fp, sep="\t", dtype={"allele1": str, "allele2": str})
    rows = "".join(
        f"<tr><td>HLA-{esc(r.locus)}</td><td>*{esc(r.allele1)} / *{esc(r.allele2)}</td><td>{r.prob:.2f}</td>"
        f"<td>{'' if r.prob >= 0.8 else '⚠ uncertain: not used for any interpretation'}</td><td>{r.snps_found}/{r.model_snps}</td></tr>"
        for r in h.itertuples()
    )
    d = h.set_index("locus")
    dqb = {d.at["DQB1", "allele1"], d.at["DQB1", "allele2"]}
    dqa = {d.at["DQA1", "allele1"], d.at["DQA1", "allele2"]}
    has_dqb02 = any(a.startswith("02") for a in dqb)
    dq8 = "03:02" in dqb and any(a.startswith("03") for a in dqa)
    dq2_5 = has_dqb02 and any(a.startswith("05") for a in dqa)
    dq2_2 = has_dqb02 and "02:01" in dqa
    b_alleles = {d.at["B", "allele1"], d.at["B", "allele2"]}
    b27 = any(a.startswith("27") for a in b_alleles)
    coeliac = (
        "neither DQ2 (DQ2.5 or DQ2.2) nor DQ8"
        if not (dq2_5 or dq2_2 or dq8)
        else ", ".join(
            n for n, v in (("DQ2.5", dq2_5), ("DQ2.2", dq2_2), ("DQ8", dq8)) if v
        )
    )
    return f"""
<section id="hla">
  <h2>8 · HLA types</h2>
  <p class="plainline"><b>In plain words:</b> HLA genes build the "ID badges" your immune system uses to tell self from foreign.
  They are the most variable genes in the genome. A few specific HLA types matter medically: some predict severe reactions to
  particular drugs, and some are required for certain autoimmune conditions to develop.</p>
  <p>Imputed with HIBAG {esc("1.48.2")} from your raw chip genotypes (not the imputed ones), using the published model for the Illumina
  Global Screening Array v2.0, European ancestry, two-field resolution. "Posterior" is the model's confidence in the pair.</p>
  <div class="tablewrap"><table class="num"><thead><tr><th>Gene</th><th>Your alleles</th><th>Posterior</th><th></th>
  <th>Model SNPs found</th></tr></thead><tbody>{rows}</tbody></table></div>
  <ul class="flags">
    <li><b>Coeliac disease:</b> you carry {coeliac}. ARUP Laboratories' interpretation table rates a result "Negative for DQ2 and
    DQ8" as <i>"Not at risk"</i>; in plain terms, coeliac disease is <b>very unlikely</b>. This rests on imputed DQB1 types
    (posterior {d.at["DQB1", "prob"]:.2f}), not a clinical HLA test. (DQ2.5 = DQA1*05 with DQB1*02; DQ2.2 = DQA1*02 with DQB1*02;
    DQ8 = DQB1*03:02, per ARUP.)</li>
    <li><b>HLA-B27</b> (ankylosing spondylitis): {"present" if b27 else "not present"}. ARUP: "80-90% of AS patients are HLA-B27 positive,
    compared to 5-10% of the general population"; a negative result "does not rule out AS".</li>
    <li><b>Drug hypersensitivity HLA types</b> (HLA-A and -B): passed to PharmCAT, which reports whether the high-risk alleles are
    <i>predicted</i>. An imputed HLA type is not accepted in place of a clinical HLA test before prescribing.</li>
  </ul>
</section>"""


# ---------- PGS ----------


def pgs_hist(sc, width=360, height=120):
    counts, bins = np.array(sc["ref_hist"]), np.array(sc["ref_bins"])
    left, right, top, bottom = 8, 8, 8, 20
    pw, ph = width - left - right, height - top - bottom
    xmin, xmax = min(bins[0], sc["self_score"]), max(bins[-1], sc["self_score"])
    sx = lambda v: left + pw * (v - xmin) / (xmax - xmin)
    bw = pw * (bins[1] - bins[0]) / (xmax - xmin)
    cmax = counts.max()
    parts = [f'<svg viewBox="0 0 {width} {height}" class="chart mini" role="img">']
    for c, x0 in zip(counts, bins[:-1]):
        h = ph * c / cmax
        parts.append(
            f'<rect x="{sx(x0) + 1:.1f}" y="{top + ph - h:.1f}" width="{max(bw - 2, 1):.1f}" height="{h:.1f}" rx="2" class="histbar"/>'
        )
    mx = sx(sc["self_score"])
    parts.append(
        f'<line x1="{mx:.1f}" x2="{mx:.1f}" y1="{top - 4}" y2="{top + ph}" class="selfline"/>'
    )
    parts.append(
        f'<g class="hit" tabindex="0" data-tip="You: {sc["percentile"]:.0f}th percentile (z = {sc["z"]:+.2f})">'
        f'<rect x="{mx - 10:.1f}" y="0" width="20" height="{height}" fill="transparent"/></g>'
    )
    parts.append(
        f'<text x="{min(mx + 4, width - 60):.1f}" y="{top + 8}" class="selflabel small">You</text>'
    )
    parts.append(
        f'<line x1="{left}" x2="{width - right}" y1="{top + ph}" y2="{top + ph}" class="axis"/>'
    )
    parts.append(
        f'<text x="{left}" y="{height - 4}" class="tick">lower score</text>'
        f'<text x="{width - right}" y="{height - 4}" class="tick" text-anchor="end">higher score</text>'
    )
    parts.append("</svg>")
    return "".join(parts)


def pgs_section() -> str:
    fp = RES / "pgs" / "pgs_summary.json"
    if not fp.exists():
        return stub(
            "9 · Polygenic scores", "scripts/pgs.py (needs imputation)"
        ).replace("<section>", '<section id="pgs">', 1)
    s = json.loads(fp.read_text())
    ordinal = lambda n: f"{n}{'th' if 10 <= n % 100 <= 20 else {1: 'st', 2: 'nd', 3: 'rd'}.get(n % 10, 'th')}"
    cards = []
    for pid, sc in s["scores"].items():
        p = round(sc["percentile"])
        lo, hi = (round(x) for x in sc.get("percentile_ci95", [p, p]))
        bp = sc.get("percentile_by_population", {})
        pop_rng = f"{round(min(bp.values()))}–{round(max(bp.values()))}" if bp else "–"
        seq = sc.get("percentile_vs_sequenced_reference")
        trunc = " Note: this score lost over 20% of its weight to missing variants, so it is not the published score." \
            if sc["weight_fraction_used"] < 0.8 else ""
        cards.append(f"""<div class="pgscard"><div class="tl">{esc(sc["trait"])}</div>
  <div class="tv">{ordinal(p)}<span class="unit">percentile</span></div>
  <div class="ts">95% interval {lo}–{hi} · range across the five European populations {pop_rng} · z = {sc["z"]:+.2f}</div>
  {pgs_hist(sc)}
  <div class="ts">{esc(pid)} · {esc(sc["publication"])} · {sc["n_used"]:,} of {sc["n_variants_catalog"]:,} variants used
  ({sc["weight_fraction_used"] * 100:.0f}% of total weight){f" · vs sequenced reference: {round(seq)}" if seq is not None else ""}.{trunc}</div></div>""")
    mult = s.get("multiplicity_note", {})
    return f"""
<section id="pgs">
  <h2>9 · Polygenic scores</h2>
  <p class="plainline"><b>In plain words:</b> a polygenic score adds up thousands to millions of tiny genetic nudges into one number.
  The percentile says where you rank against a reference group: 70th means your score is higher than 70% of the reference people's.
  <b>It is relative, not a risk.</b> Someone at the 90th percentile for a condition that affects few people still most likely
  won't get it, and lifestyle, age and family history also matter.</p>
  <aside class="caveat" role="note"><b>Read these as a set.</b> With {mult.get("n_scores", 5)} scores, the highest one is expected near the
  {ordinal(round(mult.get("expected_max_percentile_under_null", 83)))} percentile by chance alone, so one high-looking score among five is not a signal
  by itself. <b>How the comparison was made:</b> you and the {s["scores"][next(iter(s["scores"]))]["ref_n"]} 1000 Genomes Europeans were cut down
  to your chip's positions and imputed together, identically, so both sides carry the same imputation effects. Your chip-artefact calls were
  removed first. The reference has only five European populations (British, Finnish, Iberian, Tuscan, Utah residents of northern/western European ancestry); the range across them shows how much the choice of comparison group moves your percentile. The scores were developed and validated mostly in people of European ancestry. The PGS
  Catalog lists their effect sizes without a consistent unit, so per-SD odds ratios and absolute risks are not given here.</aside>
  <div class="pgsgrid">{"".join(cards)}</div>
</section>"""


# ---------- Traits ----------


def traits_section() -> str:
    fp = RES / "traits.json"
    if not fp.exists():
        return stub("10 · Traits & blood group", "scripts/traits.py")
    t = json.loads(fp.read_text())
    rows = "".join(
        f"<tr><td>{esc(r['trait'])}</td><td>{esc(r['gene'])} {esc(r['rsid'])}</td><td>{esc(r['genotype'])}</td>"
        f"<td>{r['copies'] if r['copies'] is not None else '–'} × {esc(r['effect_allele'])}</td><td>{esc(r['meaning'])}</td>"
        f"<td class='src'>{esc(r['source'])}</td></tr>"
        for r in t["single"]
    )
    tas = t["tas2r38"]
    abo = t["abo"]
    abo_txt = (
        f"Blood group <b>{esc(abo['group'])}</b> ({esc(abo['alleles'])}): one B allele from rs8176746 = {esc(abo['rs8176746'])} (measured);"
        f" O-allele deletion rs8176719 imputed (dosage {abo['rs8176719']['ds']:.2f}, genotype probability {max(map(float, abo['rs8176719']['gp'].split(','))):.2f}) → {abo['o_copies']} O allele(s)."
        f" No A2 allele (rs1053878 = {esc(abo['rs1053878'])}). Rh (+/−) can't be read from SNPs."
        if "group" in abo
        else "ABO needs the imputed O-allele deletion (rs8176719); not available yet."
    )
    act = t.get("actn3")
    act_letters = ""
    if act:
        alle = [act["ref"], act["alt"]]
        act_letters = "".join(
            sorted(alle[int(c)] for c in act["gt"].replace("|", "/").split("/"))
        )
    act_txt = (
        f"ACTN3 R577X (rs1815739, imputed with genotype probability {max(map(float, act['gp'].split(','))):.2f}): "
        f"{esc(act_letters)}. T is the X (stop) allele, so you carry {act_letters.count('T')} X allele(s)"
        f"{' (R/R, both working copies)' if act_letters.count('T') == 0 else ''}. ClinVar classes the X allele as benign ('Sprinting performance')."
        if act
        else "ACTN3 (rs1815739) needs imputation; not available yet."
    )
    return f"""
<section id="traits">
  <h2>10 · Traits &amp; blood group</h2>
  <p class="plainline"><b>In plain words:</b> a few traits depend strongly on one or two DNA letters. These are mostly for curiosity.
  Alpha-1 antitrypsin is the one with health relevance: ClinVar links the Z and S alleles to alpha-1 antitrypsin deficiency and
  chronic obstructive pulmonary disease.</p>
  <div class="tablewrap"><table><thead><tr><th>Trait</th><th>SNP</th><th>Genotype</th><th>Effect allele</th><th>Effect allele means</th>
  <th>Source</th></tr></thead><tbody>{rows}</tbody></table></div>
  <ul class="flags">
    <li><b>Bitter taste (TAS2R38):</b> {esc(tas["call"])}. You are heterozygous at all three sites
    ({", ".join(f"{k} {v}" for k, v in tas["genotypes"].items())}). PAV = rs713598 G, rs1726866 G, rs10246939 C (plus strand); the chip can't phase
    the three, so PAV/AVI assumes the two common haplotypes.</li>
    <li><b>ABO:</b> {abo_txt}</li>
    <li><b>Muscle fibre type:</b> {act_txt}</li>
  </ul>
</section>"""


# ---------- PharmCAT ----------

NORMAL = {"Normal Metabolizer", "Normal Function", "Normal function", "Normal"}


def pgx_evidence() -> dict:
    """chr:pos (GRCh38) -> (imputed?, max GP) from the lifted VCF that fed PharmCAT."""
    import subprocess

    fp = ROOT / "data" / "pgx" / "pgx_hg38.vcf.gz"
    if not fp.exists():
        return {}
    q = subprocess.run(
        [
            str(ROOT / "env" / "bin" / "bcftools"),
            "query",
            "-f",
            "%CHROM:%POS\t%INFO/IMP\t[%GP]\n",
            str(fp),
        ],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    ev = {}
    for line in q.strip().splitlines():
        key, imp, gp = line.split("\t")
        ev[key] = (imp == "1", max(map(float, gp.split(","))) if gp != "." else None)
    return ev


def cyp2d6_research() -> str:
    fp = RES / "pharmcat" / "research" / "cyp2d6.report.tsv"
    if not fp.exists():
        return ""
    rows = [
        l.split("\t") for l in fp.read_text().splitlines()[2:] if l.startswith("CYP2D6")
    ]
    if not rows:
        return ""
    r = rows[0]
    return (
        f'<p class="note"><b>CYP2D6, research grade only:</b> PharmCAT\'s research mode, which ignores copy number, reads '
        f"<b>{esc(r[1])}</b> ({esc(r[2])}). A deletion (*5) or duplication would change this, and neither is visible here, "
        f"so no drug recommendation uses it.</p>"
    )


_PHENO_CACHE: dict = {}


def pharmcat_phenotype(gene: str, diplotype: str) -> str:
    """Phenotype for a diplotype from PharmCAT's own gene tables inside pharmcat.jar."""
    import zipfile
    if gene not in _PHENO_CACHE:
        try:
            with zipfile.ZipFile(ROOT / "env" / "bin" / "pharmcat.jar") as z:
                d = json.loads(z.read(f"org/pharmgkb/pharmcat/phenotype/{gene}.json"))
            _PHENO_CACHE[gene] = {e["diplotype"]: e.get("phenotype") or e.get("generesult") for e in d.get("diplotypes", [])}
        except KeyError:
            _PHENO_CACHE[gene] = {}
    tab = _PHENO_CACHE[gene]
    a, _, b = diplotype.partition("/")
    return tab.get(diplotype) or tab.get(f"{b}/{a}") or "–"


def pgx_provenance() -> dict:
    """Per gene: counts of measured / imputed / rescued / non-panel / masked sites fed to PharmCAT, and PharmCAT's
    count of missing defining positions."""
    import subprocess
    prov = {}
    reg = pd.read_csv(ROOT / "data" / "pgx" / "regions_hg19.nochr.bed", sep="\t", header=None,
                      usecols=[0, 1, 2, 3], names=["chr", "start", "end", "gene"], dtype={"chr": str})
    reg["gene"] = reg["gene"].str.replace("PX=", "")

    def gene_of37(c, p):
        hit = reg[(reg["chr"] == str(c)) & (reg["start"] < p) & (reg["end"] >= p)]
        return hit["gene"].iloc[0] if len(hit) else None

    fp = ROOT / "data" / "pgx" / "pgx_hg38.vcf.gz"
    q = subprocess.run([str(ROOT / "env" / "bin" / "bcftools"), "query", "-f", "%ID\t%INFO/IMP\t%INFO/NONPANEL\t%INFO/PX\n",
                        str(ROOT / "data" / "pgx" / "pgx.preprocessed.vcf.bgz")], capture_output=True, text=True).stdout
    rescued = set()
    rf = ROOT / "data" / "pgx" / "rescued.vcf"
    if rf.exists():
        rescued = {l.split("\t")[2] for l in rf.read_text().splitlines() if not l.startswith("#")}
    for line in q.strip().splitlines():
        rid, imp, nonp, px = (line.split("\t") + ["."] * 4)[:4]
        gene = px if px not in (".", "") else None
        if not gene:
            continue
        d = prov.setdefault(gene, {"measured": 0, "imputed": 0, "rescued": 0, "nonpanel": 0, "masked": 0, "missing": 0})
        if nonp == "1":
            d["nonpanel"] += 1
        elif imp == "1":
            d["imputed"] += 1
        else:
            d["measured"] += 1
        if rid in rescued:
            d["rescued"] += 1
    sf = RES / "pharmcat" / "site_filter.tsv"
    if sf.exists():
        s = pd.read_csv(sf, sep="\t", dtype={"chrom": str})
        for c, pos in zip(s["chrom"], s["pos_grch37"]):
            g = gene_of37(c, int(pos))
            if g:
                prov.setdefault(g, {"measured": 0, "imputed": 0, "rescued": 0, "nonpanel": 0, "masked": 0, "missing": 0})["masked"] += 1
    mf = RES / "pharmcat" / "pharmcat.match.json"
    if mf.exists():
        for g in json.loads(mf.read_text())["results"]:
            prov.setdefault(g["gene"], {"measured": 0, "imputed": 0, "rescued": 0, "nonpanel": 0, "masked": 0, "missing": 0})
            prov[g["gene"]]["missing"] = len((g.get("matchData") or {}).get("missingPositions", []))
    return prov


SKIP_PHENO = ("Indeterminate", "No Result", "Uncertain", "n/a", "ambiguous")


def pgx_section() -> str:
    from run_config import sex
    sex_label = sex()
    fp = RES / "pharmcat" / "pharmcat.report.json"
    if not fp.exists():
        return stub("7 · Pharmacogenomics (PharmCAT)", "scripts/pharmcat.sh (needs imputation)").replace("<section>", '<section id="pgx">', 1)
    r = json.loads(fp.read_text())
    prov = pgx_provenance()
    alts = {}
    af = RES / "pharmcat" / "all" / "all.match.json"
    if af.exists():
        for g in json.loads(af.read_text())["results"]:
            alts[g["gene"]] = [d["name"] for d in g.get("diplotypes", [])]
    ambiguous_genes = set()
    rows = []
    for gene, e in sorted(r["genes"].items()):
        dips = e.get("sourceDiplotypes") or []
        if not dips or not dips[0].get("allele1"):
            continue
        d = dips[0]
        name = f'{d["allele1"]["name"]}/{d["allele2"]["name"]}' if d.get("allele2") else d["allele1"]["name"]
        pheno = "; ".join(d.get("phenotypes") or []) or "–"
        if len(dips) > 1:
            opts = [f'{x["allele1"]["name"]}/{x["allele2"]["name"]}' for x in dips if x.get("allele2")]
            name = " or ".join(opts[:3]) + (f" … ({len(opts)} possible)" if len(opts) > 3 else "")
            pheno = "ambiguous: " + " / ".join(sorted({p for x in dips for p in (x.get("phenotypes") or []) if p}))
            ambiguous_genes.add(gene)
        if e.get("callSource") == "OUTSIDE":
            pheno = pheno.replace("negative", "not predicted")
            basis = "imputed HLA type (section 8), not a clinical HLA test"
        elif pheno in ("No Result", "–") and "Unknown" in name:
            basis, name = "no call: too few defining positions available", "–"
        else:
            pv = prov.get(gene, {})
            bits = [f"{pv.get(k, 0)} {lbl}" for k, lbl in (("measured", "measured"), ("imputed", "imputed (GP ≥ 0.99)"),
                                                           ("nonpanel", "measured, not in 1000G"), ("rescued", "rescued in liftover"),
                                                           ("masked", "masked"), ("missing", "missing")) if pv.get(k)]
            basis = ", ".join(bits) or "–"
        extra = ""
        others = [a for a in alts.get(gene, []) if a != name and " or " not in name]
        if others and gene not in ambiguous_genes:
            extra = "<br><span class='note'>Unphased data also fit " + ", ".join(
                f"{esc(o)} ({esc(pharmcat_phenotype(gene, o))})" for o in others[:3]) + ".</span>"
        rows.append((gene, name, pheno, basis, extra))
    gene_rows = "".join(f'<tr><td>{esc(g)}</td><td>{esc(n)}{x}</td><td>{esc(p)}</td><td class="src">{esc(b)}</td></tr>'
                        for g, n, p, b, x in rows)
    drugs, seen, held_back = [], set(), 0
    for drug, dd in sorted(r["drugs"].get("CPIC Guideline Annotation", {}).items()):
        for gl in dd.get("guidelines", []):
            for a in gl.get("annotations", []):
                ph = a.get("phenotypes") or {}
                if not ph or all(v in NORMAL or "negative" in str(v) for v in ph.values()):
                    continue
                if (a.get("population") or "general") != "general":
                    continue
                if any(any(s in str(v) for s in SKIP_PHENO) for v in ph.values()) or any(g in ambiguous_genes for g in ph):
                    held_back += 1
                    continue
                rec = (a.get("drugRecommendation") or "").split("<h4")[0].strip()
                key = (drug, rec)
                if not rec or key in seen:
                    continue
                seen.add(key)
                drugs.append((drug, ", ".join(f"{k}: {v}" for k, v in ph.items()), a.get("classification") or "", rec))
    drug_rows = "".join(f"<tr><td>{esc(d)}</td><td>{esc(p)}</td><td>{esc(c)}</td><td>{esc(rec)}</td></tr>" for d, p, c, rec in drugs)
    np_fp = RES / "pharmcat" / "nonpanel_calls.tsv"
    n_nonpanel = sum(1 for l in np_fp.read_text().splitlines()[1:] if "passed" in l and "NOT" not in l) if np_fp.exists() else 0
    return f"""
<section id="pgx">
  <h2>7 · Pharmacogenomics (PharmCAT)</h2>
  <p class="plainline"><b>In plain words:</b> some genes control how fast your body activates or clears particular medicines.
  Clinical guidelines (CPIC) turn a <i>confirmed</i> gene result into advice like "use a lower dose". <b>Nothing here is a prescription.</b>
  These results are only worth raising with a doctor if you're prescribed one of the listed drugs, and a clinic would confirm them with a
  dedicated test first. Many non-"normal" results are common: CYP3A5 *3/*3, for example, is the genotype of 89% of 1000 Genomes Europeans (448 of 503).</p>
  <p>PharmCAT {esc(r.get("pharmcatVersion", ""))} (data {esc(r.get("dataVersion", ""))}), run on your genome lifted to GRCh38. Before PharmCAT
  sees the data, four rules apply, each listed in the "Based on" column:
  (1) imputed genotypes below 99% probability are treated as missing; (2) a measured two-copy call of an allele rarer than 1% is treated as a
  chip artefact and masked; (3) {n_nonpanel} measured positions that 1000 Genomes lacks are added when they read reference (confirming an
  allele is absent); (4) positions GRCh38 renumbered are rescued rather than dropped. Genotypes are passed <i>unphased</i>, so PharmCAT
  lists every diplotype that fits. G6PD uses your measured X calls{' (a single copy in males)' if sex_label == 'male' else ''}. <b>CYP2D6 is not in this table</b>: it depends on whole-gene
  deletions and duplications that arrays can't measure.</p>
  {cyp2d6_research()}
  <div class="tablewrap"><table><thead><tr><th>Gene</th><th>Diplotype</th><th>Phenotype</th><th>Based on (sites)</th></tr></thead>
  <tbody>{gene_rows}</tbody></table></div>
  <p class="note"><b>Notes.</b> CYP2C19 *38 is the current name for the reference allele (formerly *1). UGT1A1 *28 is a repeat the chip
  can't measure; it is only ever inferred from nearby variants. Where the unphased data fit more than one diplotype, the alternatives are listed with their phenotypes. Statistical phasing from imputation supports the listed one, but only family or long-read data can prove it.</p>
  <h3>For a prescriber, after confirmatory testing: CPIC recommendations ({len(drugs)})</h3>
  <p class="note">Only recommendations from unambiguous, non-"normal" calls for the general population, de-duplicated. {held_back} further
  CPIC rows were left out because they rest on ambiguous, indeterminate or uncertain calls.</p>
  <div class="tablewrap"><table><thead><tr><th>Drug</th><th>Result used</th><th>Strength</th><th>CPIC recommendation</th></tr></thead>
  <tbody>{drug_rows}</tbody></table></div>
  <p class="note">The full PharmCAT report (all drugs; DPWG and FDA sources; missing positions) is in
  <a href="pharmcat/pharmcat.report.html">pharmcat/pharmcat.report.html</a>.</p>
</section>"""


# ---------- summary facts + imputed ClinVar block ----------


M_DIR = RES / "membership"


def summary_facts() -> dict:
    f = {}
    if (RES / "aadr" / "proj_zoom.sscore").exists():
        b, bp = load_aadr("zoom")
        mfp = RES / "membership" / "aadr_weur_pooled.json"
        if mfp.exists():
            m = json.loads(mfp.read_text())
            f["aadr_near"] = [r["pop"].replace("_", " ") for r in m["self_top"][:2]]
            per = M_DIR / "aadr_weur.json"
            if per.exists():
                pp = json.loads(per.read_text())
                weak = [r["pop"].replace("_", " ") for r in pp["self_top"][:8] if pp["loo_accuracy"].get(r["pop"], 1) < 0.6]
                f["aadr_limit"] = (f"The data can't reliably tell apart some of the closest groups ({', '.join(weak[:4])})"
                                   if weak else "The closest groups are well separated in this data")
            f["aadr_region"] = next(iter(m["self_region_posterior"]))
    fp = RES / "pharmcat" / "pharmcat.report.json"
    if fp.exists():
        r = json.loads(fp.read_text())
        notable = []
        for gene, e in sorted(r["genes"].items()):
            dips = e.get("sourceDiplotypes") or [{}]
            d = dips[0]
            ph = "; ".join(d.get("phenotypes") or [])
            if len(dips) > 1:
                notable.append(f"{gene} (ambiguous: {len(dips)} possible diplotypes)")
                continue
            if (
                ph
                and ph not in NORMAL
                and "negative" not in ph
                and ph not in ("No Result", "n/a")
                and "Indeterminate" not in ph
                and "non-responsive" not in ph
                and "Uncertain" not in ph
            ):
                notable.append(f"{gene} ({ph})")
        f["pgx_notable"] = notable
    fp = RES / "pgs" / "pgs_summary.json"
    if fp.exists():
        f["pgs"] = {
            v["trait"]: v["percentile"]
            for v in json.loads(fp.read_text())["scores"].values()
        }
    fp = RES / "traits.json"
    if fp.exists():
        f["abo"] = json.loads(fp.read_text())["abo"].get("group")
    fp = RES / "hla" / "hla_calls.tsv"
    if fp.exists():
        h = pd.read_csv(fp, sep="\t", dtype=str).set_index("locus")
        f["dq_neg"] = not any(
            a.startswith("02") or a == "03:02"
            for a in (h.at["DQB1", "allele1"], h.at["DQB1", "allele2"])
        )
    return f


def imputed_clinvar_block() -> str:
    fp = RES / "clinvar_imputed.tsv"
    if not fp.exists():
        return ""
    h = pd.read_csv(fp, sep="\t", dtype=str).rename(columns={"RS# (dbSNP)": "rs"})
    h["max_gp"] = h["max_gp"].astype(float)
    h["stars"] = h["stars"].astype(int)
    conf = h[
        (h["max_gp"] >= 0.99)
        & ((h["stars"] >= 2) | (h["category"] == "Pathogenic / likely pathogenic"))
    ]
    conf = conf.sort_values(["category", "stars"], ascending=[True, False])
    afp = RES / "clinvar_af.tsv"
    af = {}
    if afp.exists():
        a = pd.read_csv(afp, sep="\t", dtype=str)
        af = {
            (c, p): v
            for c, p, v, s in zip(a["chrom"], a["pos"], a["af_all"], a["source"])
            if s == "imputed"
        }
    rows = "".join(
        f"<tr><td>{esc(r.GeneSymbol)}</td><td>rs{esc(r.rs)}</td><td>{esc(r.gt.replace('|', '/'))}</td><td>{r.max_gp:.2f}</td>"
        f"<td>{fmt_af(af.get((str(r.chrom), str(r.pos))))}</td>"
        f"<td>{'★' * r.stars}{'☆' * (4 - r.stars)}</td><td>{esc(r.ClinicalSignificance)}</td>"
        f"<td>{esc('; '.join(p for p in dict.fromkeys(str(r.PhenotypeList).split('|')) if p not in ('not provided', 'not specified'))[:140])}</td></tr>"
        for r in conf.itertuples()
    )
    return f"""
  <h3>Positions the chip didn't measure (imputed)</h3>
  <p class="note">The same ClinVar match on imputed genotypes ({len(h)} records where you carry the variant allele). Shown: calls with genotype
  probability ≥ 0.99 (99% accurate in the masked check) that are pathogenic or have ≥ 2 review stars. Imputation can't recover rare
  variants reliably, so an empty pathogenic list here does not rule anything out. Genotype 0/1 = one copy, 1/1 = two copies of the variant allele.</p>
  <div class="tablewrap"><table><thead><tr><th>Gene</th><th>rsID</th><th>Genotype</th><th>Probability</th><th>1000G freq</th><th>Review</th>
  <th>ClinVar classification</th><th>Condition / drug</th></tr></thead><tbody>{rows}</tbody></table></div>"""
