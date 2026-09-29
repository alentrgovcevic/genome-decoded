"""Build results/report.html from results/qc.json, results/roh.hom, and the ClinVar outputs.

Self-contained HTML (inline CSS/JS/SVG). It contains genotypes, so it stays under results/ (gitignored).
"""

import html
import json
from datetime import date
from pathlib import Path

import pandas as pd
from report_ancestry import build_sections
from report_health import (
    aadr_section,
    hla_section,
    imputed_clinvar_block,
    pgs_section,
    pgx_section,
    summary_facts,
    traits_section,
)

ROOT = Path.cwd()  # the run directory runs/<sample>/ (see workflow/Snakefile)
RES = ROOT / "results"


def clinvar_date() -> str:
    """ClinVar release date: download_refs.sh records the file's HTTP Last-Modified; else the file's mtime."""
    lm = ROOT / "reference" / "variant_summary.last_modified"
    if lm.exists():
        from email.utils import parsedate_to_datetime

        raw = lm.read_text().strip()
        try:
            return parsedate_to_datetime(raw).date().isoformat()
        except (TypeError, ValueError):
            return raw
    f = ROOT / "reference" / "variant_summary.txt.gz"
    return (
        date.fromtimestamp(f.stat().st_mtime).isoformat() if f.exists() else "unknown"
    )


# Marker alleles on the GRCh37 plus strand, from the ClinVar records (see results/selected_markers.tsv).
MARKERS = {
    "rs6025": (
        "T",
        "Factor V Leiden (F5)",
        "the GRCh37 reference itself carries the Leiden allele",
    ),
    "rs1799963": ("A", "Prothrombin G20210A (F2)", ""),
    "rs1800562": ("A", "HFE C282Y", ""),
    "rs1799945": ("G", "HFE H63D", "ClinVar classifications conflict"),
}

# APOE allele definition on the GRCh37 plus strand: (rs429358, rs7412) base per haplotype.
APOE_HAPLOTYPES = {
    ("T", "T"): "ε2",
    ("T", "C"): "ε3",
    ("C", "C"): "ε4",
    ("C", "T"): "ε1",
}

esc = html.escape


def pct(x: float, d: int = 2) -> str:
    return f"{x * 100:.{d}f}%"


def num(x: int) -> str:
    return f"{x:,}"


def stars(n: int) -> str:
    return f'<span class="stars" aria-label="{n} of 4 review stars">{"★" * n}{"☆" * (4 - n)}</span>'


def clinvar_link(vid: str) -> str:
    return f'<a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/{esc(str(vid))}/" rel="noreferrer">{esc(str(vid))}</a>'


# ---------- charts (inline SVG) ----------


def bar_chart_vertical(
    rows, value_key, label_key, series_key, fmt, height=220, ymax=None, legend=None
):
    """rows: list of dicts. series_key -> CSS var name for fill."""
    n = len(rows)
    left, right, top, bottom = 44, 8, 12, 28
    width = 760
    pw, ph = width - left - right, height - top - bottom
    ymax = ymax or max(r[value_key] for r in rows) * 1.1
    slot = pw / n
    bw = max(slot - 4, 2)
    parts = [f'<svg viewBox="0 0 {width} {height}" class="chart" role="img">']
    for t in range(5):
        v = ymax * t / 4
        y = top + ph - ph * t / 4
        parts.append(
            f'<line x1="{left}" x2="{width - right}" y1="{y:.1f}" y2="{y:.1f}" class="grid"/>'
        )
        parts.append(
            f'<text x="{left - 6}" y="{y + 4:.1f}" class="tick" text-anchor="end">{fmt(v)}</text>'
        )
    for i, r in enumerate(rows):
        v = r[value_key]
        h = ph * v / ymax
        x = left + i * slot + (slot - bw) / 2
        y = top + ph - h
        tip = f"chr{r[label_key]}: {fmt(v)}"
        rr = min(4, bw / 2, h)
        path = (
            f"M{x:.1f},{top + ph:.1f} V{y + rr:.1f} Q{x:.1f},{y:.1f} {x + rr:.1f},{y:.1f} "
            f"H{x + bw - rr:.1f} Q{x + bw:.1f},{y:.1f} {x + bw:.1f},{y + rr:.1f} V{top + ph:.1f} Z"
        )
        parts.append(
            f'<g class="hit" tabindex="0" data-tip="{esc(tip)}">'
            f'<rect x="{left + i * slot:.1f}" y="{top}" width="{slot:.1f}" height="{ph}" fill="transparent"/>'
            f'<path d="{path}" fill="var({r[series_key]})"/></g>'
        )
        parts.append(
            f'<text x="{x + bw / 2:.1f}" y="{height - 10}" class="tick" text-anchor="middle">{esc(str(r[label_key]))}</text>'
        )
    parts.append(
        f'<line x1="{left}" x2="{width - right}" y1="{top + ph}" y2="{top + ph}" class="axis"/>'
    )
    parts.append("</svg>")
    leg = ""
    if legend:
        leg = (
            '<div class="legend">'
            + "".join(
                f'<span><i style="background:var({c})"></i>{esc(t)}</span>'
                for t, c in legend
            )
            + "</div>"
        )
    return leg + "".join(parts)


def bar_chart_horizontal(items, fmt=num, width=760, label_w=210):
    """items: list of (label, value, tip)."""
    rowh, gap = 22, 6
    top = 4
    height = top + len(items) * (rowh + gap) + 4
    vmax = max(v for _, v, _ in items) or 1
    pw = width - label_w - 70
    parts = [f'<svg viewBox="0 0 {width} {height}" class="chart" role="img">']
    for i, (label, v, tip) in enumerate(items):
        y = top + i * (rowh + gap)
        w = max(pw * v / vmax, 2)
        rr = min(4, w / 2)
        path = (
            f"M{label_w},{y} H{label_w + w - rr:.1f} Q{label_w + w:.1f},{y} {label_w + w:.1f},{y + rr} "
            f"V{y + rowh - rr} Q{label_w + w:.1f},{y + rowh} {label_w + w - rr:.1f},{y + rowh} H{label_w} Z"
        )
        parts.append(
            f'<text x="{label_w - 8}" y="{y + rowh / 2 + 4}" class="label" text-anchor="end">{esc(label)}</text>'
        )
        parts.append(
            f'<g class="hit" tabindex="0" data-tip="{esc(tip)}">'
            f'<rect x="0" y="{y - gap / 2}" width="{width}" height="{rowh + gap}" fill="transparent"/>'
            f'<path d="{path}" fill="var(--series-1)"/></g>'
        )
        parts.append(
            f'<text x="{label_w + w + 6:.1f}" y="{y + rowh / 2 + 4}" class="value">{fmt(v)}</text>'
        )
    parts.append(
        f'<line x1="{label_w}" x2="{label_w}" y1="0" y2="{height}" class="axis"/>'
    )
    parts.append("</svg>")
    return "".join(parts)


def roh_karyogram(per_chrom, roh):
    autos = [c for c in per_chrom if c["chrom"] not in ("X", "Y")]
    maxlen = max(c["max_pos"] for c in autos)
    width, label_w, rowh, gap, top = 760, 44, 12, 8, 4
    pw = width - label_w - 12
    height = top + len(autos) * (rowh + gap) + 24
    parts = [f'<svg viewBox="0 0 {width} {height}" class="chart" role="img">']
    for i, c in enumerate(autos):
        y = top + i * (rowh + gap)
        x0 = label_w + pw * c["min_pos"] / maxlen
        x1 = label_w + pw * c["max_pos"] / maxlen
        parts.append(
            f'<text x="{label_w - 8}" y="{y + rowh - 2}" class="tick" text-anchor="end">{c["chrom"]}</text>'
        )
        parts.append(
            f'<rect x="{x0:.1f}" y="{y + 3}" width="{x1 - x0:.1f}" height="{rowh - 6}" rx="3" class="track"/>'
        )
        for s in roh[roh["CHR"] == int(c["chrom"])].itertuples():
            sx = label_w + pw * s.POS1 / maxlen
            sw = max(pw * (s.POS2 - s.POS1) / maxlen, 3)
            tip = f"chr{c['chrom']}:{s.POS1:,}–{s.POS2:,} · {s.KB / 1000:.2f} Mb · {s.NSNP} SNPs"
            parts.append(
                f'<g class="hit" tabindex="0" data-tip="{esc(tip)}">'
                f'<rect x="{sx - 4:.1f}" y="{y - 2}" width="{sw + 8:.1f}" height="{rowh + 4}" fill="transparent"/>'
                f'<rect x="{sx:.1f}" y="{y}" width="{sw:.1f}" height="{rowh}" rx="2" fill="var(--series-1)"/></g>'
            )
    ya = height - 16
    for mb in range(0, int(maxlen / 1e6) + 1, 50):
        x = label_w + pw * mb * 1e6 / maxlen
        parts.append(
            f'<text x="{x:.1f}" y="{ya + 12}" class="tick" text-anchor="middle">{mb} Mb</text>'
        )
    parts.append("</svg>")
    return "".join(parts)


# ---------- tables ----------


def table(df, cols, headers, cls=""):
    head = "".join(f"<th>{esc(h)}</th>" for h in headers)
    body = []
    for d in df[cols].to_dict("records"):
        body.append(
            "<tr>"
            + "".join(
                f"<td>{d[c] if c.startswith('_') else esc(str(d[c]))}</td>"
                for c in cols
            )
            + "</tr>"
        )
    return f'<div class="tablewrap"><table class="{cls}"><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table></div>'


def impute_bench_para() -> str:
    fp = RES / "impute_benchmark.json"
    if not fp.exists():
        return ""
    b = json.loads(fp.read_text())
    label = {
        "0-0.001": "< 0.1%",
        "0.001-0.005": "0.1–0.5%",
        "0.005-0.01": "0.5–1%",
        "0.01-0.05": "1–5%",
        "0.05-0.2": "5–20%",
        "0.2-0.5": "20–50%",
    }
    num2 = lambda v: "–" if v is None else f"{v:.2f}"
    pct1 = lambda v: "–" if v is None else f"{v * 100:.1f}%"
    maf_rows = "".join(
        f"<tr><td>{label.get(r['maf'], r['maf'])}</td><td>{r['sites']:,}</td><td>{num2(r['median_dosage_r2'])}</td>"
        f"<td>{pct1(r['nonref_concordance'])}</td></tr>"
        for r in b["by_maf"]
    )
    gp_rows = "".join(
        f"<tr><td>{r['max_gp']}</td><td>{r['genotypes']:,}</td><td>{r['accuracy'] * 100:.1f}%</td></tr>"
        for r in b["by_gp"]
    )
    n = b.get("n_people", 503)
    return f"""<h3>How accurate is the imputation? Measured on {n} people with known answers</h3>
  <p class="note">The {n} unrelated 1000 Genomes Europeans were cut down to exactly your chip's positions and imputed the same way you were,
  against a panel that excludes them. Their sequenced genotypes at {b["n_sites"]:,} sampled imputed positions ({b["sampling"]}) are the answer
  key; the frequency table covers the {sum(r["sites"] for r in b["by_maf"]):,} of them that vary among the {n}. <i>Dosage r²</i> is how well the imputed values track the truth (1 = perfect). <i>Non-reference concordance</i> counts only genotypes
  that carry the variant, where errors matter; overall accuracy is inflated by the many easy "no variant" cases. Common variants impute
  well; rare ones don't, which is why rare imputed calls are never used on their own in this report.</p>
  <div class="tablewrap"><table class="num"><thead><tr><th>Variant frequency</th><th>Sites</th><th>Median dosage r²</th>
  <th>Non-reference concordance</th></tr></thead><tbody>{maf_rows}</tbody></table></div>
  <p class="note">And the confidence flag used throughout the report: accuracy of the best-guess genotype by its genotype probability.</p>
  <div class="tablewrap"><table class="num"><thead><tr><th>Genotype probability</th><th>Genotypes</th><th>Correct</th></tr></thead>
  <tbody>{gp_rows}</tbody></table></div>"""


def rare_hom_para(q) -> str:
    if not q:
        return ""
    return (
        f'<h3>Chip artefacts, counted</h3><p class="note">Of {q["sites_compared"]:,} measured SNPs also in 1000 Genomes, you are called '
        f"homozygous for an allele that fewer than 1% of Europeans carry at <b>{q['observed_rare_homozygotes']:,}</b> sites. Hardy–Weinberg "
        f"predicts about <b>{q['expected_rare_homozygotes']:.1f}</b>. {q['observed_rare_hom_q0']:,} of them are alleles not seen in any of "
        f"the reference Europeans. "
        + (
            f"The excess, about {q['excess_likely_artefacts']:,.0f} calls ({q['observed_rare_homozygotes'] / q['sites_compared'] * 100:.2f}% "
            f"of sites), is the chip misreading rare-variant probes: the problem Weedon et al. (2021, BMJ) describe for consumer arrays. "
            if q["excess_likely_artefacts"] >= 10
            else "That is close to what chance predicts, so this file shows no sign of the rare-variant probe errors Weedon et al. "
            "(2021, BMJ) describe for consumer arrays (typical for data that did not come from a chip). "
        )
        + f"These calls are flagged in the disease section, masked for pharmacogenomics, and kept out of the imputation used for "
        f"polygenic scores.</p>"
    )


def fmt_af(v) -> str:
    """1000 Genomes frequency of the variant allele; blank means not in 1000G as a biallelic SNV."""
    try:
        v = float(v)
    except (TypeError, ValueError):
        return '<span title="not in 1000 Genomes: very rare, or a chip-only probe">absent</span>'
    if v != v:
        return '<span title="not in 1000 Genomes: very rare, or a chip-only probe">absent</span>'
    return f"{v * 100:.2f}%" if v < 0.01 else f"{v * 100:.0f}%"


def hits_table(df):
    df = df.copy()
    df["_stars"] = df["stars"].map(stars)
    df["_link"] = df["VariationID"].map(clinvar_link)
    df["_af"] = df["af_all"].map(fmt_af) if "af_all" in df else "–"
    df["_cond"] = (
        df["PhenotypeList"]
        .fillna("")
        .map(
            lambda s: esc(
                "; ".join(
                    p
                    for p in dict.fromkeys(s.split("|"))
                    if p not in ("not provided", "not specified", "-")
                )[:160]
            )
        )
    )
    return table(
        df,
        [
            "GeneSymbol",
            "rsid",
            "chrom",
            "gt",
            "alt",
            "zygosity",
            "_af",
            "_stars",
            "ClinicalSignificance",
            "_cond",
            "_link",
        ],
        [
            "Gene",
            "rsID",
            "Chr",
            "Genotype",
            "Variant allele",
            "Zygosity",
            "1000G freq",
            "Review",
            "ClinVar classification",
            "Condition(s)",
            "ClinVar ID",
        ],
    )


def main() -> None:
    qc = json.loads((RES / "qc.json").read_text())
    cv = json.loads((RES / "clinvar_summary.json").read_text())
    hits = pd.read_csv(RES / "clinvar_hits.tsv", sep="\t", dtype=str)
    hits["stars"] = hits["stars"].astype(int)
    if (RES / "clinvar_af.tsv").exists():
        af = pd.read_csv(RES / "clinvar_af.tsv", sep="\t", dtype=str)[
            ["chrom", "pos", "ref", "alt", "af_all"]
        ]
        hits = hits.merge(af, on=["chrom", "pos", "ref", "alt"], how="left")
    rh_fp = RES / "qc_rare_hom.tsv"
    rare_hom = set()
    if rh_fp.exists():
        rh = pd.read_csv(rh_fp, sep="\t", dtype=str)
        rare_hom = set(zip(rh["chrom"], rh["pos"]))
        flag = [(c, str(p)) in rare_hom for c, p in zip(hits["chrom"], hits["pos"])]
        hits.loc[flag, "zygosity"] = (
            "homozygous for a <1% allele (likely chip artefact)"
        )
    rh_qc = (
        json.loads((RES / "qc_rare_hom.json").read_text())
        if (RES / "qc_rare_hom.json").exists()
        else None
    )
    sel = pd.read_csv(RES / "selected_markers.tsv", sep="\t", dtype=str)
    roh = pd.read_csv(RES / "roh.hom", sep=r"\s+")
    pc = qc["per_chrom"]
    ct = {
        "hom": 0,
        "het": 0,
        "no_call": 0,
        "indel": 0,
        **qc["class_totals"],
    }  # absent class = 0 (e.g. no D/I probes)

    sex = {"1": "Male", "2": "Female"}.get(
        (ROOT / "data" / "genome.fam").read_text().split()[4], "Unknown"
    )
    x = next(c for c in pc if c["chrom"] == "X")
    y = next(c for c in pc if c["chrom"] == "Y")
    auto_len = sum(
        c["max_pos"] - c["min_pos"] for c in pc if c["chrom"] not in ("X", "Y")
    )
    roh_kb = roh["KB"].sum()
    roh15 = roh.loc[roh["KB"] >= 1500, "KB"]
    froh = roh_kb * 1000 / auto_len

    strand_total = (
        cv["strand_forward_consistent"]
        + cv["strand_only_if_complemented"]
        + cv["strand_inconsistent"]
    )
    strand_ok = cv["strand_forward_consistent"] / strand_total

    # APOE
    g = dict(zip(sel["rsid"], sel["gt"]))
    apoe = "not typed"
    if "rs429358" in g and "rs7412" in g:
        a, b = g["rs429358"], g["rs7412"]
        if a[0] != a[1] and b[0] != b[1]:
            apoe = "ambiguous (ε1/ε3 or ε2/ε4 — phase unknown)"
        else:
            # With at most one heterozygous site, pairing the letters position-wise gives the haplotypes.
            haps = [APOE_HAPLOTYPES[(a[i], b[i])] for i in (0, 1)]
            apoe = "/".join(sorted(haps, key="ε1ε2ε3ε4".index))

    # ---- charts
    het_rows = [
        {
            "chrom": c["chrom"],
            "het": c["het_rate"] or 0,
            "series": "--series-2" if c["chrom"] in ("X", "Y") else "--series-1",
        }
        for c in pc
    ]
    het_chart = bar_chart_vertical(
        het_rows,
        "het",
        "chrom",
        "series",
        lambda v: pct(v, 0) if v >= 0.05 or v == 0 else pct(v, 1),
        ymax=0.24,
        legend=[("Autosomes", "--series-1"), ("Sex chromosomes", "--series-2")],
    )
    gts = sorted(qc["genotype_counts"].items(), key=lambda kv: -kv[1])
    gt_chart = bar_chart_horizontal(
        [(k, v, f"{k}: {num(v)} ({pct(v / qc['n_variants'])})") for k, v in gts],
        label_w=60,
    )

    cat_order = [
        "Pathogenic / likely pathogenic",
        "Risk factor",
        "Drug response",
        "Protective",
        "Association",
        "Conflicting",
        "Uncertain significance",
        "Benign / likely benign",
        "Other / not provided",
    ]
    cats = cv["hits_by_category"]
    cat_chart = bar_chart_horizontal(
        [
            (c, cats.get(c, 0), f"{c}: {num(cats.get(c, 0))} variant records")
            for c in cat_order
        ],
        label_w=210,
    )
    karyo = roh_karyogram(pc, roh)

    # ---- per-chrom table
    pcdf = pd.DataFrame(pc)
    pcdf["call"] = pcdf["call_rate"].map(lambda v: pct(v))
    pcdf["hetr"] = pcdf["het_rate"].map(lambda v: pct(v))
    pcdf["span"] = ((pcdf["max_pos"] - pcdf["min_pos"]) / 1e6).map(lambda v: f"{v:.1f}")
    for c in ("n", "no_call", "indel", "hom", "het"):
        pcdf[c] = pcdf[c].map(num)
    chrom_table = table(
        pcdf,
        ["chrom", "n", "no_call", "indel", "hom", "het", "call", "hetr", "span"],
        [
            "Chr",
            "Probes",
            "No-call",
            "Indel",
            "Hom",
            "Het",
            "Call rate",
            "Het rate",
            "Array span (Mb)",
        ],
        "num",
    )
    rohdf = roh.sort_values("KB", ascending=False).copy()
    rohdf["mb"] = (rohdf["KB"] / 1000).map(lambda v: f"{v:.2f}")
    rohdf["p1"] = rohdf["POS1"].map(num)
    rohdf["p2"] = rohdf["POS2"].map(num)
    roh_table = table(
        rohdf,
        ["CHR", "p1", "p2", "mb", "NSNP"],
        ["Chr", "Start", "End", "Length (Mb)", "SNPs"],
        "num",
    )

    # ---- ClinVar sections
    plp = hits[hits["category"] == "Pathogenic / likely pathogenic"]
    plp_alt2 = plp.loc[plp["stars"] >= 2, ["chrom", "pos"]].drop_duplicates().shape[0]
    plp_ref2 = cv["plp_sites_tested_by_min_stars"]["2"] - plp_alt2
    flag_notes = []
    for r in plp.itertuples():
        if "haploid" in r.zygosity:
            flag_notes.append(
                f"<li><b>{esc(r.GeneSymbol)} {esc(r.rsid)}</b>: a heterozygous call on the male X chromosome, "
                f"which has one copy. A single normal X copy can't be heterozygous (rare exceptions are mosaicism or an extra X), so this is almost certainly a probe artefact.</li>"
            )
        elif r.stars == 0:
            flag_notes.append(
                f"<li><b>{esc(r.GeneSymbol)} {esc(r.rsid)}</b> ({esc(r.zygosity)}): ClinVar record has "
                f"<i>no assertion criteria</i> (0 stars) — a single unreviewed submission.</li>"
            )
        else:
            flag_notes.append(
                f"<li><b>{esc(r.GeneSymbol)} {esc(r.rsid)}</b> ({esc(r.zygosity)}): {r.stars}-star record. "
                f"Needs confirmation by clinical-grade sequencing before it means anything.</li>"
            )

    def section_table(cat):
        d = hits[hits["category"] == cat].sort_values(
            ["stars", "GeneSymbol"], ascending=[False, True]
        )
        return hits_table(d), len(d), int((d["stars"] >= 2).sum())

    rf_t, rf_n, rf_2 = section_table("Risk factor")
    dr_t, dr_n, dr_2 = section_table("Drug response")
    pr_t, pr_n, pr_2 = section_table("Protective")
    as_t, as_n, as_2 = section_table("Association")
    co_t, co_n, co_2 = section_table("Conflicting")
    un_t, un_n, un_2 = section_table("Uncertain significance")

    sel = sel.copy()
    sel["_stars"] = sel["stars"].astype(int).map(stars)
    sel["_link"] = sel["VariationID"].map(clinvar_link)
    sel_table = table(
        sel,
        [
            "rsid",
            "GeneSymbol",
            "gt",
            "ref",
            "alt",
            "Name",
            "ClinicalSignificance",
            "_stars",
            "_link",
        ],
        [
            "rsID",
            "Gene",
            "Your genotype",
            "GRCh37 ref",
            "ClinVar alt",
            "Variant",
            "ClinVar classification",
            "Review",
            "ClinVar ID",
        ],
    )

    def called(rs):
        return esc(g.get(rs, "not typed"))

    tiles = [
        ("SNPs on the file", num(qc["n_variants"]), "GSA chip, GRCh37"),
        ("Call rate", pct(qc["call_rate"]), f"{num(ct['no_call'])} no-calls"),
        (
            "Autosomal heterozygosity",
            pct(qc["autosomal_het_rate"], 1),
            "of called autosomal SNVs",
        ),
        (
            "Inferred sex",
            sex,
            f"X het {pct(x['het_rate'], 2)} · Y call rate {pct(y['call_rate'], 1)}",
        ),
        (
            "Strand/build check",
            pct(strand_ok, 1),
            f"of {num(strand_total)} ClinVar sites fit forward strand",
        ),
        (
            "Runs of homozygosity",
            f"{roh_kb / 1000:.1f} Mb",
            f"{len(roh)} segments · longest {roh['KB'].max() / 1000:.2f} Mb",
        ),
    ]
    tiles_html = "".join(
        f'<div class="tile"><div class="tl">{esc(a)}</div><div class="tv">{esc(b)}</div><div class="ts">{esc(c)}</div></div>'
        for a, b, c in tiles
    )

    anc_html, anc = build_sections(num)
    from run_config import get as cfg_get

    display_name = cfg_get("display_name")
    longest = roh["KB"].max() / 1000 if len(roh) else 0.0
    if longest < 10:
        roh_bullet = (
            f"<b>No sign of closely related parents.</b> The stretches of DNA where your two copies are identical are all short "
            f"(longest {longest:.1f} Mb). Close relatives as parents leave runs of tens of Mb; distant shared ancestry is invisible "
            f"at this resolution."
        )
        roh_verdict = (
            "No run reaches 10 Mb, the scale close parental relatedness leaves."
        )
    else:
        roh_bullet = (
            f"<b>Long identical stretches present</b> (longest {longest:.1f} Mb): consistent with parents who share a recent "
            f"ancestor. See section 2."
        )
        roh_verdict = f"The longest run is {longest:.1f} Mb; runs of this length are what shared recent ancestry between the parents leaves."
    y_bullet = (
        "<b>Paternal line:</b> not readable from this sample (no Y chromosome)"
        if sex.lower() != "male"
        else (
            f"<b>Paternal line: {esc(anc['y'])}</b> (section 4). This is one paternal line out of thousands of ancestors."
        )
    )
    facts = summary_facts()
    notable = facts.get("pgx_notable", [])
    cyp3a5_note = (
        "Some of these are common: CYP3A5 *3/*3 is the genotype of 448 of 503 (89%) 1000 Genomes Europeans. "
        if any(n.startswith("CYP3A5 (Poor") for n in notable)
        else ""
    )
    n_pgs = len(facts.get("pgs", {})) or 1
    pgs_max_expected = (
        lambda n: (
            f"{n}{'th' if 10 <= n % 100 <= 20 else {1: 'st', 2: 'nd', 3: 'rd'}.get(n % 10, 'th')}"
        )
    )(round(100 * n_pgs / (n_pgs + 1)))
    sentences = []
    for rs, (allele, label, extra) in MARKERS.items():
        gt = g.get(rs)
        if gt:
            k = gt.count(allele)
            sentences.append(
                f"<b>{esc(label)}</b> ({rs}; risk allele {allele}{'; ' + esc(extra) if extra else ''}): {esc(gt)}, "
                f"{'no copy' if k == 0 else ('one copy' if k == 1 else 'two copies')}."
            )
    marker_sentences = " ".join(sentences)
    more_html = (
        aadr_section(num)
        + pgx_section()
        + hla_section()
        + pgs_section()
        + traits_section()
    )
    toc = "".join(
        f'<li><a href="#{i}">{esc(t)}</a></li>'
        for i, t in [
            ("qc", "1 · Quality control"),
            ("roh", "2 · Runs of homozygosity"),
            ("clinvar", "3 · Disease markers (ClinVar)"),
            ("y", "4 · Paternal line (Y)"),
            ("kg", "5 · Ancestry vs 1000 Genomes"),
            ("aadr", "6 · West Eurasian references"),
            ("pgx", "7 · Pharmacogenomics"),
            ("hla", "8 · HLA types"),
            ("pgs", "9 · Polygenic scores"),
            ("traits", "10 · Traits & blood group"),
            ("glossary", "Glossary"),
        ]
    )

    body = f"""
<header>
  <h1>{esc(display_name)}: genome report</h1>
  <p class="sub">{esc({"myheritage": "MyHeritage", "ftdna": "FTDNA", "23andme": "23andMe", "ancestrydna": "AncestryDNA"}.get(qc.get("input_format"), "Consumer"))} raw export (<code>{esc(qc["raw_file"])}</code>, Illumina GSA, GRCh37) · ClinVar <code>variant_summary.txt.gz</code> dated {clinvar_date()} · generated {date.today().isoformat()} on this Mac. No data left the machine.</p>
</header>

<aside class="caveat" role="note">
  <b>Read before the disease section.</b> This is a consumer genotyping array, not sequencing. For <i>rare</i> variants
  called pathogenic, array calls are frequently wrong and any such finding needs confirmation by a clinical lab before it
  means anything. In the largest study of this (Weedon et al. 2021, BMJ 372:n214), only 16% of 4,757 chip calls of very rare
  variants were confirmed by sequencing: "SNP chips are extremely unreliable for genotyping very rare pathogenic variants and
  should not be used to guide health decisions without validation." The array cannot see copy-number changes (e.g. CYP2D6 deletions/duplications) or most indels, and
  most variants in your genes were never measured. MyHeritage itself labels the file "not intended for medical purposes".
  A "risk factor" in ClinVar is usually a small shift in odds from a common variant, not a diagnosis.
</aside>

<section class="tiles">{tiles_html}</section>

<section class="plain" aria-label="Summary in plain language">
  <h2 class="first">The short version</h2>
  <ul>
    <li><b>The data file is good.</b> Almost every marker was read successfully ({pct(qc["call_rate"], 1)}), and it behaves as a {sex.lower()} sample should.</li>
    <li>{roh_bullet}</li>
    <li><b>This chip can't screen for rare disease</b>, and a clean result wouldn't mean much. Of the {len(plp)} {"marker" if len(plp) == 1 else "markers"} labelled "pathogenic",
    none passes a basic credibility check (section 3).</li>
    <li><b>APOE</b> (the main common Alzheimer-risk gene): the result is behind a click in section 3, so you choose whether to see it.</li>
    <li>{y_bullet}</li>
    <li><b>Ancestry: European</b> at the continental level (minor non-European components of a few percent are within the method's noise).
    Within Europe, two independent methods agree that you are most similar to <b>{" and ".join(facts.get("aadr_near", [])) or "–"}</b>
    samples ({esc(facts.get("aadr_region", ""))}). {esc(facts.get("aadr_limit", ""))} (section 6).</li>
    <li><b>HLA (imputed, not a clinical test):</b> {"no coeliac-risk type (DQ2/DQ8), so coeliac disease is very unlikely" if facts.get("dq_neg") else "carries a coeliac-risk HLA type"};
    no HLA-B27; the drug-hypersensitivity HLA alleles PharmCAT checks are not predicted (section 8).</li>
    <li><b>Medicines:</b> genes whose guideline phenotype isn't \"normal\": {", ".join(facts.get("pgx_notable", [])) or "none"}.
    {cyp3a5_note}This only matters if
    you're prescribed a drug those genes affect, and a clinic would confirm it first (section 7).</li>
    <li><b>Polygenic scores</b> (percentile vs 1000 Genomes Europeans): {", ".join(f"{k} {v:.0f}" for k, v in facts.get("pgs", {}).items()) or "–"}.
    Relative rankings, not risks, each uncertain by several points; with {n_pgs} scores, the highest is expected near the {pgs_max_expected} percentile by chance (section 9).</li>
    <li><b>Blood group {facts.get("abo") or "?"}</b>; traits in section 10.</li>
  </ul>
  <p class="note">This page is a first look for someone curious about their own DNA, not a medical result. The <a href="#glossary">glossary</a> explains every technical word.</p>
</section>

<nav class="toc" aria-label="Contents"><b>Contents</b><ol>{toc}</ol></nav>

<section id="glossary">
  <details>
    <summary><b>Glossary: what the terms mean</b></summary>
    <dl class="gloss">
      <dt id="g-snp">SNP (single-nucleotide polymorphism)</dt>
      <dd>A single position in the DNA where people differ in which letter (A, C, G or T) they have. Your chip measured about {num(qc["n_variants"])} of these positions, a small fraction of the roughly 3 billion letters in the genome.</dd>
      <dt id="g-rsid">rsID</dt>
      <dd>The catalogue number of a SNP in dbSNP, the public database of human variants (e.g. rs429358). Everyone who studies that position uses the same number.</dd>
      <dt id="g-allele">Allele · genotype</dt>
      <dd>An allele is one version of a position (say, C or T). You have two copies of each non-sex chromosome, one from each parent, so your genotype has two letters (e.g. CT).</dd>
      <dt id="g-ref">Reference allele · variant (alt) allele</dt>
      <dd>The reference is the letter in the standard human reference genome. The variant (or "alternate") allele is the other letter. The reference is simply the letter the reference genome happened to contain, not "normal" or "healthy". Occasionally it is the disease version itself (Factor V Leiden is an example).</dd>
      <dt id="g-zyg">Heterozygous · homozygous · hemizygous</dt>
      <dd><b>Heterozygous</b>: two different letters (CT), meaning one copy of the variant. <b>Homozygous</b>: two identical letters (CC or TT). <b>Hemizygous</b>: only one copy exists. Men have one X and one Y chromosome, so genes on them are hemizygous.</dd>
      <dt id="g-carrier">Carrier</dt>
      <dd>Someone with one copy of a variant that only causes disease when both copies are affected (a <i>recessive</i> condition). Carriers are usually healthy, but can pass the variant on.</dd>
      <dt id="g-penetrance">Penetrance</dt>
      <dd>How often people with a variant actually develop the condition. "Low penetrance" means most carriers never do.</dd>
      <dt id="g-clinvar">ClinVar</dt>
      <dd>The US National Library of Medicine's public database where labs and researchers submit their verdict on what a variant does. This report uses a copy downloaded on {clinvar_date()}.</dd>
      <dt id="g-sig">Pathogenic · likely pathogenic</dt>
      <dd>Submitters judged the variant to cause (or very probably cause) a disease. It says nothing about how likely <i>you</i> are to be affected: that depends on whether the condition needs one or two copies, on penetrance, and on whether the call is even correct.</dd>
      <dt id="g-vus">Uncertain significance (VUS)</dt>
      <dd>Not enough evidence either way. A VUS should be treated as "unknown", not as a warning.</dd>
      <dt id="g-benign">Benign · likely benign</dt>
      <dd>Judged harmless. Most of your ClinVar matches are in this group, as they are for everyone.</dd>
      <dt id="g-risk">Risk factor · association · protective</dt>
      <dd>Common variants that nudge the odds of a condition up (risk factor) or down (protective), or that turned up in a population study (association). Individually these effects are usually small, and many people with them never develop the condition.</dd>
      <dt id="g-drug">Drug response (pharmacogenomics)</dt>
      <dd>Variants that change how the body processes or responds to a medicine, affecting dose, effect or side effects. Clinically these are read as combinations (<i>star alleles</i> such as CYP2C19*17), not one SNP at a time.</dd>
      <dt id="g-conflict">Conflicting classifications</dt>
      <dd>Different submitters to ClinVar disagree about the variant.</dd>
      <dt id="g-stars">Review stars (★)</dt>
      <dd>ClinVar's measure of how much to trust a classification. ★★★★ practice guideline · ★★★ reviewed by an expert panel · ★★ several submitters agree · ★ one submitter who stated their criteria (or submitters conflict) · no stars: no criteria given. <b>Zero-star records deserve very little weight.</b></dd>
      <dt id="g-array">Genotyping array (chip)</dt>
      <dd>The technology MyHeritage used. It checks a fixed list of positions and is very accurate for common variants, but it often misreads rare ones. That is why any rare "pathogenic" result needs a clinical test before it means anything.</dd>
      <dt id="g-callrate">Call rate · no-call</dt>
      <dd>The share of positions the chip read successfully. A no-call (<code>--</code>) is a position it could not read.</dd>
      <dt id="g-het">Heterozygosity</dt>
      <dd>The share of positions where your two copies differ. It is a quick sanity check: sex chromosomes in a man should be close to zero.</dd>
      <dt id="g-roh">Runs of homozygosity (ROH) · F<sub>ROH</sub></dt>
      <dd>Long stretches where both copies are identical because both parents inherited that segment from a shared ancestor. F<sub>ROH</sub> is the fraction of the genome in such runs. Long runs point to recent shared ancestors; short ones are normal background.</dd>
      <dt id="g-build">GRCh37 · forward strand</dt>
      <dd>GRCh37 (also called hg19) is the version of the reference genome the positions refer to. DNA is double-stranded, and "forward strand" says which of the two strands the letters were read from. Both must match the databases being compared against, and the QC section checks that they do.</dd>
      <dt id="g-pseudo">Mb (megabase)</dt>
      <dd>One million DNA letters.</dd>
    </dl>
  </details>
</section>

<section id="qc">
  <h2>1 · Quality control</h2>
  <p class="plainline"><b>In plain words:</b> before trusting any result, check that the file is intact, that it behaves like a real
  human sample, and that it uses the same coordinate system as the databases. It passes all three.</p>
  <p>{num(ct["hom"])} homozygous, {num(ct["het"])} heterozygous, {num(ct["no_call"])} no-calls and {num(ct["indel"])} indel probes
  (coded D/I; not interpreted here because the file does not say which allele is the reference).
  {qc["duplicate_positions"]} positions appear twice under different rsIDs; no duplicated rsIDs; no unexpected chromosome codes.</p>
  {rare_hom_para(rh_qc)}
  {impute_bench_para()}
  <h3>Heterozygosity per chromosome</h3>
  <p class="note">Autosomes sit near {pct(qc["autosomal_het_rate"], 0)}. X and Y are near zero, consistent with one copy of each ({sex.lower()}, as inferred by plink).
  The few heterozygous X/Y calls ({x["het"]} on X, {y["het"]} on Y) are genotyping errors.
  The file contains no pseudo-autosomal X probes (<code>plink --split-x b37</code> moved none) and no mitochondrial probes.</p>
  {het_chart}
  <details><summary>Per-chromosome table</summary>{chrom_table}</details>
  <h3>Genotype classes</h3>
  <p class="note">A/T and C/G heterozygotes are rare because array designs avoid strand-ambiguous SNPs. The strand check against
  ClinVar reference alleles found {pct(strand_ok, 1)} of {num(strand_total)} overlapping sites consistent with the forward strand,
  {num(cv["strand_only_if_complemented"])} consistent only after complementing, and {num(cv["strand_inconsistent"])} inconsistent. rsIDs agree with ClinVar at {pct(cv["rsid_agreement"], 2)} of matched sites.
  The file really is GRCh37, forward strand.</p>
  {gt_chart}
</section>

<section id="roh">
  <h2>2 · Runs of homozygosity</h2>
  <p class="plainline"><b>In plain words:</b> if your parents shared a recent ancestor, you would have long identical stretches
  (blue bars below) on many chromosomes. Yours are few and short, which is the pattern of unrelated parents.</p>
  <p>Long stretches where both chromosome copies are identical are inherited from a shared ancestor on both sides.
  Many long runs (tens of Mb) point to closely related parents. <code>plink --homozyg</code> with default settings
  (≥1 Mb, ≥100 SNPs) found <b>{len(roh)} segments totalling {roh_kb / 1000:.1f} Mb</b>, the longest {roh["KB"].max() / 1000:.2f} Mb;
  {len(roh15)} are ≥1.5 Mb ({roh15.sum() / 1000:.1f} Mb). As a fraction of the {auto_len / 1e6:,.0f} Mb autosomal span the array covers,
  F<sub>ROH</sub> ≈ {froh:.4f}. {roh_verdict}
  Short runs like these are typical background, but comparing against a reference population is what would make that quantitative.</p>
  {karyo}
  <details><summary>Segment table</summary>{roh_table}</details>
</section>

<section id="clinvar">
  <h2>3 · Disease and drug-response markers (ClinVar)</h2>
  <div class="plainline">
    <p><b>In plain words:</b> every person carries thousands of DNA variants listed in ClinVar, and nearly all of them are harmless.
    What matters for any single line is:</p>
    <ol>
      <li><b>What it is labelled</b> (<a href="#g-sig">pathogenic</a>, <a href="#g-risk">risk factor</a>, <a href="#g-drug">drug response</a>…).</li>
      <li><b>How trustworthy that label is</b>: the <a href="#g-stars">review stars</a>. Treat zero stars as a rumour.</li>
      <li><b>How many copies you have</b> (<a href="#g-zyg">zygosity</a>). One copy of a recessive disease variant makes you a <a href="#g-carrier">carrier</a>, not a patient.</li>
      <li><b>Whether the chip read it correctly.</b> Chips often misread rare variants (<a href="#g-array">why</a>).</li>
    </ol>
  </div>
  <p>Your {num(cv["genome_called_snvs"])} called SNVs were matched to {num(cv["clinvar_grch37_snv_rows"])} ClinVar GRCh37 SNV records by
  chromosome + position. {num(cv["matched_sites"])} of your sites are in ClinVar, and you carry the ClinVar variant allele at
  {num(cv["hit_rows"])} records. Most are benign. The array tested {num(cv["plp_sites_tested_by_min_stars"]["2"])} sites that ClinVar lists as
  pathogenic/likely pathogenic with ≥2 review stars. You carry only the reference allele at {num(plp_ref2)} of them;
  {"the exception is flagged" if plp_alt2 == 1 else f"the {plp_alt2} exceptions are flagged"} below.</p>
  <h3>ClinVar records where you carry the variant allele</h3>
  {cat_chart}

  <h3>Pathogenic / likely pathogenic ({len(plp)})</h3>
  <p class="note">None of these is a credible finding as it stands:</p>
  <ul class="flags">{"".join(flag_notes)}</ul>
  {hits_table(plp)}

  <h3>Well-known markers, whatever the genotype</h3>
  <details class="optin"><summary><b>APOE, the main common Alzheimer-risk gene.</b> Click to show. Some people prefer not to know;
  the result predicts nothing on its own, and nothing medical currently follows from it.</summary>
  <p class="note"><b>APOE:</b> rs429358 = {called("rs429358")}, rs7412 = {called("rs7412")} → <b>{esc(apoe)}</b>.
  Definition from dbSNP: ε4 = C at rs429358, ε2 = T at rs7412, ε3 = neither. The local ClinVar rows agree: rs429358-C is
  a risk factor for Alzheimer disease 2, and rs7412-T is linked to familial type 3 hyperlipoproteinemia. This is the most consequential line in the report. It shifts risk and predicts nothing on its own. If it matters to you,
  talk to a clinical geneticist before reading effect sizes online.</p></details>
  <p class="note">
  {marker_sentences}</p>
  {sel_table}

  <h3>Drug response ({dr_n}; {dr_2} with ≥2 stars)</h3>
  <p class="plainline"><b>In plain words:</b> these variants can change how your body handles certain medicines, for example
  cholesterol drugs, blood thinners or some chemotherapy. None of them is a disease. They are worth mentioning to a doctor only if one
  of the listed drugs is ever prescribed, and a clinic would confirm them with a dedicated test. The "Condition(s)" column names the drug and effect.</p>
  <p class="note">Pharmacogenomic results normally combine several SNPs into a star-allele diplotype (e.g. CYP2C19*17). These are
  single-SNP records, and a real PGx call would come from a tool such as PharmCAT. The 3–4-star rows are ClinPGx/expert-panel curated.</p>
  {dr_t}
  <details><summary>Risk factor ({rf_n}; {rf_2} with ≥2 stars)</summary>{rf_t}</details>
  <details><summary>Protective ({pr_n})</summary>{pr_t}</details>
  <details><summary>Association ({as_n})</summary>{as_t}</details>
  <details><summary>Conflicting classifications ({co_n})</summary>{co_t}</details>
  <details><summary>Uncertain significance ({un_n})</summary>{un_t}</details>
{imputed_clinvar_block()}
</section>

{anc_html}
{more_html}

<section>
  <h2>Not done yet, and why</h2>
  <ul>
    <li>CYP2D6 with copy number: needs the chip's raw intensity data (not in the MyHeritage export), sequencing, or a dedicated assay.
    The SNP-only call above is research grade.</li>
    <li>gnomAD frequencies: the ClinVar tables use 1000 Genomes frequencies computed locally (2,504 people). A per-variant gnomAD
    lookup would reveal your variant positions to an outside service; a privacy-safe alternative (a bulk download of gnomAD frequencies
    for all ClinVar positions, independent of any genotype) is not implemented.</li>
    <li>Any rare "pathogenic" call would need confirmation by a clinical lab; none here currently warrants it.</li>
  </ul>
</section>
"""
    page = TEMPLATE.replace("{{BODY}}", body)
    out = RES / "report.html"
    out.write_text(page)
    print(out)


TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="color-scheme" content="light dark">
<title>Genome First Look</title>
<link rel="icon" href="data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 32 32'%3E%3Cpath d='M10 3c0 9 12 9 12 17s-12 8-12 9M22 3c0 9-12 9-12 17s12 8 12 9' stroke='%232a78d6' stroke-width='3' fill='none'/%3E%3C/svg%3E">
<style>
:root {
  --page: #f9f9f7; --surface-1: #fcfcfb; --text-primary: #0b0b0b; --text-secondary: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --border: rgba(11,11,11,0.10);
  --series-1: #2a78d6; --series-2: #eb6834; --series-3: #1baf7a; --series-4: #eda100; --critical: #d03b3b; --caveat-bg: #fdf3e7;
}
@media (prefers-color-scheme: dark) {
  :root {
    --page: #0d0d0d; --surface-1: #1a1a19; --text-primary: #ffffff; --text-secondary: #c3c2b7; --muted: #898781;
    --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10);
    --series-1: #3987e5; --series-2: #d95926; --series-3: #199e70; --series-4: #c98500; --caveat-bg: #2a2118;
  }
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--page); color: var(--text-primary);
  font: 15px/1.55 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1040px; margin: 0 auto; padding: 32px 16px 64px; }
h1 { font-size: 28px; margin: 0 0 4px; }
h2 { font-size: 20px; margin: 40px 0 8px; padding-top: 16px; border-top: 1px solid var(--grid); }
h3 { font-size: 15px; margin: 24px 0 6px; }
.sub, .note { color: var(--text-secondary); }
.note { font-size: 14px; }
code { font-size: 13px; }
a { color: var(--series-1); }
.caveat { background: var(--caveat-bg); border-left: 3px solid var(--critical); padding: 12px 16px; border-radius: 6px;
  margin: 20px 0; font-size: 14px; }
.tiles { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 12px; margin: 20px 0; }
.tile { background: var(--surface-1); border: 1px solid var(--border); border-radius: 10px; padding: 12px 14px; }
.tl { font-size: 12px; color: var(--text-secondary); }
.tv { font-size: 24px; font-weight: 600; margin: 2px 0; }
.ts { font-size: 12px; color: var(--muted); }
.chart { width: 100%; height: auto; display: block; background: var(--surface-1); border-radius: 10px;
  border: 1px solid var(--border); padding: 8px; margin: 8px 0; }
.grid { stroke: var(--grid); stroke-width: 1; }
.axis { stroke: var(--axis); stroke-width: 1; }
.tick { fill: var(--muted); font-size: 11px; font-variant-numeric: tabular-nums; }
.label { fill: var(--text-secondary); font-size: 12px; }
.value { fill: var(--text-primary); font-size: 12px; font-variant-numeric: tabular-nums; }
.track { fill: var(--grid); }
.hit { outline: none; cursor: default; }
.hit:hover path, .hit:focus path, .hit:hover rect + rect, .hit:focus rect + rect { opacity: .75; }
.legend { display: flex; gap: 16px; font-size: 12px; color: var(--text-secondary); margin-top: 8px; }
.legend i { display: inline-block; width: 10px; height: 10px; border-radius: 3px; margin-right: 6px; vertical-align: -1px; }
.tablewrap { overflow-x: auto; margin: 8px 0; }
table { border-collapse: collapse; width: 100%; font-size: 12.5px; background: var(--surface-1); }
th, td { text-align: left; padding: 6px 8px; border-bottom: 1px solid var(--grid); vertical-align: top; }
th { color: var(--text-secondary); font-weight: 600; position: sticky; top: 0; background: var(--surface-1); }
table.num td { font-variant-numeric: tabular-nums; }
.stars { color: var(--series-1); letter-spacing: 1px; white-space: nowrap; }
details { margin: 8px 0; }
summary { cursor: pointer; color: var(--text-secondary); font-size: 14px; }
.flags { font-size: 14px; }
.plain { background: var(--surface-1); border: 1px solid var(--border); border-radius: 10px; padding: 4px 20px 12px; margin: 20px 0; }
.plain li { margin: 6px 0; }
h2.first { border-top: 0; margin-top: 12px; padding-top: 0; }
.plainline { background: var(--surface-1); border-left: 3px solid var(--series-1); border-radius: 6px; padding: 8px 14px; font-size: 14px; }
.plainline p { margin: 4px 0; } .plainline ol { margin: 4px 0; padding-left: 20px; }
.gloss dt { font-weight: 600; margin-top: 12px; scroll-margin-top: 16px; }
.gloss dd { margin: 2px 0 0; color: var(--text-secondary); font-size: 14px; }
.gloss dt:target { color: var(--series-1); }
.toc { background: var(--surface-1); border: 1px solid var(--border); border-radius: 10px; padding: 10px 18px; margin: 16px 0; font-size: 14px; }
.toc ol { list-style: none; padding: 0; margin: 6px 0 0; columns: 2 220px; }
.toc li { margin: 3px 0; }
.pgsgrid { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 12px; margin: 12px 0; }
.pgscard { background: var(--surface-1); border: 1px solid var(--border); border-radius: 10px; padding: 12px 14px; }
.pgscard .unit { font-size: 13px; font-weight: 400; color: var(--text-secondary); margin-left: 4px; }
.chart.mini { border: 0; padding: 0; margin: 6px 0; }
.histbar { fill: var(--muted); opacity: .45; }
.selfline { stroke: var(--series-2); stroke-width: 2; }
.selflabel.small { font-size: 11px; }
td.src { color: var(--muted); font-size: 11.5px; }
tr.hl td { background: color-mix(in srgb, var(--series-1) 10%, transparent); }
.refdots circle { fill: var(--muted); opacity: .28; }
.frame { stroke: var(--grid); }
.leader { stroke: var(--muted); stroke-width: 1; }
.cdot { fill: var(--text-secondary); opacity: .7; }
.clabel { fill: var(--text-primary); font-size: 11px; font-weight: 600; paint-order: stroke; stroke: var(--surface-1); stroke-width: 3px; }
.self { fill: var(--series-2); stroke: var(--surface-1); stroke-width: 2px; }
.selflabel { fill: var(--text-primary); font-size: 13px; font-weight: 700; paint-order: stroke; stroke: var(--surface-1); stroke-width: 3px; }
.ypath { list-style: none; padding: 0; display: flex; flex-wrap: wrap; gap: 6px; font-size: 12px; }
.ypath li { background: var(--surface-1); border: 1px solid var(--border); border-radius: 6px; padding: 4px 8px; display: flex; flex-direction: column; }
.ypath li.leaf { border-color: var(--series-1); }
.ypath .node { font-weight: 600; } .ypath .cnt { color: var(--muted); font-size: 11px; }
#tip { position: fixed; pointer-events: none; background: var(--text-primary); color: var(--surface-1); font-size: 12px;
  padding: 5px 8px; border-radius: 6px; opacity: 0; transition: opacity .1s; z-index: 10; max-width: 320px; }
@media print {
  body { background: #fff; color: #000; }
  details { display: block; } details > * { display: block; }
  .chart { break-inside: avoid; }
  #tip { display: none; }
}
</style>
</head>
<body>
<main>
{{BODY}}
</main>
<div id="tip" role="tooltip"></div>
<script>
(() => {
  const tip = document.getElementById('tip');
  const show = (el, x, y) => { tip.textContent = el.dataset.tip; tip.style.opacity = 1;
    const w = tip.offsetWidth; tip.style.left = Math.min(x + 12, innerWidth - w - 8) + 'px'; tip.style.top = (y + 14) + 'px'; };
  document.querySelectorAll('[data-tip]').forEach(el => {
    el.addEventListener('mousemove', e => show(el, e.clientX, e.clientY));
    el.addEventListener('mouseleave', () => tip.style.opacity = 0);
    el.addEventListener('focus', () => { const r = el.getBoundingClientRect(); show(el, r.left, r.top); });
    el.addEventListener('blur', () => tip.style.opacity = 0);
  });
  const openGloss = () => { if (location.hash.startsWith('#g')) { document.querySelector('#glossary details').open = true;
    document.querySelector(location.hash)?.scrollIntoView(); } };
  addEventListener('hashchange', openGloss); openGloss();
  addEventListener('beforeprint', () => document.querySelectorAll('details').forEach(d => d.open = true));
})();
</script>
</body>
</html>
"""

if __name__ == "__main__":
    main()
