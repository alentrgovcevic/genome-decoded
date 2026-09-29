"""Y-haplogroup and ancestry sections for report.py. Returns HTML fragments; reads results/yhaplo and results/ancestry."""

import html
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path.cwd()  # the run directory runs/<sample>/ (see workflow/Snakefile)
RES = ROOT / "results"
esc = html.escape

SUPERPOP_NAMES = {
    "AFR": "African",
    "AMR": "Admixed American",
    "EAS": "East Asian",
    "EUR": "European",
    "SAS": "South Asian",
}
# Adjacent-validated categorical slots 1–4 (dataviz reference palette), for the ADMIXTURE stack.
STACK_VARS = ["--series-1", "--series-2", "--series-3", "--series-4"]


def load_y():
    hg = (RES / "yhaplo" / "haplogroups.y.txt").read_text().split()
    path = (RES / "yhaplo" / "paths.y.txt").read_text().split("|")[1].split()
    derived = (RES / "yhaplo" / "derived.snps.y.txt").read_text().split()
    return {
        "hg_snp": hg[2],
        "ycc": hg[3],
        "path": [p.split(":") for p in path],
        "derived": [d.split(":") for d in derived if ":" in d],
    }


def y_section(y, n_y_calls):
    steps = "".join(
        f'<li class="{"leaf" if i == len(y["path"]) - 1 else ""}"><span class="node">{esc(n)}</span>'
        f'<span class="cnt">{esc(c)} derived</span></li>'
        for i, (n, c) in enumerate(y["path"])
    )
    ycc_to_snps = {}
    for node, snp in y["derived"]:
        ycc_to_snps.setdefault(node, []).append(snp)
    rows = "".join(
        f"<tr><td>{esc(k)}</td><td>{esc(', '.join(v[:12]))}{' …' if len(v) > 12 else ''}</td><td>{len(v)}</td></tr>"
        for k, v in ycc_to_snps.items()
        if k[:1] == y["ycc"][:1]  # branches within the called major haplogroup
    )
    return steps, rows


def load_pca(label):
    s = pd.read_csv(RES / "ancestry" / f"proj_{label}.sscore", sep="\t")
    s = s.rename(columns={"#FID": "FID", "#IID": "IID"})
    pcs = [c for c in s.columns if c.endswith("_AVG")]
    pops = pd.read_csv(
        RES / "ancestry" / "pops.tsv", sep="\t", names=["IID", "super", "pop"]
    )
    s = s.merge(pops, on="IID", how="left")
    return s, pcs


def nearest(s, pcs, group_col, k=4, n_pc=4):
    ref = s[s["IID"] != "SELF"]
    me = s.loc[s["IID"] == "SELF", pcs[:n_pc]].to_numpy()[0]
    cent = ref.groupby(group_col)[pcs[:n_pc]].mean()
    spread = ref.groupby(group_col)[pcs[:n_pc]].std().mean(axis=1)
    d = np.sqrt(((cent.to_numpy() - me) ** 2).sum(axis=1))
    out = pd.DataFrame(
        {"group": cent.index, "dist": d, "spread": spread.to_numpy()}
    ).sort_values("dist")
    return out.head(k)


def label_short(g: str) -> str:
    return g.replace("_", " ")


def scatter(
    s,
    pcs,
    group_col,
    label_fn,
    title_x="PC1",
    title_y="PC2",
    width=760,
    height=460,
    robust=False,
):
    """robust=True clips the axes to the central 98% of reference points (outliers drawn on the frame)."""
    x, y = pcs[0], pcs[1]
    left, right, top, bottom = 48, 16, 16, 36
    pw, ph = width - left - right, height - top - bottom
    ref0 = s[s["IID"] != "SELF"]
    q = (0.01, 0.99) if robust else (0, 1)
    xmin, xmax = ref0[x].quantile(q[0]), ref0[x].quantile(q[1])
    ymin, ymax = ref0[y].quantile(q[0]), ref0[y].quantile(q[1])
    me0 = s[s["IID"] == "SELF"].iloc[0]
    xmin, xmax = min(xmin, me0[x]), max(xmax, me0[x])
    ymin, ymax = min(ymin, me0[y]), max(ymax, me0[y])
    padx, pady = (xmax - xmin) * 0.05, (ymax - ymin) * 0.05
    xmin, xmax, ymin, ymax = xmin - padx, xmax + padx, ymin - pady, ymax + pady
    clip = lambda v, a, b: min(max(v, a), b)
    sx = lambda v: left + pw * (clip(v, xmin, xmax) - xmin) / (xmax - xmin)
    sy = lambda v: top + ph - ph * (clip(v, ymin, ymax) - ymin) / (ymax - ymin)
    parts = [
        f'<svg viewBox="0 0 {width} {height}" class="chart" role="img" aria-label="PCA scatter">'
    ]
    parts.append(
        f'<rect x="{left}" y="{top}" width="{pw}" height="{ph}" fill="none" class="frame"/>'
    )
    ref = s[s["IID"] != "SELF"]
    dots = "".join(
        f'<circle cx="{sx(a):.1f}" cy="{sy(b):.1f}" r="2.2"/>'
        for a, b in zip(ref[x], ref[y])
    )
    parts.append(f'<g class="refdots">{dots}</g>')
    me = s[s["IID"] == "SELF"].iloc[0]
    mx, my = sx(me[x]), sy(me[y])
    # Greedy label placement: boxes are (x0, y0, x1, y1); the self marker is reserved first.
    taken = [(mx - 9, my - 9, mx + 9, my + 9)]

    def place(cx, cy, w, h, offsets):
        for dx, dy in offsets:
            b = (cx + dx - w / 2, cy + dy - h / 2, cx + dx + w / 2, cy + dy + h / 2)
            inside = (
                b[0] >= left
                and b[2] <= width - right
                and b[1] >= top
                and b[3] <= top + ph
            )
            if inside and not any(
                b[0] < t[2] and b[2] > t[0] and b[1] < t[3] and b[3] > t[1]
                for t in taken
            ):
                taken.append(b)
                return cx + dx, cy + dy
        return None

    you_x, you_y = place(
        mx, my, 34, 16, [(28, -12), (-28, -12), (28, 12), (-28, 12), (0, -22), (0, 22)]
    ) or (mx + 28, my - 12)
    offsets = [
        (0, 0),
        (0, -18),
        (0, 18),
        (-44, 0),
        (44, 0),
        (0, -36),
        (0, 36),
        (-44, -18),
        (44, 18),
    ]
    # Label populations nearest to SELF first; a label that cannot be placed without overlap is skipped
    # (its hover target stays, and the nearest-population table lists it).
    cents = ref.groupby(group_col)[[x, y]].mean()
    order = ((cents[x] - me[x]) ** 2 + (cents[y] - me[y]) ** 2).sort_values().index
    for g in order:
        sub = ref[ref[group_col] == g]
        cx, cy = sx(sub[x].mean()), sy(sub[y].mean())
        spot = place(cx, cy, 8 + 6.6 * len(label_short(g)), 15, offsets)
        tip = f"{label_fn(g)} · {len(sub)} samples"
        if spot is None:
            parts.append(
                f'<g class="hit" tabindex="0" data-tip="{esc(tip)}"><circle cx="{cx:.1f}" cy="{cy:.1f}" r="4" class="cdot"/></g>'
            )
            continue
        lx, ly = spot
        leader = (
            f'<line x1="{cx:.1f}" y1="{cy:.1f}" x2="{lx:.1f}" y2="{ly:.1f}" class="leader"/>'
            if (lx, ly) != (cx, cy)
            else ""
        )
        parts.append(
            f'<g class="hit" tabindex="0" data-tip="{esc(tip)}">{leader}<rect x="{lx - 26:.1f}" y="{ly - 10:.1f}" width="52" height="20" fill="transparent"/>'
            f'<text x="{lx:.1f}" y="{ly + 4:.1f}" class="clabel" text-anchor="middle">{esc(label_short(g))}</text></g>'
        )
    parts.append(
        f'<g class="hit" tabindex="0" data-tip="You · {title_x} {me[x]:.4f}, {title_y} {me[y]:.4f}">'
        f'<circle cx="{mx:.1f}" cy="{my:.1f}" r="12" fill="transparent"/>'
        f'<circle cx="{mx:.1f}" cy="{my:.1f}" r="6" class="self"/></g>'
    )
    parts.append(
        f'<text x="{you_x:.1f}" y="{you_y + 5:.1f}" class="selflabel" text-anchor="middle">You</text>'
    )
    parts.append(
        f'<text x="{left + pw / 2}" y="{height - 8}" class="tick" text-anchor="middle">{title_x}</text>'
    )
    parts.append(
        f'<text x="14" y="{top + ph / 2}" class="tick" text-anchor="middle" transform="rotate(-90 14 {top + ph / 2})">{title_y}</text>'
    )
    parts.append("</svg>")
    return "".join(parts)


def load_admixture():
    cols = (RES / "ancestry" / "admixture_columns.txt").read_text().split()
    q = pd.read_csv(
        RES / "ancestry" / "admixture_Q.tsv", sep=r"\s+", names=["IID"] + cols
    )
    pops = pd.read_csv(
        RES / "ancestry" / "pops.tsv", sep="\t", names=["IID", "super", "pop"]
    )
    q = q.merge(pops, on="IID", how="left")
    me = q[q["IID"] == "SELF"].iloc[0]
    ref_mean = q[q["IID"] != "SELF"].groupby("super")[cols].mean()
    return cols, me, ref_mean


def stack_bar(cols, values, width=760, height=44):
    order = sorted(range(len(cols)), key=lambda i: cols[i])
    parts = [
        f'<svg viewBox="0 0 {width} {height}" class="chart" role="img" aria-label="Ancestry proportions">'
    ]
    x = 8
    pw = width - 16
    for i in order:
        w = pw * values[i]
        if w <= 0:
            continue
        gap = 2 if w > 4 else 0
        tip = f"{SUPERPOP_NAMES.get(cols[i], cols[i])}: {values[i] * 100:.1f}%"
        parts.append(
            f'<g class="hit" tabindex="0" data-tip="{esc(tip)}"><rect x="{x:.1f}" y="8" width="{max(w - gap, 1):.1f}" height="{height - 16}" rx="4" '
            f'fill="var({STACK_VARS[order.index(i)]})"/></g>'
        )
        x += w
    parts.append("</svg>")
    legend = (
        '<div class="legend">'
        + "".join(
            f'<span><i style="background:var({STACK_VARS[k]})"></i>{esc(SUPERPOP_NAMES.get(cols[i], cols[i]))} {values[i] * 100:.1f}%</span>'
            for k, i in enumerate(order)
        )
        + "</div>"
    )
    return legend + "".join(parts)


POP_NAMES = {  # IGSR README_populations.md
    "CEU": "Utah residents with Northern & Western European ancestry",
    "TSI": "Tuscans (Italy)",
    "GBR": "British (England & Scotland)",
    "FIN": "Finnish",
    "IBS": "Iberians (Spain)",
}


def y_html() -> str:
    """Section 4. Female samples (no Y) get a stub; notes come from haplogroup_notes.json when a key SNP is derived."""
    from run_config import sex

    title = "4 · Paternal line: Y-chromosome haplogroup"
    if sex() != "male" or not (RES / "yhaplo" / "haplogroups.y.txt").exists():
        why = (
            "This sample has no Y chromosome (inferred female), so the paternal line can't be read from it. A father's, brother's "
            "or paternal uncle's DNA would carry it."
            if sex() == "female"
            else "Not run: scripts/y_haplogroup.py"
        )
        return f'<section id="y"><h2>{title}</h2><p class="note">{why}</p></section>'
    y = load_y()
    steps, yrows = y_section(y, None)
    derived_snps = {s for _, s in y["derived"]}
    notes = json.loads(
        (Path(__file__).resolve().parent / "haplogroup_notes.json").read_text()
    )["notes"]
    shown = [n for n in notes if n["key_snp"] in derived_snps]
    note_html = "".join(
        f'<p class="note"><b>About {esc(n["clade"])}</b> ({esc(n["source"])}): {esc(n["text"])} This is the lineage of one '
        f"paternal ancestor, not a statement about the rest of your ancestry.</p>"
        for n in shown
    )
    last_node, last_count = y["path"][-1]
    thin = (
        f" The last step rests on {last_count} SNP, so treat that level as provisional."
        if last_count in ("1", 1)
        else ""
    )
    return f"""<section id="y">
  <h2>{title}</h2>
  <p class="plainline"><b>In plain words:</b> the Y chromosome passes from father to son almost unchanged, so it traces a single line:
  your father's father's father… back thousands of years. Mutations along the way sort men into branches of one big family tree
  (haplogroups). It covers one ancestor per generation out of many, so it says nothing about the rest of your ancestry.</p>
  <div class="tiles"><div class="tile"><div class="tl">Y haplogroup</div><div class="tv">{esc(y["hg_snp"])}</div>
  <div class="ts">ISOGG 2016 name {esc(y["ycc"])}</div></div></div>
  <p>Called with 23andMe's <code>yhaplo</code> (ISOGG 2016 tree) from the Y-chromosome calls. The path runs from the root of the tree to
  the called branch; the number is how many SNPs confirm each step.{thin}</p>
  <ol class="ypath">{steps}</ol>
  {note_html}
  <details><summary>Derived SNPs observed along the path</summary><div class="tablewrap"><table>
  <thead><tr><th>Branch</th><th>SNPs</th><th>Count</th></tr></thead><tbody>{yrows}</tbody></table></div></details>
</section>"""


def build_sections(num):
    y_block = y_html()
    g, gpcs = load_pca("global")
    e, epcs = load_pca("eur")
    e = e[
        (e["super"] == "EUR") | (e["IID"] == "SELF")
    ]  # everyone was scored; show the EUR axes for EUR only
    g_near = nearest(g, gpcs, "super")
    e_near = nearest(e, epcs, "pop")
    cols, me, ref_mean = load_admixture()
    se_fp = RES / "ancestry" / "admixture_Q_se.tsv"
    se_row = (
        [l.split() for l in se_fp.read_text().splitlines() if l.startswith("SELF")]
        if se_fp.exists()
        else []
    )
    se = dict(zip(cols, map(float, se_row[0][1:]))) if se_row else {}
    se_of = lambda c: se.get(c, float("nan"))
    vals = [float(me[c]) for c in cols]
    top = cols[int(np.argmax(vals))]

    def near_table(df, namefn):
        rows = "".join(
            f"<tr><td>{esc(r.group)}</td><td>{esc(namefn(r.group))}</td><td>{r.dist:.4f}</td><td>{r.spread:.4f}</td></tr>"
            for r in df.itertuples()
        )
        return (
            '<div class="tablewrap"><table class="num"><thead><tr><th>Code</th><th>Population</th>'
            "<th>Distance from you (PC1–4)</th><th>Typical spread within population</th></tr></thead>"
            f"<tbody>{rows}</tbody></table></div>"
        )

    ref_rows = "".join(
        "<tr><td>"
        + esc(sp)
        + "</td>"
        + "".join(f"<td>{ref_mean.loc[sp, c] * 100:.1f}%</td>" for c in cols)
        + "</tr>"
        for sp in ref_mean.index
    )
    ref_table = (
        '<div class="tablewrap"><table class="num"><thead><tr><th>Reference group</th>'
        + "".join(f"<th>{esc(c)}</th>" for c in cols)
        + f"</tr></thead><tbody>{ref_rows}</tbody></table></div>"
    )
    n_snps = int(g.loc[g["IID"] == "SELF", "ALLELE_CT"].iloc[0]) // 2
    e_cent = e[e["IID"] != "SELF"].groupby("pop")[epcs].mean()
    e_me = e[e["IID"] == "SELF"].iloc[0]
    e_ratio = e_near.iloc[0]["dist"] / e_near.iloc[0]["spread"]
    kg_m = json.loads((RES / "membership" / "kg_eur.json").read_text())

    html_out = f"""
{y_block}

<section id="kg">
  <h2>5 · Ancestry compared with 1000 Genomes</h2>
  <p class="plainline"><b>In plain words:</b> PCA arranges people on a map according to how similar their DNA is. People from the same part of
  the world cluster together. We draw {
        num(len(g) - 1)
    } unrelated reference people from 26 world populations, then place you on the same map.
  ADMIXTURE then estimates what share of your DNA looks most like each continental group.</p>
  <aside class="caveat" role="note"><b>Limitation that matters here.</b> The 1000 Genomes European panel contains only five populations:
  Utah residents of Northern/Western European ancestry (CEU), Tuscans (TSI), British (GBR), Finns (FIN) and Iberians (IBS).
  Anyone whose ancestry lies elsewhere in Europe will land "between" the available groups or next to whichever one happens to be
  closest. That tells you where the reference samples are, not where your ancestors lived; section 6 uses a much larger reference.</aside>
  <p>{
        num(n_snps)
    } LD-pruned SNPs shared between your chip and the reference. Principal components were computed on the reference only;
  everyone, including you, was then projected with the same weights.</p>

  <h3>World map (PC1 vs PC2)</h3>
  <p class="note">Grey dots are reference individuals; labels mark the centre of each continental group. Nearest group: <b>{
        esc(g_near.iloc[0]["group"])
    }</b>.</p>
  {scatter(g, gpcs, "super", lambda s: SUPERPOP_NAMES.get(s, s))}
  {near_table(g_near, lambda s: SUPERPOP_NAMES.get(s, s))}

  <h3>Within Europe (PC1 vs PC2, European reference only)</h3>
  <p class="note">The same idea, zoomed in to the five European populations. In this data PC1 separates Finns
  (centre {
        e_cent.loc[
            "FIN", epcs[0]
        ]:+.3f}) from everyone else, and PC2 separates the north-west (CEU {
        e_cent.loc["CEU", epcs[1]]:+.3f},
  GBR {e_cent.loc["GBR", epcs[1]]:+.3f}) from the south (IBS {
        e_cent.loc["IBS", epcs[1]]:+.3f}, TSI {e_cent.loc["TSI", epcs[1]]:+.3f}).
  On a PCA computed with you included (no projection bias), the most similar of these five is <b>{
        esc(kg_m["self_top"][0]["pop"])
    }</b>
  ({esc(POP_NAMES.get(kg_m["self_top"][0]["pop"], ""))}). A test of fit gives p = {
        kg_m["fit_p_best"]:.3f} (a small p means "none of these"):
  {
        "you are not typical of any of these five, which is what an ancestry the panel doesn't cover looks like"
        if kg_m["fit_p_best"] < 0.05
        else "you fit within the range of that population"
    }.</p>
  {scatter(e, epcs, "pop", lambda s: POP_NAMES.get(s, s))}

  <h3>Continental proportions (supervised ADMIXTURE)</h3>
  <p class="note">Each reference person is fixed to their own continental group: 20 random people from each of that continent's populations,
  100 per group. The Admixed American group and the admixed African-Caribbean (ACB) and African-American (ASW) samples are left out because they
  are themselves mixtures. Largest component: <b>{
        esc(SUPERPOP_NAMES.get(top, top))
    }</b>, {vals[cols.index(top)] * 100:.1f}%
  (bootstrap SE {
        se_of(top)
        * 100:.1f} points). Components below about 1–2% are within the method's noise.
  An earlier version took the first 100 reference people by ID, which were mostly African-Caribbean; it reported 1.9% "African".</p>
  {stack_bar(cols, vals)}
  <details><summary>Sanity check: the same model on the reference groups</summary>{
        ref_table
    }</details>
</section>
"""
    return html_out, {
        "y": load_y()["hg_snp"]
        if (RES / "yhaplo" / "haplogroups.y.txt").exists()
        else None,  # only read for males
        "global": g_near.iloc[0]["group"],
        "eur": e_near.iloc[0]["group"],
        "adm": dict(zip(cols, vals)),
    }
