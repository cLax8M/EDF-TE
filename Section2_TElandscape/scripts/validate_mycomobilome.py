#!/usr/bin/env python3
"""validate_mycomobilome.py — Validate custom TE library annotations.

Compares EarlGrey annotations using MycoMobilome custom library against
published literature values for two well-studied reference genomes:
  - Cryptococcus neoformans H99
  - Candida albicans SC5314

Outputs:
  - Summary table (TSV)
  - TE composition comparison bar chart (PDF)
  - Superfamily-level breakdown (PDF)

Literature references:
  - Cneo H99: Loftus et al. 2005 (Science), Janbon et al. 2014 (Genome Biol)
    Standard RepeatMasker: ~5.0% total repeat, ~3.6% TIR
  - Calb SC5314: Jones et al. 2004 (PNAS), Goodwin et al. 2015 (Genome Res)
    Standard RepeatMasker: ~3.1% total repeat
"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Paths ──────────────────────────────────────────────────────────────────
ANNOTATED = Path("/data/carlos/EDF-TE/TE_analysis/annotated")
OUTDIR = Path("/legserv/Temp/Carlos/EDF-TE/results/Validation_MycoMobilome")

SPECIES = {
    "VL_CneoH99": {
        "label": "C. neoformans H99",
        "dir": ANNOTATED / "VL_CneoH99_EarlGrey",
        "summary_dir": "VL_CneoH99_summaryFiles",
        "prefix": "VL_CneoH99",
        # Published values (standard RepeatMasker / Dfam-only)
        "lit_te_pct": 5.0,         # Loftus et al. 2005 / Janbon et al. 2014
        "lit_source": "Janbon et al. 2014",
        "lit_ltr_pct": 3.6,        # Predominantly Ty3/gypsy
        "lit_dna_pct": 0.3,
        "lit_line_pct": 0.2,
        "lit_unclass_pct": 0.9,
    },
    "VL_CalbSC5314": {
        "label": "C. albicans SC5314",
        "dir": ANNOTATED / "VL_CalbSC5314_EarlGrey",
        "summary_dir": "VL_CalbSC5314_summaryFiles",
        "prefix": "VL_CalbSC5314",
        "lit_te_pct": 3.1,         # Jones et al. 2004 / Goodwin et al. 2015
        "lit_source": "Goodwin et al. 2015",
        "lit_ltr_pct": 0.9,
        "lit_dna_pct": 0.5,
        "lit_line_pct": 0.3,
        "lit_unclass_pct": 1.4,
    },
}

TE_CLASSES = ["DNA", "Rolling Circle", "Penelope", "LINE", "SINE", "LTR", "Unclassified"]
TE_COLORS = {
    "DNA": "#377eb8", "Rolling Circle": "#ff7f00", "Penelope": "#66c2a5",
    "LINE": "#4daf4a", "SINE": "#f1c232", "LTR": "#e41a1c",
    "Unclassified": "#bdbdbd", "Other": "#999999",
}


def parse_high_level(path):
    """Parse highLevelCount.txt → dict of class: {bp, copies, pct, families}."""
    data = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            cls = parts[0]
            if cls.endswith("-nested") or cls.startswith("Total") or cls.startswith("Non-Repeat"):
                continue
            if cls.startswith("Other"):
                cls = "Other"
            bp = int(parts[1])
            copies = int(parts[2])
            pct = float(parts[3])
            genome = int(parts[4])
            fam = parts[5]
            fam_count = int(fam) if fam != "NA" else 0
            data[cls] = {"bp": bp, "copies": copies, "pct": pct,
                         "genome": genome, "families": fam_count}
    return data


def parse_divergence_table(path):
    """Parse divergence_summary_table.tsv → DataFrame."""
    return pd.read_csv(path, sep="\t")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows = []

    for sp_id, info in SPECIES.items():
        print(f"\n{'='*60}")
        print(f"  {info['label']} ({sp_id})")
        print(f"{'='*60}")

        hl_path = info["dir"] / info["summary_dir"] / f"{info['prefix']}.highLevelCount.txt"
        div_path = info["dir"] / info["summary_dir"] / f"{info['prefix']}_divergence_summary_table.tsv"

        hl = parse_high_level(hl_path)
        div_df = parse_divergence_table(div_path)

        # Total TE content
        total_bp = sum(v["bp"] for v in hl.values())
        genome_size = list(hl.values())[0]["genome"]
        total_pct = total_bp / genome_size * 100
        total_copies = sum(v["copies"] for v in hl.values())
        total_families = sum(v["families"] for v in hl.values())

        lit_pct = info["lit_te_pct"]
        improvement = total_pct - lit_pct
        improvement_rel = (improvement / lit_pct) * 100

        print(f"  Genome size: {genome_size:,} bp")
        print(f"  EarlGrey+MycoMobilome: {total_pct:.2f}% TE ({total_bp:,} bp)")
        print(f"  Literature (standard): {lit_pct:.1f}%")
        print(f"  Improvement: +{improvement:.2f} pp ({improvement_rel:+.1f}%)")
        print(f"  Total copies: {total_copies:,}")
        print(f"  Total families: {total_families:,}")

        # Per-class breakdown
        print(f"\n  {'Class':<20} {'EG+MM %':>8} {'Lit %':>8} {'Δ':>8} {'Copies':>8} {'Fam':>6}")
        print(f"  {'-'*58}")
        for cls in TE_CLASSES:
            eg_pct = hl.get(cls, {}).get("pct", 0)
            lit_key = f"lit_{cls.lower().replace(' ', '_')}_pct"
            lit_val = info.get(lit_key, 0)
            delta = eg_pct - lit_val
            copies = hl.get(cls, {}).get("copies", 0)
            fam = hl.get(cls, {}).get("families", 0)
            print(f"  {cls:<20} {eg_pct:>8.3f} {lit_val:>8.3f} {delta:>+8.3f} {copies:>8,} {fam:>6,}")

        # Superfamily-level summary
        if not div_df.empty:
            n_superfam = div_df.groupby(["subclass", "superfamily"]).ngroups
            n_families = len(div_df["family"].unique())
            print(f"\n  Superfamilies detected: {n_superfam}")
            print(f"  Unique families: {n_families}")

            # Top superfamilies by total bp
            sf_bp = div_df.groupby(["subclass", "superfamily"])["total_bp"].sum()
            sf_bp = sf_bp.sort_values(ascending=False)
            print(f"\n  Top 10 superfamilies by coverage:")
            for (sub, sf), bp in sf_bp.head(10).items():
                print(f"    {sub}/{sf}: {bp:,} bp")

        # Store for table
        rows.append({
            "species": info["label"],
            "genome_Mb": genome_size / 1e6,
            "EG_MM_pct": round(total_pct, 2),
            "lit_pct": lit_pct,
            "lit_source": info["lit_source"],
            "delta_pp": round(improvement, 2),
            "delta_rel_pct": round(improvement_rel, 1),
            "total_copies": total_copies,
            "total_families": total_families,
            "DNA_pct": round(hl.get("DNA", {}).get("pct", 0), 3),
            "LTR_pct": round(hl.get("LTR", {}).get("pct", 0), 3),
            "LINE_pct": round(hl.get("LINE", {}).get("pct", 0), 3),
            "SINE_pct": round(hl.get("SINE", {}).get("pct", 0), 3),
            "Unclassified_pct": round(hl.get("Unclassified", {}).get("pct", 0), 3),
        })

    # ── Save summary table ────────────────────────────────────────────────
    df_out = pd.DataFrame(rows)
    tsv_path = OUTDIR / "validation_summary.tsv"
    df_out.to_csv(tsv_path, sep="\t", index=False)
    print(f"\n[SAVED] {tsv_path}")

    # ── Figure 1: Comparison bar chart ────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)

    for i, (sp_id, info) in enumerate(SPECIES.items()):
        ax = axes[i]
        hl_path = info["dir"] / info["summary_dir"] / f"{info['prefix']}.highLevelCount.txt"
        hl = parse_high_level(hl_path)

        classes = ["LTR", "DNA", "LINE", "SINE", "Unclassified"]
        eg_vals = [hl.get(c, {}).get("pct", 0) for c in classes]
        lit_vals = [info.get(f"lit_{c.lower()}_pct", 0) for c in classes]

        x = np.arange(len(classes))
        w = 0.35
        bars1 = ax.bar(x - w/2, eg_vals, w, label="EarlGrey + MycoMobilome",
                       color=[TE_COLORS.get(c, "#999") for c in classes],
                       edgecolor="black", linewidth=0.5, alpha=0.9)
        bars2 = ax.bar(x + w/2, lit_vals, w, label="Literature (standard RM)",
                       color=[TE_COLORS.get(c, "#999") for c in classes],
                       edgecolor="black", linewidth=0.5, alpha=0.35,
                       hatch="//")

        ax.set_xticks(x, classes, fontsize=8)
        ax.set_ylabel("% genome", fontsize=9)
        ax.set_title(info["label"], fontsize=10, fontweight="bold")
        ax.legend(fontsize=7, loc="upper right")

        # Total annotation
        total_eg = sum(eg_vals)
        total_lit = sum(lit_vals)
        ax.text(0.98, 0.85, f"Total TE: {total_eg:.1f}% vs {total_lit:.1f}%\n"
                f"Δ = +{total_eg-total_lit:.1f} pp",
                transform=ax.transAxes, ha="right", fontsize=7,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow",
                         edgecolor="gray", alpha=0.8))

    fig.suptitle("TE annotation: EarlGrey+MycoMobilome vs published (standard RM)",
                 fontsize=11, y=1.02)
    fig.tight_layout()
    out1 = OUTDIR / "validation_comparison.pdf"
    fig.savefig(out1, dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "validation_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[SAVED] {out1}")

    # ── Figure 2: Superfamily breakdown (EarlGrey+MM only) ───────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for i, (sp_id, info) in enumerate(SPECIES.items()):
        ax = axes[i]
        div_path = info["dir"] / info["summary_dir"] / f"{info['prefix']}_divergence_summary_table.tsv"
        div_df = parse_divergence_table(div_path)

        sf_bp = div_df.groupby("subclass")["total_bp"].sum().sort_values(ascending=True)
        colors = [TE_COLORS.get(c, "#999") for c in sf_bp.index]
        sf_bp_kb = sf_bp / 1000

        ax.barh(range(len(sf_bp)), sf_bp_kb.values, color=colors, edgecolor="none")
        ax.set_yticks(range(len(sf_bp)))
        ax.set_yticklabels(sf_bp.index, fontsize=8)
        ax.set_xlabel("Total coverage (kb)", fontsize=9)
        ax.set_title(info["label"], fontsize=10, fontweight="bold")

        # Annotate top 3
        for j, (cls, val) in enumerate(sf_bp.tail(3).items()):
            ax.text(val/1000 + sf_bp_kb.max()*0.02, len(sf_bp)-3+j,
                    f"{val/1000:.1f} kb", fontsize=7, va="center")

    fig.suptitle("TE superfamily coverage (EarlGrey + MycoMobilome)",
                 fontsize=11, y=1.02)
    fig.tight_layout()
    out2 = OUTDIR / "validation_superfamily_breakdown.pdf"
    fig.savefig(out2, dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "validation_superfamily_breakdown.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[SAVED] {out2}")

    # ── Figure 3: Kimura divergence distribution ─────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for i, (sp_id, info) in enumerate(SPECIES.items()):
        ax = axes[i]
        gff_path = info["dir"] / f"{info['prefix']}_summaryFiles" / f"{info['prefix']}.filteredRepeats.gff"

        kimura_vals = []
        classes_for_color = []
        with open(gff_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                te_class = parts[2]
                attrs = parts[8]
                # Extract kimura80 value
                for attr in attrs.split(";"):
                    if attr.startswith("kimura80="):
                        try:
                            k = float(attr.split("=")[1])
                            kimura_vals.append(k)
                            # Map to broad class
                            if "LTR" in te_class:
                                classes_for_color.append("LTR")
                            elif "DNA" in te_class or "TIR" in te_class:
                                classes_for_color.append("DNA")
                            elif "LINE" in te_class:
                                classes_for_color.append("LINE")
                            elif "SINE" in te_class:
                                classes_for_color.append("SINE")
                            elif "RC" in te_class or "Helitron" in te_class:
                                classes_for_color.append("Rolling Circle")
                            elif "Unknown" in te_class:
                                classes_for_color.append("Unclassified")
                            else:
                                classes_for_color.append("Other")
                        except ValueError:
                            pass
                        break

        if not kimura_vals:
            continue

        # Stacked histogram by class
        bins = np.arange(0, 0.55, 0.02)
        class_order = ["LTR", "DNA", "LINE", "SINE", "Rolling Circle", "Unclassified", "Other"]
        bottom = np.zeros(len(bins) - 1)
        for cls in class_order:
            vals = [k for k, c in zip(kimura_vals, classes_for_color) if c == cls]
            if vals:
                counts, _ = np.histogram(vals, bins=bins)
                ax.bar(bins[:-1], counts, width=0.018, bottom=bottom,
                       color=TE_COLORS.get(cls, "#999"), label=cls,
                       edgecolor="none", alpha=0.85)
                bottom += counts

        ax.set_xlabel("Kimura divergence", fontsize=9)
        ax.set_ylabel("Copy count", fontsize=9)
        ax.set_title(info["label"], fontsize=10, fontweight="bold")
        ax.legend(fontsize=6, loc="upper right")
        ax.set_xlim(0, 0.55)

    fig.suptitle("Kimura divergence distribution (EarlGrey + MycoMobilome)",
                 fontsize=11, y=1.02)
    fig.tight_layout()
    out3 = OUTDIR / "validation_kimura_distribution.pdf"
    fig.savefig(out3, dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "validation_kimura_distribution.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[SAVED] {out3}")

    print(f"\n{'='*60}")
    print(f"  Validation complete. All outputs in:\n  {OUTDIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
