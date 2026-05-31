#!/usr/bin/env python3
"""
Report: Key findings from TE landscape analysis
Sections: Kimura/Shannon diversity, Synteny & TE influence
"""

import io
from pathlib import Path
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib.colors import HexColor, white, black
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, PageBreak, Table, TableStyle
)
from reportlab.platypus.flowables import Flowable
from pypdf import PdfWriter, PdfReader

# ── Paths ─────────────────────────────────────────────────────────────────────
FIG_DIR   = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figures")
SYN_DIR   = Path("/legserv/Temp/Carlos/EDF-TE/results/Synteny/figures_filtered_test")
OUT_DIR   = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape")
OUT_PDF   = OUT_DIR / "TE_findings_report.pdf"

# ── Colours ───────────────────────────────────────────────────────────────────
DARK      = HexColor("#1B2A3B")
MID       = HexColor("#2E6DA4")
LIGHT     = HexColor("#D6E8F7")
ACCENT    = HexColor("#E86E2D")

W, H = A4         # 595.27 x 841.89 pts

# ── Styles ────────────────────────────────────────────────────────────────────
styles = getSampleStyleSheet()

def make_styles():
    s = {}
    s["cover_title"] = ParagraphStyle(
        "cover_title", fontName="Helvetica-Bold", fontSize=26, leading=32,
        textColor=white, spaceAfter=12, alignment=1,
    )
    s["cover_sub"] = ParagraphStyle(
        "cover_sub", fontName="Helvetica", fontSize=14, leading=18,
        textColor=LIGHT, spaceAfter=8, alignment=1,
    )
    s["sec_title"] = ParagraphStyle(
        "sec_title", fontName="Helvetica-Bold", fontSize=18, leading=22,
        textColor=white, spaceAfter=6, alignment=0,
    )
    s["body"] = ParagraphStyle(
        "body", fontName="Helvetica", fontSize=10, leading=14,
        textColor=DARK, spaceAfter=6,
    )
    s["caption"] = ParagraphStyle(
        "caption", fontName="Helvetica-Oblique", fontSize=8.5, leading=11,
        textColor=HexColor("#555555"), spaceAfter=4, alignment=1,
    )
    s["bullet"] = ParagraphStyle(
        "bullet", fontName="Helvetica", fontSize=10, leading=14,
        textColor=DARK, spaceAfter=4, leftIndent=12, bulletIndent=0,
    )
    return s

ST = make_styles()

# ── Helper: coloured banner flowable ─────────────────────────────────────────
class Banner(Flowable):
    def __init__(self, text, color=MID, height=1.1*cm, text_style=None):
        Flowable.__init__(self)
        self.text  = text
        self.color = color
        self.height = height
        self.text_style = text_style or ST["sec_title"]
        self.width = W - 4*cm   # will be overridden by wrap

    def wrap(self, avail_w, avail_h):
        self.width = avail_w
        return avail_w, self.height

    def draw(self):
        self.canv.setFillColor(self.color)
        self.canv.roundRect(0, 0, self.width, self.height, 6, fill=1, stroke=0)
        self.canv.setFillColor(white)
        self.canv.setFont(self.text_style.fontName, self.text_style.fontSize)
        self.canv.drawString(10, (self.height - self.text_style.fontSize) / 2 + 2, self.text)


# ── Cover page (reportlab) ────────────────────────────────────────────────────
def build_cover_pdf() -> bytes:
    buf = io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=A4,
                            leftMargin=2.5*cm, rightMargin=2.5*cm,
                            topMargin=2*cm, bottomMargin=2*cm)

    def draw_bg(canvas, doc):
        canvas.setFillColor(DARK)
        canvas.rect(0, 0, W, H, fill=1, stroke=0)
        canvas.setFillColor(MID)
        canvas.rect(0, H*0.58, W, H*0.42, fill=1, stroke=0)
        # accent stripe
        canvas.setFillColor(ACCENT)
        canvas.rect(0, H*0.57, W, 6, fill=1, stroke=0)

    story = [
        Spacer(1, 2.8*cm),
        Paragraph("TE Landscape Analysis", ST["cover_title"]),
        Paragraph("Key Findings Report", ST["cover_sub"]),
        Spacer(1, 0.3*cm),
        Paragraph("Early Diverging Fungi — 238 species", ST["cover_sub"]),
        Spacer(1, 3.5*cm),
    ]

    # Sections box
    sections = [
        ["Section 1", "Kimura distance & TE activity"],
        ["Section 2", "Shannon diversity of TE superfamilies"],
        ["Section 3", "Microsynteny & TE influence"],
    ]
    tbl = Table(sections, colWidths=[3.5*cm, 9*cm])
    tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (0, -1), ACCENT),
        ("BACKGROUND", (1, 0), (1, -1), HexColor("#203245")),
        ("TEXTCOLOR", (0, 0), (-1, -1), white),
        ("FONTNAME",  (0, 0), (0, -1), "Helvetica-Bold"),
        ("FONTNAME",  (1, 0), (1, -1), "Helvetica"),
        ("FONTSIZE",  (0, 0), (-1, -1), 11),
        ("ROWBACKGROUNDS", (0, 0), (-1, -1), [HexColor("#203245"), HexColor("#25394E")]),
        ("GRID",      (0, 0), (-1, -1), 0.5, HexColor("#3A5570")),
        ("LEFTPADDING",  (0,0), (-1,-1), 10),
        ("RIGHTPADDING", (0,0), (-1,-1), 10),
        ("TOPPADDING",   (0,0), (-1,-1), 8),
        ("BOTTOMPADDING",(0,0), (-1,-1), 8),
        ("VALIGN",    (0, 0), (-1, -1), "MIDDLE"),
    ]))
    story.append(tbl)
    story.append(Spacer(1, 1.5*cm))
    story.append(Paragraph("April 2026", ST["cover_sub"]))

    doc.build(story, onFirstPage=draw_bg, onLaterPages=draw_bg)
    return buf.getvalue()


# ── Section intro page ────────────────────────────────────────────────────────
def build_section_page(title: str, subtitle: str, bullets: list[str],
                       color=MID) -> bytes:
    buf = io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=A4,
                            leftMargin=2.5*cm, rightMargin=2.5*cm,
                            topMargin=2.5*cm, bottomMargin=2*cm)

    def draw_stripe(canvas, doc):
        canvas.setFillColor(color)
        canvas.rect(0, H - 4.5*cm, W, 4.5*cm, fill=1, stroke=0)
        canvas.setFillColor(ACCENT)
        canvas.rect(0, H - 4.5*cm - 4, W, 4, fill=1, stroke=0)

    story = [
        Spacer(1, 1.8*cm),
        Paragraph(title, ST["cover_title"]),
        Paragraph(subtitle, ParagraphStyle(
            "ss", fontName="Helvetica-Oblique", fontSize=12, leading=16,
            textColor=LIGHT, alignment=1)),
        Spacer(1, 2.2*cm),
    ]
    for b in bullets:
        story.append(Paragraph(f"▸  {b}", ST["bullet"]))
        story.append(Spacer(1, 0.15*cm))

    # white background for text area
    def draw_bg(canvas, doc):
        draw_stripe(canvas, doc)
        canvas.setFillColor(white)
        canvas.rect(0, 0, W, H - 4.5*cm - 4, fill=1, stroke=0)

    doc.build(story, onFirstPage=draw_bg, onLaterPages=draw_bg)
    return buf.getvalue()


# ── Assemble: merge PDFs ───────────────────────────────────────────────────────
def append_pdf(writer: PdfWriter, pdf_bytes: bytes):
    reader = PdfReader(io.BytesIO(pdf_bytes))
    for page in reader.pages:
        writer.add_page(page)

def append_pdf_file(writer: PdfWriter, path: Path):
    if not path.exists():
        print(f"  [skip] not found: {path}")
        return
    reader = PdfReader(str(path))
    for page in reader.pages:
        writer.add_page(page)
    print(f"  [ok] {path.name}")


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    writer = PdfWriter()

    # ── COVER ─────────────────────────────────────────────────────────────────
    print("Building cover...")
    append_pdf(writer, build_cover_pdf())

    # ══════════════════════════════════════════════════════════════════════════
    # SECTION 1 — Kimura distance & TE activity
    # ══════════════════════════════════════════════════════════════════════════
    print("\nSection 1: Kimura / TE age...")
    append_pdf(writer, build_section_page(
        title    = "Section 1 · Kimura Distance & TE Activity",
        subtitle = "TE age distribution and recent burst dynamics across 238 fungal/holozoan species",
        bullets  = [
            "Kimura distance (CpG-corrected) between TE copies and their consensus is used as a "
            "proxy for insertion age: lower distance = more recent activity.",
            "The Recency Index (RI) summarises the fraction of TE content at Kimura < 5 "
            "(recently active copies).",
            "Most Mucoromycota (MU) and Chytridiomycota (CH) show elevated RI, indicating "
            "ongoing or recent TE bursts — especially LTR retrotransposons.",
            "Ascomycota (AS) and Basidiomycota (BA) show older, more diverged TE landscapes "
            "consistent with stronger silencing (RIP/methylation).",
            "A positive correlation between RI and Shannon diversity suggests that species "
            "experiencing recent bursts also harbour more diverse TE superfamilies.",
        ],
        color = MID,
    ))

    kimura_figs = [
        (FIG_DIR / "kimura_landscape_stacked.pdf",
         "Figure 1.1 — Kimura distance landscape (stacked) across all 238 species. "
         "Each column is a species; colours represent TE classes. "
         "Peaks at low Kimura values indicate recent TE activity."),
        (FIG_DIR / "kimura_clade_heatmap.pdf",
         "Figure 1.2 — Heatmap of mean Kimura distance per TE superfamily across clades. "
         "Darker colours indicate more recently active families."),
        (FIG_DIR / "kimura_ri_by_clade_boxplot.pdf",
         "Figure 1.3 — Recency Index (RI) distribution per clade. "
         "MU and CH clades show the highest recent TE activity."),
        (FIG_DIR / "kimura_ri_vs_shannon.pdf",
         "Figure 1.4 — Correlation between RI and Shannon diversity (H'). "
         "Species with more recent TE bursts tend to show higher superfamily diversity."),
        (FIG_DIR / "kimura_ri_vs_mean_div.pdf",
         "Figure 1.5 — RI vs. mean Kimura divergence. "
         "The inverse trend confirms that high-RI species have globally younger TE insertions."),
    ]

    for fig_path, caption in kimura_figs:
        append_pdf_file(writer, fig_path)

    # ══════════════════════════════════════════════════════════════════════════
    # SECTION 2 — Shannon diversity
    # ══════════════════════════════════════════════════════════════════════════
    print("\nSection 2: Shannon diversity...")
    append_pdf(writer, build_section_page(
        title    = "Section 2 · Shannon Diversity of TE Superfamilies",
        subtitle = "Evenness and richness of TE composition across the fungal tree of life",
        bullets  = [
            "Shannon entropy H' measures the diversity of TE superfamilies by abundance; "
            "Pielou's J measures evenness (how uniformly copies are distributed).",
            "Higher H' in early-diverging lineages (CH, MU) reflects a broader and more "
            "even repertoire of TE superfamilies.",
            "Ascomycota show low H' driven by dominance of a single class (often "
            "Unclassified or DNA transposons after RIP-mediated degeneration).",
            "A strong positive correlation exists between total TE content (%) and H', "
            "suggesting that larger TE loads are accompanied by greater diversity.",
            "The dominant superfamily per species (barplot) reveals clade-level signatures: "
            "LTR/Gypsy in MU, DNA/TIR in CH, LINE-dominant in early holozoans.",
        ],
        color = HexColor("#1A6B45"),
    ))

    shannon_figs = [
        (FIG_DIR / "shannon_h_by_clade_boxplot.pdf",
         "Figure 2.1 — Shannon H' per clade. "
         "Mucoromycota and Chytridiomycota harbour the most diverse TE landscapes."),
        (FIG_DIR / "shannon_j_by_clade_boxplot.pdf",
         "Figure 2.2 — Pielou's J (evenness) per clade. "
         "Higher evenness in MU indicates multiple co-existing superfamilies rather than "
         "dominance by a single type."),
        (FIG_DIR / "shannon_h_vs_total_te.pdf",
         "Figure 2.3 — Shannon H' vs. total TE content (%). "
         "Species with higher TE burden tend to harbour more diverse superfamilies."),
        (FIG_DIR / "shannon_h_vs_n_superfamilies.pdf",
         "Figure 2.4 — H' vs. number of detected superfamilies. "
         "Confirms that diversity increases with superfamily richness, though evenness varies."),
        (FIG_DIR / "shannon_dominant_sf_barplot.pdf",
         "Figure 2.5 — Dominant TE superfamily per species. "
         "Bar colours indicate the most abundant superfamily; ordered by clade and total TE%."),
    ]

    for fig_path, caption in shannon_figs:
        append_pdf_file(writer, fig_path)

    # ══════════════════════════════════════════════════════════════════════════
    # SECTION 3 — Microsynteny & TE influence
    # ══════════════════════════════════════════════════════════════════════════
    print("\nSection 3: Synteny...")
    append_pdf(writer, build_section_page(
        title    = "Section 3 · Microsynteny & TE Influence",
        subtitle = "Conservation of gene neighbourhood and disruption by transposable elements",
        bullets  = [
            "Microsynteny was assessed by computing the fraction of shared gene neighbours "
            "(collinear gene pairs within a sliding window) across pairwise species comparisons.",
            "Within-clade synteny conservation is generally high in AS and BA, "
            "intermediate in CH, and lowest in MU — consistent with elevated TE activity.",
            "Species with high TE content (%) show significantly lower synteny scores, "
            "supporting the hypothesis that TE insertions drive genome rearrangements.",
            "Heatmaps show pairwise microsynteny conservation within each clade; "
            "diagonal blocks reveal subgroup structure matching the species phylogeny.",
            "MU shows the highest within-clade heterogeneity, with some species retaining "
            "near-zero synteny to close relatives — coinciding with the TE outliers "
            "(Jimfl, Jimlac1, Endsp1 with >65% TE).",
        ],
        color = HexColor("#7B2D8B"),
    ))

    synteny_figs = [
        (SYN_DIR / "MU_synteny_heatmap.pdf",
         "Figure 3.1 — Pairwise microsynteny heatmap for Mucoromycota (MU). "
         "TE-rich species cluster as outliers with reduced synteny to all relatives."),
        (SYN_DIR / "CH_synteny_heatmap.pdf",
         "Figure 3.2 — Pairwise microsynteny heatmap for Chytridiomycota (CH). "
         "Intermediate synteny conservation consistent with moderate TE loads."),
        (SYN_DIR / "BA_synteny_heatmap.pdf",
         "Figure 3.3 — Pairwise microsynteny heatmap for Basidiomycota (BA). "
         "High synteny consistent with strong TE silencing (RIP, heterochromatin)."),
        (SYN_DIR / "GL_synteny_heatmap.pdf",
         "Figure 3.4 — Pairwise microsynteny heatmap for Glomeromycota (GL). "
         "Despite large genomes, most species retain moderate synteny — reflecting "
         "ancient TE insertions rather than ongoing disruption."),
    ]

    for fig_path, caption in synteny_figs:
        append_pdf_file(writer, fig_path)

    # ── Write output ───────────────────────────────────────────────────────────
    print(f"\nWriting {OUT_PDF} ...")
    with open(OUT_PDF, "wb") as f:
        writer.write(f)
    print(f"Done — {OUT_PDF.stat().st_size / 1024:.0f} KB")


if __name__ == "__main__":
    main()
