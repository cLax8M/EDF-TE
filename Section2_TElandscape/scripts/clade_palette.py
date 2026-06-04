"""Canonical clade color palette and full-name labels.

Single source of truth used by all figures in the EDF-TE manuscript.
Colors match the legend approved by the user (April 2026):
18 clades including NE (Neocallimastigomycetes) and HO (Holozoa+Amoebozoa).
"""

CLADE_COLORS = {
    # Dikarya
    "AS": "#e8887a",  # Ascomycota              (salmon)
    "BA": "#4472c4",  # Basidiomycota           (steel blue)
    # Mucoromycota subphyla
    "MU": "#e07840",  # Mucoromycotina          (orange)
    "MO": "#f4a8b0",  # Mortierellomycotina     (light pink)
    "GL": "#8b5e48",  # Glomeromycotina         (brown)
    # Zoopagomycota subphyla
    "KX": "#f0a820",  # Kickxellomycotina       (amber)
    "ZO": "#c8dc70",  # Zoopagomycotina         (lime green)
    "EN": "#e01010",  # Entomophthoromycotina   (vivid red)
    "BS": "#d0b4e0",  # Basidiobolomycota       (lavender)
    # Early-diverging fungi
    "OL": "#3668b0",  # Olpidiomycota           (medium blue)
    "CH": "#e0d040",  # Chytridiomycetes        (muted yellow)
    "NE": "#9060b0",  # Neocallimastigomycetes  (medium purple)
    "BL": "#7030a0",  # Blastocladiomycota      (dark purple)
    "MI": "#3aab98",  # Microsporidia           (teal)
    "CR": "#a8c8e8",  # Cryptomycota            (powder blue)
    "AP": "#5c9e3a",  # Aphelidiomycota         (forest green)
    # Outgroups
    "NU": "#1a2850",  # Nucleariida             (dark navy)
    "HO": "#b4d0e8",  # Holozoa+Amoebozoa       (pale blue)
}

CLADE_NAMES = {
    "AS": "Ascomycota",
    "BA": "Basidiomycota",
    "MU": "Mucoromycotina",
    "MO": "Mortierellomycotina",
    "GL": "Glomeromycotina",
    "KX": "Kickxellomycotina",
    "ZO": "Zoopagomycotina",
    "EN": "Entomophthoromycotina",
    "BS": "Basidiobolomycota",
    "OL": "Olpidiomycota",
    "CH": "Chytridiomycetes",
    "NE": "Neocallimastigomycetes",
    "BL": "Blastocladiomycota",
    "MI": "Microsporidia",
    "CR": "Cryptomycota",
    "AP": "Aphelidiomycota",
    "NU": "Nucleariida",
    "HO": "Holozoa+Amoebozoa",
}

# Phylogenetic order (root → Dikarya), used for boxplots / heatmaps
CLADE_ORDER = [
    "HO", "NU", "AP", "CR", "MI", "BL", "OL", "CH", "NE",
    "BS", "EN", "ZO", "KX", "GL", "MO", "MU", "BA", "AS",
]

# Neocallimastigomycetes species are stored with a CH_ prefix in EarlGrey but
# phylogenetically belong to the NE clade. Override here.
NE_SPECIES = {
    "CH_Anasp1", "CH_Neocon1", "CH_Neolan1",
    "CH_Orpsp1", "CH_Pecora1", "CH_Caecom1", "CH_Pirfi3",
}


def get_clade(species_name: str) -> str:
    """Return the canonical 2-letter clade code for a species name,
    applying NE_SPECIES override."""
    if species_name in NE_SPECIES:
        return "NE"
    return species_name[:2]


# Canonical TE class palette (matches cladogram_te legend)
TE_CLASS_COLORS = {
    "LTR":            "#e41a1c",  # red
    "DNA":            "#377eb8",  # blue
    "LINE":           "#4daf4a",  # green
    "SINE":           "#f1c232",  # gold
    "Rolling Circle": "#ff7f00",  # orange
    "Penelope":       "#66c2a5",  # teal
    "PLE":            "#66c2a5",  # teal (alias)
    "Unclassified":   "#bdbdbd",  # grey
    "Unknown":        "#bdbdbd",
    "Other":          "#bdbdbd",
    "Total TE":       "#000000",  # black for total row
}


# Species display names: short binomial in italics, e.g. "MU_Absrep1" → "A. repens".
# Centralised here so figures can use a single map.
SPECIES_DISPLAY = {
    # Mucoromycotina
    "MU_Abscae1":   "A. caerulea",   "MU_Absrep1":   "A. repens",
    "MU_Chabre1":   "C. brevicaulis","MU_Chocucu1":  "C. cucurbitarum",
    "MU_Cunech1":   "C. echinulata", "MU_Cirumb1":   "C. umbellata",
    "MU_Cokrec1":   "C. recurvatus", "MU_Ellano2":   "E. lanuginosa",
    "MU_Fenlin":    "F. linderi",    "MU_Gilper1":   "G. persicaria",
    "MU_Hesve2":    "H. vesiculosa", "MU_Halrad1":   "H. radiata",
    "MU_Helpul1":   "H. pulchra",    "MU_Kircor1":   "K. cordense",
    "MU_Lenpar1":   "L. parva",      "MU_Lenzyc1":   "L. zycha",
    "MU_Mucam1":    "M. ambiguus",   "MU_Mucci3":    "M. circinelloides",
    "MU_Muccirc1":  "M. circinelloides","MU_Mucend1": "M. endophyticus",
    "MU_Mucmuc1":   "M. mucedo",     "MU_Mucpir1":   "M. piriformis",
    "MU_Mycafr1":   "M. africanus",  "MU_Mycind1":   "M. indica",
    "MU_Phaart1":   "P. articulata", "MU_Phybl4":    "P. blakesleeanus",
    "MU_PhyblU21":  "P. blakesleeanus U21",
    "MU_Pilano1":   "P. anomala",    "MU_Pilumb1":   "P. umbellatus",
    "MU_PnitS607":  "P. nitens S607","MU_PnitSJTU":  "P. nitens SJTU",
    "MU_Rhihy3":    "R. hyalina",    "MU_Rhimic1":   "R. micropora",
    "MU_Rhirh1":    "R. rhizopodiformis","MU_Rhiar1": "R. arrhizus",
    # Mortierellomycotina
    "MO_Moramb1":   "M. ambigua",    "MO_Morwol1":   "M. wolfii",
    # Microsporidia
    "MI_Enche1":    "E. cuniculi",   "MI_Encro1":    "E. romaleae",
    # Neocallimastigomycetes (CH-prefixed in our pipeline)
    "CH_Anasp1":    "A. robustus",
    "CH_Caecom1":   "C. churrovis",
    "CH_Neocon1":   "N. constans",
    "CH_Neolan1":   "N. lanati",
    "CH_Orpsp1":    "Orpinomyces sp.",
    "CH_Pecora1":   "Pecoramyces sp.",
    "CH_Pirfi3":    "P. finnis",
    # Outgroup (yeasts as VL/control) — EXCLUDED from all main figures
    # (kept here only so that legacy code does not crash on lookup)
    "VL_CalbSC5314": "C. albicans",
    "VL_CneoH99":    "C. neoformans",
    # ── User-supplied curated names (2026-04-21) ─────────────────────────
    "BA_Maglo":     "M. globosa",
    "BS_Basme":     "B. meristosporus",
    "CH_Batde5":    "B. dendrobatidis",
    "CH_Caupr1":    "C. protostelioides",
    "CR_Rozal1":    "R. allomycis",
    "EN_Conth1":    "C. thromboides",
    "HO_Drome":     "D. melanogaster",
    "MI_MicG_I_3":  "O. pajunii",
    "MO_MoalAD071": "M. alpina AD071",
    "MO_MoamBC1065":"M. amoeboidea",
    "MO_MochAD033": "M. chlamydospora",
    "MO_MoclAM1000":"M. claussenii",
    "MO_MoelGBA34": "M. elongata",
    "MO_MozyN2592": "M. zychae",
    "MU_Jimfl":     "J. flammicorona",
    "MU_Rhimi59":   "R. microsporus ATCC 11559",
    "MU_Rhior3":    "R. delemar",
    "ZO_Pipcy3":    "P. cylindrospora",
    # Disambiguations (the auto-parser truncates these too aggressively)
    "AS_AspfumUD":  "A. fumigatus",
    "AS_AspfumUD1": "A. fumigatus",
    "CH_Entlut1":   "E. luteoviride",
    "MU_Endsp1":    "Endogone sp.",
    "MU_Rhisto1":   "R. stolonifer",
    # Lunasporangiospora selenospora (auto-parser falls back to "M. sel"
    # because the MycoCosm short ID column carries an NCBI accession)
    "MO_Morsel1":   "L. selenospora",
    # ── User-supplied curated names (2026-04-22) ─────────────────────────
    "AP_Amoocci":    "A. occidentale",
    "AP_Amoprot":    "A. protococcarum",
    "BL_Allomarbu":  "A. arbusculus",
    "CH_Dinkinn":    "D. kinnereticum",
    "CH_Germichi":   "G. michiganensis",
    "CH_Lobangu":    "L. angularis",
    "CH_Phybulla":   "P. bullatum",
    "CH_Rhirose":    "R. rosea",
    "CH_Terramyc":   "Terramyces sp.",
    "HO_Acancaste":  "A. castellanii",
    "HO_Capowcza":   "C. owczarzaki",
    "HO_Ciointe":    "C. intestinalis",
    "HO_Dicdisco":   "D. discoideum",
    "HO_Homosap":    "H. sapiens",
    "HO_Hydravulg":  "H. vulgaris",
    "HO_Magisepte":  "M. septendecula",
    "HO_Mneleidy":   "M. leidyi",
    "HO_Monbrevi":   "M. brevicollis",
    "HO_Naegrub":    "N. gruberi",
    "HO_Salpirose":  "S. rosetta",
    "HO_Sphaearti":  "S. arctica",
    "HO_Thecatra":   "T. trahens",
    "HO_Trichadha":  "T. adhaerens",
    "KX_Dispsimple": "D. simplex",
    "KX_Spiroaspi":  "S. aspiralis",
    "MU_Apoossi":    "A. ossiformis",
    "MU_Blatri":     "B. trispora",
    "MU_Chlapade":   "C. padenii",
    "NU_Fontalba":   "F. alba",
    "NU_Parvatl":    "P. atlantica",
    # Add more as needed; default fallback handles unknowns
}


# Species/clades to EXCLUDE globally from main figures (control yeasts).
# Filter helpers below honour this set; per-script consumers can also use it
# directly via `EXCLUDE_SPECIES`/`EXCLUDE_CLADES`.
EXCLUDE_CLADES  = {"VL"}
EXCLUDE_SPECIES = {
    "VL_CalbSC5314", "VL_CneoH99",
    # Bare short IDs (in case scripts use them without prefix)
    "CalbSC5314", "CneoH99",
}


def is_excluded(species_id: str) -> bool:
    """Return True if a species/clade should be hidden from main figures."""
    if species_id in EXCLUDE_SPECIES:
        return True
    if species_id[:2] in EXCLUDE_CLADES:
        return True
    return False


# ── Auto-load full Genus species names from Genome_List.txt ─────────────
# This is the master mapping table maintained by the user. Columns:
#   col 4 = clade code  (e.g. "MU")
#   col 5 = full organism name ("Absidia repens")
#   col 9 = genus,  col 10 = species
#   col 15 = MycoCosm short ID ("Absrep1") used in our species_id strings
#
# We load it once at import; values populate SPECIES_DISPLAY without
# overriding any existing curated entries.
def _load_genome_list():
    import os as _os
    path = "/legserv/Temp/Carlos/EDF-TE/data/Genome_List.txt"
    if not _os.path.isfile(path):
        return
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            header = fh.readline()  # skip header
            for ln in fh:
                cols = ln.rstrip("\n").split("\t")
                if len(cols) < 15:
                    continue
                clade  = cols[3].strip()
                organ  = cols[4].strip()
                genus  = cols[8].strip() if len(cols) > 8 else ""
                spec   = cols[9].strip() if len(cols) > 9 else ""
                short  = cols[14].strip() if len(cols) > 14 else ""
                # Strip trailing "_NNN" suffix some MycoCosm IDs carry
                # (e.g. "Orpsp1_1" → "Orpsp1", "AmpWSBS2006_1" → "AmpWSBS2006")
                import re as _re_local
                short = _re_local.sub(r'_\d+$', '', short)
                if not clade or not short:
                    continue
                # Build "G. species" abbreviation
                if genus and spec and spec.lower() not in ("sp", "sp.", "nan", ""):
                    pretty = f"{genus[0].upper()}. {spec.lower()}"
                elif organ:
                    parts = organ.split()
                    if len(parts) >= 2:
                        pretty = f"{parts[0][0].upper()}. {parts[1].lower()}"
                    else:
                        pretty = organ
                else:
                    continue
                key = f"{clade}_{short}"
                # Don't overwrite curated entries
                if key not in SPECIES_DISPLAY:
                    SPECIES_DISPLAY[key] = pretty
                # Also accept the bare short ID without prefix
                if short not in SPECIES_DISPLAY:
                    SPECIES_DISPLAY[short] = pretty
    except Exception as exc:
        import sys as _sys
        print(f"[clade_palette] WARN: could not parse Genome_List.txt: {exc}",
              file=_sys.stderr)


_load_genome_list()


def species_display(species_id: str, italic: bool = True) -> str:
    """Return human-readable species name (italics-ready) for a given ID.

    Strategy:
      1. Lookup in SPECIES_DISPLAY (curated overrides).
      2. Otherwise apply MycoCosm convention: first 3 letters of genus +
         first 3 letters of species (e.g. "Absrep1" → "A. repens").
      3. Final fallback: return the suffix unchanged.
    """
    import re as _re
    if species_id in SPECIES_DISPLAY:
        nm = SPECIES_DISPLAY[species_id]
    else:
        body = species_id.split("_", 1)[-1] if "_" in species_id else species_id
        # Strip trailing version suffix: digits, "U21", "_c", "JTU", "S607" etc.
        base = _re.sub(r'(?:[A-Z]+\d*|_.*|\d+)$', '', body)
        if len(base) >= 6:
            genus_letters = base[:3]
            species_letters = base[3:].lower()
            # Take up to 6 chars of species portion
            species_letters = species_letters[:8]
            genus_initial = genus_letters[0].upper() + "."
            nm = f"{genus_initial} {species_letters}"
        else:
            nm = body
    if italic:
        # Escape characters that mathtext would otherwise interpret
        # specially (underscore → subscript, etc.).
        nm_tex = (nm.replace("\\", r"\backslash ")
                    .replace("_", r"\_")
                    .replace("^", r"\^{}")
                    .replace(" ", r"\ "))
        return f"$\\mathit{{{nm_tex}}}$"
    return nm
