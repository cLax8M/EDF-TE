#!/bin/bash
#SBATCH --job-name=replot_colors
#SBATCH --partition=normal.1000h
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/data/carlos/EDF-TE/logs/replot_colors_%j.log

eval "$(micromamba shell hook --shell bash)"
micromamba activate syntenyEDFTE

cd /data/carlos/EDF-TE/EDF-TE/Section2_TElandscape/scripts

# Re-run figF2 (in vivo, 4 metrics) with canonical colors
python - << 'EOF'
import sys; sys.path.insert(0,".")
import pandas as pd
from chromatin_by_class_multimetric import run_invivo, make_figure2
df_iv = run_invivo()
make_figure2(df_iv)

# Re-run figF1 (all fungi, 3 metrics) using saved TSV + new colors
from chromatin_by_class_multimetric import make_figure1
df_all = pd.read_csv(
    "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/invivo_landscape/figF_summary_allfungi.tsv",
    sep="\t")
make_figure1(df_all)
print("Done F1+F2")
EOF

# Re-run figG1 (all fungi NRL) using saved TSV + new colors
python - << 'EOF'
import sys; sys.path.insert(0,".")
import pandas as pd
from nrl_allfungi_nuclpred import make_figure
df = pd.read_csv(
    "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/invivo_landscape/figG1_summary_allfungi_nrl.tsv",
    sep="\t")
# map phylum_code to full names using clade_palette
from clade_palette import CLADE_NAMES, get_clade
df["phylum_code"] = df["species"].apply(get_clade)
df["phylum"]      = df["phylum_code"].map(CLADE_NAMES).fillna(df["phylum_code"])
make_figure(df)
print("Done G1")
EOF
