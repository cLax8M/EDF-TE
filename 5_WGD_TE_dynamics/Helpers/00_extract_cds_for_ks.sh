#!/usr/bin/env bash

set -euo pipefail

# Extract CDS FASTA per species for Ks calculation.


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../core/config_wgd.sh"

PARALOG_TSV="${RESULTS_DIR}/ks/paralog_pairs_for_ks.tsv"
OUT_DIR="${RESULTS_DIR}/ks/cds"

if [ ! -f "${PARALOG_TSV}" ]; then
	echo "[WGD][Ks] ERROR: paralog_pairs_for_ks.tsv not found at ${PARALOG_TSV}" >&2
	exit 1
fi

mkdir -p "${OUT_DIR}"

# Resolver binario de gffread:
#   1) variable GFFREAD_BIN (si está definida)
#   2) $HOME/micromamba/envs/wgd_tools/bin/gffread (instalación típica aquí)
#   3) buscar gffread en PATH
GFFREAD_BIN="${GFFREAD_BIN:-"$HOME/micromamba/envs/wgd_tools/bin/gffread"}"
if [ ! -x "${GFFREAD_BIN}" ]; then
	if command -v gffread >/dev/null 2>&1; then
		GFFREAD_BIN="$(command -v gffread)"
	else
		echo "[WGD][Ks] ERROR: gffread binary not found. Set GFFREAD_BIN or ensure gffread is in PATH." >&2
		exit 1
	fi
fi

echo "[WGD][Ks] Extracting CDS per species using GFF+genome" >&2

# Get unique species list from the first column (header 'species').
species_list=$(cut -f1 "${PARALOG_TSV}" | tail -n +2 | sort -u)

for sp in ${species_list}; do
	# Quitar sufijo .aa si existe para mapear a los nombres de genoma/GFF
	base="${sp%.aa}"
	base_no_ver="${base%[0-9]}"  # e.g. AS_AspfumUD1 -> AS_AspfumUD

	# Resolver genoma probando distintas bases y extensiones
	genome=""
	for b in "${base}" "${base_no_ver}"; do
		[ -n "${genome}" ] && break
		# Probamos primero las extensiones "normales" y también las versiones
		# preprocesadas que terminan en .prep (p. ej. HO_Drome.fasta.prep)
		for ext in ".fasta" ".fa" ".fna" ".fasta.prep" ".fa.prep" ".fna.prep"; do
			cand="${GENOME_DIR}/${b}${ext}"
			if [ -f "${cand}" ]; then
				genome="${cand}"
				break
			fi
		done
	done

	# Resolver GFF probando .gff3 y .gff
	gff=""
	for b in "${base}" "${base_no_ver}"; do
		[ -n "${gff}" ] && break
		for ext in ".gff3" ".gff"; do
			cand="${GFF_DIR}/${b}${ext}"
			if [ -f "${cand}" ]; then
				gff="${cand}"
				break
			fi
		done
	done

	out_cds="${OUT_DIR}/${sp}.cds.fa"

	if [ -z "${genome}" ] || [ -z "${gff}" ]; then
		echo "[WGD][Ks] WARNING: missing genome or GFF for species ${sp}; tried bases '${base}', '${base_no_ver}' in ${GENOME_DIR} and ${GFF_DIR}. Skipping." >&2
		continue
	fi

	if [ -s "${out_cds}" ]; then
		echo "[WGD][Ks] CDS FASTA already exists for ${sp} at ${out_cds}, skipping." >&2
		continue
	fi

	echo "[WGD][Ks] Extracting CDS for ${sp} using genome=${genome}, gff=${gff} (${GFFREAD_BIN})" >&2
	if ! "${GFFREAD_BIN}" "${gff}" -g "${genome}" -x "${out_cds}"; then
		echo "[WGD][Ks] WARNING: gffread failed for ${sp}; skipping this species." >&2
		rm -f "${out_cds}"
		continue
	fi

	# A veces gffread puede terminar sin error pero no extraer
	# ninguna CDS (por ejemplo, en GFF con atributos no estándar
	# o sin transcritos bien definidos). En ese caso el fichero
	# resultante queda con tamaño 0; intentamos una ruta de
	# rescate específica para GFF de JGI basados en "proteinId".
	if [ ! -s "${out_cds}" ]; then
		echo "[WGD][Ks] WARNING: gffread produced empty CDS FASTA for ${sp}; trying JGI-style fallback extractor." >&2
		prot_fa=""
		for ext in .faa .fa .fasta .aa.fasta .aa.fa; do
			cand="${PROTEIN_DIR}/${sp}${ext}"
			if [ -f "${cand}" ]; then
				prot_fa="${cand}"
				break
			fi
		done
		if [ -z "${prot_fa}" ]; then
			echo "[WGD][Ks] WARNING: no protein FASTA found for ${sp} when trying JGI fallback; deleting empty CDS and skipping." >&2
			rm -f "${out_cds}"
			continue
		fi
		if ! python "${SCRIPT_DIR}/extract_cds_from_jgi_gff.py" \
			--genome "${genome}" \
			--gff "${gff}" \
			--proteins "${prot_fa}" \
			--out-cds "${out_cds}"; then
			echo "[WGD][Ks] WARNING: JGI-style CDS extractor failed for ${sp}; deleting empty CDS and skipping." >&2
			rm -f "${out_cds}"
			continue
		fi
		if [ ! -s "${out_cds}" ]; then
			echo "[WGD][Ks] WARNING: JGI-style CDS extractor produced empty CDS for ${sp}; deleting file and skipping." >&2
			rm -f "${out_cds}"
			continue
		fi
	fi
done

echo "[WGD][Ks] CDS extraction step finished" >&2
