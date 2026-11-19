#!/usr/bin/env bash

set -euo pipefail

# Prepare ParaAT + KaKs_Calculator inputs per species for Ks estimation.


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../core/config_wgd.sh"

# If RUN_PARAAT=1, this script will also execute ParaAT+KaKs_Calculator
# for each species instead of only preparing the pairs.
RUN_PARAAT="${RUN_PARAAT:-0}"

PARALOG_TSV="${RESULTS_DIR}/ks/paralog_pairs_for_ks.filtered.tsv"
CDS_DIR="${RESULTS_DIR}/ks/cds"
OUT_ROOT="${RESULTS_DIR}/ks/paraAT"

if [ ! -f "${PARALOG_TSV}" ]; then
	echo "[WGD][Ks] ERROR: filtered paralog_pairs_for_ks TSV not found at ${PARALOG_TSV}" >&2
	exit 1
fi

if [ ! -d "${CDS_DIR}" ]; then
	echo "[WGD][Ks] ERROR: CDS directory not found at ${CDS_DIR}. Run 02b_extract_cds_for_ks.sh first." >&2
	exit 1
fi

mkdir -p "${OUT_ROOT}"

echo "[WGD][Ks] Preparing ParaAT inputs per species (pairs + paths)." >&2

# Unique species from paralog_pairs_for_ks.tsv (header 'species').
species_list=$(cut -f1 "${PARALOG_TSV}" | tail -n +2 | sort -u)

for sp in ${species_list}; do
	cds="${CDS_DIR}/${sp}.cds.fa"
	sp_out="${OUT_ROOT}/${sp}"
	pairs_file="${sp_out}/${sp}.pairs.txt"

	mkdir -p "${sp_out}"

	# Resolver el FASTA de proteínas probando varios sufijos habituales.
	prot=""
	for ext in .faa .fa .fasta .aa.fasta .aa.fa; do
		cand="${PROTEIN_DIR}/${sp}${ext}"
		if [ -f "${cand}" ]; then
			prot="${cand}"
			break
		fi
	done

	if [ -z "${prot}" ]; then
		echo "[WGD][Ks] WARNING: protein FASTA not found for ${sp} under ${PROTEIN_DIR} (tried .faa/.fa/.fasta/.aa.fasta/.aa.fa); skipping." >&2
		continue
	fi

	if [ ! -f "${cds}" ]; then
		echo "[WGD][Ks] WARNING: CDS FASTA not found for ${sp} at ${cds}; skipping." >&2
		continue
	fi

	# Build pair list for this species
	awk -v sp="${sp}" 'NR>1 && $1==sp {print $2"\t"$3}' "${PARALOG_TSV}" > "${pairs_file}"

	if [ ! -s "${pairs_file}" ]; then
		echo "[WGD][Ks] No paralog pairs for species ${sp}; produced empty ${pairs_file}, skipping." >&2
		rm -f "${pairs_file}"
		continue
	fi

	echo "[WGD][Ks] Prepared ParaAT pair list for ${sp} at ${pairs_file}" >&2
	echo "[WGD][Ks] Example ParaAT command for ${sp}:" >&2
	echo "  cd ${sp_out} && ParaAT.pl -a ${prot} -n ${cds} -p ${pairs_file} -k /path/to/KaKs_Calculator -o ${sp_out}" >&2

	if [ "${RUN_PARAAT}" = "1" ]; then
		if ! command -v ParaAT.pl >/dev/null 2>&1; then
			echo "[WGD][Ks] WARNING: ParaAT.pl not found in PATH; skipping execution for ${sp}." >&2
			continue
		fi

		KAKS_DIR_RESOLVED="${KAKS_DIR:-}"
		if [ -z "${KAKS_DIR_RESOLVED}" ]; then
			if command -v KaKs_Calculator >/dev/null 2>&1; then
				KAKS_DIR_RESOLVED="$(dirname "$(command -v KaKs_Calculator)")"
			else
				echo "[WGD][Ks] WARNING: KaKs_Calculator not found (and KAKS_DIR not set); skipping execution for ${sp}." >&2
				continue
			fi
		fi

		if [ -f "${sp_out}/.paraat_kaks.done" ]; then
			echo "[WGD][Ks] ParaAT+KaKs already run for ${sp} (marker .paraat_kaks.done); skipping." >&2
			continue
		fi

		echo "[WGD][Ks] Running ParaAT+KaKs for ${sp}..." >&2
		(
			cd "${sp_out}" && \
			ParaAT.pl -a "${prot}" -n "${cds}" -p "${pairs_file}" -k "${KAKS_DIR_RESOLVED}" -o "${sp_out}"
		) && touch "${sp_out}/.paraat_kaks.done"
	fi
done

echo "[WGD][Ks] ParaAT input preparation finished. Run ParaAT manually for selected species using the commands above." >&2
