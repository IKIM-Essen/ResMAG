#!/usr/bin/env bash
set -euo pipefail

RESULTS_DIR="${1:?Usage: cleanup_results.sh <results_dir>}"

echo "Cleaning up in ${RESULTS_DIR}"

rm -rf "${RESULTS_DIR}"/binning/
rm -rf "${RESULTS_DIR}"/binning_prep/
rm -rf "${RESULTS_DIR}"/genomad/
rm -rf "${RESULTS_DIR}"/human_filtering/
rm -rf "${RESULTS_DIR}"/megahit/
rm -rf "${RESULTS_DIR}"/qc/
rm -rf "${RESULTS_DIR}"/trimmed/
rm -rf "${RESULTS_DIR}"/cleanup/

rm -f "${RESULTS_DIR}"/output/classification/*/*/kaiju.out
rm -f "${RESULTS_DIR}"/output/classification/*/*/kaiju.out.krona

rm -rf "${RESULTS_DIR}"/output/classification/bins/*/gtdbtk/

rm -f "${RESULTS_DIR}"/output/filtered_reads/*_R*.fastq

rm -f "${RESULTS_DIR}"/output/proteins/*/*.gff
rm -f "${RESULTS_DIR}"/output/proteins/*/*.fna
rm -f "${RESULTS_DIR}"/output/proteins/*/*.faa

rm -f "${RESULTS_DIR}"/output/report/prerequisites/assembly/*.bam
rm -f "${RESULTS_DIR}"/output/report/prerequisites/assembly/*.bam.bai

rm -rf "${RESULTS_DIR}"/output/report/*/bin/
rm -rf "${RESULTS_DIR}"/output/report/*/mags/
rm -rf "${RESULTS_DIR}"/output/report/*/taxonomy/

rm -f "${RESULTS_DIR}"/output/resistance/uniCARD/assembly/*_out.tsv
rm -f "${RESULTS_DIR}"/output/resistance/uniCARD/assembly/*_out.tsv.gz

rm -f "${RESULTS_DIR}"/../CARD_annotation.done
rm -f "${RESULTS_DIR}"/../CARD_load_DB_for_reads.done
rm -f "${RESULTS_DIR}"/../GTDB_prep.done

echo "Cleanup complete."
