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

echo "Cleanup complete."
