#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

BED_PREFIX="VCF/imputed/Broussonetia_RADseq.10_missing60_mac1_site_missing20_ID_chrnum.chr1.beagle4.1"
BED_FILE="${REPO_DIR}/${BED_PREFIX}.bed"
OUT_DIR="${REPO_DIR}/reports/admixture/chr1_beagle4.1_K1_20_cv10"
THREADS="${THREADS:-24}"
CV_FOLDS="${CV_FOLDS:-10}"
K_MIN="${K_MIN:-1}"
K_MAX="${K_MAX:-20}"

ADMIXTURE_BIN="${ADMIXTURE_BIN:-${ADMIXTURE:-}}"
if [[ -z "${ADMIXTURE_BIN}" ]]; then
  if command -v admixture >/dev/null 2>&1; then
    ADMIXTURE_BIN="$(command -v admixture)"
  elif [[ -x "${HOME}/bin/admixture" ]]; then
    ADMIXTURE_BIN="${HOME}/bin/admixture"
  elif [[ -x "${REPO_DIR}/tools/admixture-download/admixture_linux-1.4.0/admixture" ]]; then
    ADMIXTURE_BIN="${REPO_DIR}/tools/admixture-download/admixture_linux-1.4.0/admixture"
  elif [[ -x "/mnt/c/Users/kuowe/Physaria_SSR/tools/admixture-download/admixture_linux-1.4.0/admixture" ]]; then
    ADMIXTURE_BIN="/mnt/c/Users/kuowe/Physaria_SSR/tools/admixture-download/admixture_linux-1.4.0/admixture"
  fi
fi

if [[ -z "${ADMIXTURE_BIN}" || ! -x "${ADMIXTURE_BIN}" ]]; then
  echo "ERROR: ADMIXTURE was not found on PATH." >&2
  echo "Set ADMIXTURE_BIN=/path/to/admixture, then rerun this script." >&2
  exit 127
fi

for ext in bed bim fam; do
  if [[ ! -f "${REPO_DIR}/${BED_PREFIX}.${ext}" ]]; then
    echo "ERROR: Missing PLINK input file: ${REPO_DIR}/${BED_PREFIX}.${ext}" >&2
    exit 1
  fi
done

mkdir -p "${OUT_DIR}"

echo "ADMIXTURE: ${ADMIXTURE_BIN}"
echo "Input: ${BED_FILE}"
echo "Output directory: ${OUT_DIR}"
echo "K range: ${K_MIN}-${K_MAX}"
echo "CV folds: ${CV_FOLDS}"
echo "Threads: ${THREADS}"

cd "${OUT_DIR}"
for K in $(seq "${K_MIN}" "${K_MAX}"); do
  LOG_FILE="Broussonetia_RADseq.chr1.beagle4.1.K${K}.cv${CV_FOLDS}.log"
  echo "Running K=${K}; log: ${OUT_DIR}/${LOG_FILE}"
  "${ADMIXTURE_BIN}" \
    --cv="${CV_FOLDS}" \
    -j"${THREADS}" \
    "${BED_FILE}" \
    "${K}" \
    2>&1 | tee "${LOG_FILE}"
done

grep -h "CV error" ./*.log > "cv_errors.tsv" || true
echo "Finished. CV error summary: ${OUT_DIR}/cv_errors.tsv"
