#!/usr/bin/env bash
set -euo pipefail

VCF="VCF/filtered/Broussonetia_RADseq.09_missing50_mac1_site_missing20.vcf.gz"
OUT_DIR="reports/missing20_100kb_windows_ld/ld_decay_plink2"
PREFIX="${OUT_DIR}/missing20"
LD_PREFIX="${PREFIX}.ld_0_1000kb"
PLINK1_PREFIX="${PREFIX}.ld_0_1000kb_plink1"

PLINK2_BIN="${PLINK2_BIN:-}"
if [[ -z "${PLINK2_BIN}" ]]; then
  if command -v plink2 >/dev/null 2>&1; then
    PLINK2_BIN="$(command -v plink2)"
  elif [[ -x "/Applications/bin/plink2" ]]; then
    PLINK2_BIN="/Applications/bin/plink2"
  fi
fi

if [[ -z "${PLINK2_BIN}" || ! -x "${PLINK2_BIN}" ]]; then
  echo "ERROR: plink2 was not found on PATH." >&2
  echo "You can also set PLINK2_BIN=/path/to/plink2." >&2
  echo "Install plink2, then rerun: bash scripts/calculate_ld_decay_missing20_plink2.sh" >&2
  exit 127
fi

mkdir -p "${OUT_DIR}"

"${PLINK2_BIN}" \
  --vcf "${VCF}" \
  --double-id \
  --allow-extra-chr \
  --set-all-var-ids '@:#:$r:$a' \
  --make-pgen \
  --out "${PREFIX}"

if "${PLINK2_BIN}" --help 2>/dev/null | grep -q -- "--r2"; then
  "${PLINK2_BIN}" \
    --pfile "${PREFIX}" \
    --allow-extra-chr \
    --r2-unphased cols=chrom,pos,dist,r2 \
    --ld-window-kb 1000 \
    --ld-window-r2 0 \
    --out "${LD_PREFIX}"
else
  PLINK1_BIN="${PLINK1_BIN:-}"
  if [[ -z "${PLINK1_BIN}" ]]; then
    if command -v plink >/dev/null 2>&1; then
      PLINK1_BIN="$(command -v plink)"
    elif [[ -x "/Applications/bin/plink" ]]; then
      PLINK1_BIN="/Applications/bin/plink"
    fi
  fi

  if [[ -z "${PLINK1_BIN}" || ! -x "${PLINK1_BIN}" ]]; then
    echo "ERROR: this plink2 build cannot emit pairwise r2, and plink 1.9 was not found for fallback LD calculation." >&2
    exit 127
  fi

  "${PLINK1_BIN}" \
    --vcf "${VCF}" \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --r2 gz \
    --ld-window-kb 1000 \
    --ld-window 999999 \
    --ld-window-r2 0 \
    --out "${PLINK1_PREFIX}"
fi

Rscript scripts/plot_ld_decay_missing20_plink2.R
