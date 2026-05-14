#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

VCF="VCF/imputed/Broussonetia_RADseq.10_missing60_mac1_site_missing20_ID_chrnum.chr1.beagle4.1.vcf.gz"
OUT_DIR="${REPO_DIR}/reports/iqtree/chr1_beagle4.1_radseq_ufboot"
THREADS="${THREADS:-24}"
BOOTSTRAPS="${BOOTSTRAPS:-1000}"
MODEL="${MODEL:-GTR+ASC}"
CONVERT_ONLY="${CONVERT_ONLY:-0}"
IQTREE_BIN="${IQTREE_BIN:-${IQTREE:-}}"

if [[ -z "${IQTREE_BIN}" ]]; then
  if command -v iqtree2 >/dev/null 2>&1; then
    IQTREE_BIN="$(command -v iqtree2)"
  elif command -v iqtree >/dev/null 2>&1; then
    IQTREE_BIN="$(command -v iqtree)"
  elif [[ -x "${HOME}/bin/iqtree2" ]]; then
    IQTREE_BIN="${HOME}/bin/iqtree2"
  elif [[ -x "${HOME}/bin/iqtree" ]]; then
    IQTREE_BIN="${HOME}/bin/iqtree"
  elif [[ -x "${REPO_DIR}/tools/iqtree-download/iqtree-3.1.2-Linux/bin/iqtree3" ]]; then
    IQTREE_BIN="${REPO_DIR}/tools/iqtree-download/iqtree-3.1.2-Linux/bin/iqtree3"
  fi
fi

if [[ -z "${IQTREE_BIN}" || ! -x "${IQTREE_BIN}" ]]; then
  echo "ERROR: IQ-TREE was not found on PATH." >&2
  echo "Set IQTREE_BIN=/path/to/iqtree2, then rerun this script." >&2
  exit 127
fi

if [[ ! -f "${REPO_DIR}/${VCF}" ]]; then
  echo "ERROR: Input VCF not found: ${REPO_DIR}/${VCF}" >&2
  exit 1
fi

mkdir -p "${OUT_DIR}"

if ! command -v python3 >/dev/null 2>&1; then
  echo "ERROR: python3 is required to convert VCF genotypes to PHYLIP alignment." >&2
  exit 127
fi

ALIGNMENT="${OUT_DIR}/Broussonetia_RADseq.chr1.beagle4.1.snps.phy"
PREFIX="${OUT_DIR}/Broussonetia_RADseq.chr1.beagle4.1.${MODEL//+/_}.ufboot${BOOTSTRAPS}"

python3 - "${REPO_DIR}/${VCF}" "${ALIGNMENT}" <<'PY'
import gzip
import sys

vcf_path, out_path = sys.argv[1], sys.argv[2]
iupac = {
    frozenset(("A", "C")): "M",
    frozenset(("A", "G")): "R",
    frozenset(("A", "T")): "W",
    frozenset(("C", "G")): "S",
    frozenset(("C", "T")): "Y",
    frozenset(("G", "T")): "K",
}
states = {
    "A": set("A"),
    "C": set("C"),
    "G": set("G"),
    "T": set("T"),
    "M": set("AC"),
    "R": set("AG"),
    "W": set("AT"),
    "S": set("CG"),
    "Y": set("CT"),
    "K": set("GT"),
}

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

def gt_to_base(sample_field, gt_index, ref, alt):
    parts = sample_field.split(":")
    if gt_index >= len(parts):
        return "N"
    gt = parts[gt_index].replace("|", "/")
    if "." in gt:
        return "N"
    alleles = gt.split("/")
    bases = []
    for allele in alleles:
        if allele == "0":
            bases.append(ref)
        elif allele == "1":
            bases.append(alt)
        else:
            return "N"
    unique = frozenset(bases)
    if len(unique) == 1:
        return bases[0]
    return iupac.get(unique, "N")

samples = []
columns = []
kept = 0
skipped = 0
constant = 0

with open_text(vcf_path) as handle:
    for line in handle:
        line = line.rstrip("\n")
        if line.startswith("##"):
            continue
        fields = line.split("\t")
        if line.startswith("#CHROM"):
            samples = fields[9:]
            continue
        if not samples:
            raise SystemExit("VCF header with samples was not found.")

        ref = fields[3].upper()
        alt = fields[4].upper()
        if len(ref) != 1 or len(alt) != 1 or "," in alt or ref not in "ACGT" or alt not in "ACGT":
            skipped += 1
            continue

        fmt = fields[8].split(":")
        try:
            gt_index = fmt.index("GT")
        except ValueError:
            raise SystemExit("VCF FORMAT does not contain GT.")

        column = [gt_to_base(sample_field, gt_index, ref, alt) for sample_field in fields[9:]]
        possible_states = [states[base] for base in column if base in states]
        if not possible_states:
            skipped += 1
            continue

        shared_state = set.intersection(*possible_states)
        if shared_state:
            constant += 1
            continue

        columns.append(column)
        kept += 1

if not samples or kept == 0:
    raise SystemExit("No biallelic SNPs were available for the PHYLIP alignment.")

seqs = [[] for _ in samples]
for column in columns:
    for seq, base in zip(seqs, column):
        seq.append(base)

with open(out_path, "w") as out:
    out.write(f"{len(samples)} {kept}\n")
    for sample, seq in zip(samples, seqs):
        out.write(f"{sample} {''.join(seq)}\n")

print(f"Wrote PHYLIP alignment: {out_path}")
print(f"Samples: {len(samples)}")
print(f"Variable SNPs retained: {kept}")
print(f"Invariant SNPs removed: {constant}")
print(f"Sites skipped: {skipped}")
PY

echo "IQ-TREE: ${IQTREE_BIN}"
echo "Input alignment: ${ALIGNMENT}"
echo "Output prefix: ${PREFIX}"
echo "Model: ${MODEL}"
echo "Ultrafast bootstraps: ${BOOTSTRAPS}"
echo "Threads: ${THREADS}"

if [[ "${CONVERT_ONLY}" == "1" ]]; then
  echo "CONVERT_ONLY=1; skipping IQ-TREE run."
  exit 0
fi

"${IQTREE_BIN}" \
  -s "${ALIGNMENT}" \
  -st DNA \
  -m "${MODEL}" \
  -B "${BOOTSTRAPS}" \
  -T "${THREADS}" \
  --prefix "${PREFIX}"

echo "Finished. ML tree: ${PREFIX}.treefile"
