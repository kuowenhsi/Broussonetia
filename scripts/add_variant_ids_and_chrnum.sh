#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 input.vcf.gz output.vcf.gz" >&2
  exit 1
fi

input_vcf="$1"
output_vcf="$2"
map_file="VCF/filtered/scaffold_to_chrnum.tsv"

if [[ ! -f "$input_vcf" ]]; then
  echo "Input VCF not found: $input_vcf" >&2
  exit 1
fi

mkdir -p "$(dirname "$output_vcf")"
mkdir -p "$(dirname "$map_file")"

bcftools view -h "$input_vcf" |
  awk -F'[=,]' '/^##contig/ && $3 ~ /^scaffold_/ {split($3,a,"_"); print $3"\t"a[2]}' \
  > "$map_file"

bcftools annotate --set-id '%CHROM\_%POS' -Ou "$input_vcf" |
  bcftools annotate --rename-chrs "$map_file" -Oz -o "$output_vcf"

bcftools index -f "$output_vcf"

echo "Wrote: $output_vcf"
echo "Wrote: ${output_vcf}.csi"
echo "Chromosome rename map: $map_file"
