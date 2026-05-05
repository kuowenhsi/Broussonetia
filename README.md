# Broussonetia RADseq QC

This repository contains R scripts and generated reports for exploring SNP filtering
conditions and chromosome-scale variant distributions in a Broussonetia RADseq VCF.

## Contents

- `scripts/render_snpfiltr_radseq_qc.R`: renders the SNPfiltR RADseq QC report.
- `reports/snpfiltr_radseq_qc.Rmd`: R Markdown source for the SNPfiltR report.
- `reports/snpfiltr_radseq_qc.html`: rendered exploratory SNP filtering report.
- `scripts/plot_variant_distribution_scaffolds.R`: bins variant sites across the first
  18 scaffolds and plots chromosome-style density tracks.
- `reports/variant_distribution_scaffolds/`: plots and CSV summaries from the
  scaffold variant-density analysis.
- `PRJNA437223_Broussonetia_meta.csv`: sample metadata used for population labels.

## Large Local Inputs

Large genomics inputs are intentionally excluded from normal Git history because GitHub
rejects files over 100 MB and Git LFS is not installed in this workspace. The scripts
expect these files to be present locally:

- `VCF/Broussonetia_RADseq.07_tags_refilled.vcf.gz`
- `VCF/Broussonetia_RADseq.07_tags_refilled.vcf.gz.csi`
- `genome_AS/B.papyrifera_yahs.out_scaffolds.fa`
- `genome_AS/B.papyrifera_yahs.out_scaffolds.fa.fai`

The `.fai` and `.csi` index files are small enough to be tracked; the VCF and FASTA
are ignored.

## Re-run Reports

Install required packages into a local library if needed:

```r
dir.create("R_libs", showWarnings = FALSE)
.libPaths(c("R_libs", .libPaths()))
install.packages(c("SNPfiltR", "vcfR", "tidyverse", "rmarkdown", "knitr", "DT", "patchwork", "scales"))
```

Render the SNPfiltR QC report:

```sh
Rscript scripts/render_snpfiltr_radseq_qc.R
```

Plot variant density across the first 18 scaffolds:

```sh
Rscript scripts/plot_variant_distribution_scaffolds.R
```
