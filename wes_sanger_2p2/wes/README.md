
# WES Variant Analysis (Section 2.2 — WES, Prioritization)

Pipeline templates for proband-only **WES on GRCh38** plus prioritization filters aligned to the manuscript.

## Quick start
```bash
conda env create -f env/environment.yml
bash pipeline/01_wes/run_wes.sh pipeline/config.yaml
bash pipeline/02_annotation/run_vep.sh results/04_vcf/PROBAND.norm.vcf.gz
bash pipeline/03_interpretation/prioritize.sh results/05_vep/PROBAND.vep.vcf
```

> Exact versions and dataset identifiers are recorded in Supplementary S1–S3 (see `../sanger_validation_2p3/supplementary/` or the paper’s repository root).
