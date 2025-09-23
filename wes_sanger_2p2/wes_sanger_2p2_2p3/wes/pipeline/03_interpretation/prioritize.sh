#!/usr/bin/env bash
set -euo pipefail
VCF_IN=${1:-results/05_vep/PROBAND.vep.vcf}
OUT=results/06_prioritized
mkdir -p ${OUT}

bcftools view -i 'INFO/CSQ[*] ~ "HIGH|MODERATE"' ${VCF_IN} | \
bcftools view -i 'MIN(INFO/gnomADg_AF, INFO/gnomADe_AF) <= 0.001 || INFO/gnomADg_AF="." || INFO/gnomADe_AF="."' | \
bcftools view -i 'FMT/GT="1/1"' -Oz -o ${OUT}/PROBAND.recessive_pass.vcf.gz
tabix -p vcf ${OUT}/PROBAND.recessive_pass.vcf.gz
