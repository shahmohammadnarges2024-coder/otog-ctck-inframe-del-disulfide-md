#!/usr/bin/env bash
set -euo pipefail
VCF_IN=${1:-results/04_vcf/PROBAND.norm.vcf.gz}
OUT=results/05_vep
CACHE=/path/to/.vep
FASTA=/path/to/GRCh38.fa
mkdir -p ${OUT}

vep --cache --offline --dir_cache ${CACHE} --assembly GRCh38 \
    --vcf --everything --fasta ${FASTA} \
    -i ${VCF_IN} -o ${OUT}/PROBAND.vep.vcf
