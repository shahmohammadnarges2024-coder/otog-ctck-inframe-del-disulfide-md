#!/usr/bin/env bash
set -euo pipefail

CFG=${1:-pipeline/config.yaml}
echo "[INFO] Using config: $CFG"

need() { command -v "$1" >/dev/null 2>&1 || { echo "[ERR] $1 not found"; exit 127; }; }
for t in fastp bwa-mem2 samtools gatk bcftools tabix awk grep; do need "$t"; done

get_yaml() { awk -F': *' -v k="$1" '$1==k {print $2}' "$CFG"; }

sample=$(get_yaml sample_id)
R1=$(get_yaml fastq_R1)
R2=$(get_yaml fastq_R2)
REF=$(get_yaml reference_fasta)
OUT=$(get_yaml outdir)
threads=$(get_yaml threads)

mkdir -p "${OUT}"/{00_fastp,01_align,02_bqsr,03_gvcf,04_vcf}

echo "[STEP] fastp"
fastp -i "$R1" -I "$R2" \
  -o "${OUT}/00_fastp/${sample}_R1.trim.fq.gz" \
  -O "${OUT}/00_fastp/${sample}_R2.trim.fq.gz" \
  -h "${OUT}/00_fastp/${sample}.fastp.html" \
  -j "${OUT}/00_fastp/${sample}.fastp.json" -w "$threads"

echo "[STEP] bwa-mem2 + sort + index"
bwa-mem2 mem -t "$threads" -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" "$REF" \
  "${OUT}/00_fastp/${sample}_R1.trim.fq.gz" "${OUT}/00_fastp/${sample}_R2.trim.fq.gz" \
  | samtools sort -@ "$threads" -o "${OUT}/01_align/${sample}.sorted.bam"
samtools index "${OUT}/01_align/${sample}.sorted.bam"

echo "[STEP] MarkDuplicatesSpark"
gatk MarkDuplicatesSpark \
  -I "${OUT}/01_align/${sample}.sorted.bam" \
  -O "${OUT}/01_align/${sample}.md.bam" \
  -M "${OUT}/01_align/${sample}.md.metrics"
samtools index "${OUT}/01_align/${sample}.md.bam"

echo "[STEP] BQSR"
# known_sites از YAML: خطوط پس از کلید known_sites را بگیر
mapfile -t KS < <(awk '/^known_sites:/{f=1;next} f&&NF{print $2} /^$/{f=0}' "$CFG")
ks_args=()
for k in "${KS[@]}"; do ks_args+=("--known-sites" "$k"); done

gatk BaseRecalibrator -R "$REF" -I "${OUT}/01_align/${sample}.md.bam" \
  "${ks_args[@]}" -O "${OUT}/02_bqsr/${sample}.recal.table"
gatk ApplyBQSR -R "$REF" -I "${OUT}/01_align/${sample}.md.bam" \
  --bqsr-recal-file "${OUT}/02_bqsr/${sample}.recal.table" \
  -O "${OUT}/02_bqsr/${sample}.recal.bam"
samtools index "${OUT}/02_bqsr/${sample}.recal.bam"

echo "[STEP] HaplotypeCaller → GenotypeGVCFs"
gatk HaplotypeCaller -R "$REF" -I "${OUT}/02_bqsr/${sample}.recal.bam" -ERC GVCF \
  -O "${OUT}/03_gvcf/${sample}.g.vcf.gz"
gatk GenotypeGVCFs -R "$REF" -V "${OUT}/03_gvcf/${sample}.g.vcf.gz" \
  -O "${OUT}/04_vcf/${sample}.raw.vcf.gz"

echo "[STEP] Normalize/split/left-align + index"
bcftools norm -m -any -f "$REF" "${OUT}/04_vcf/${sample}.raw.vcf.gz" -Oz \
  -o "${OUT}/04_vcf/${sample}.norm.vcf.gz"
tabix -p vcf "${OUT}/04_vcf/${sample}.norm.vcf.gz"
bcftools stats -s - "${OUT}/04_vcf/${sample}.norm.vcf.gz" > "${OUT}/04_vcf/${sample}.stats.txt" || true

echo "[NEXT] Annotation → pipeline/02_annotation/run_vep.sh"



