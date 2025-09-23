#!/usr/bin/env bash
set -euo pipefail
CFG=${1:-pipeline/config.yaml}
echo "Using config: $CFG"

sample=$(grep sample_id $CFG | awk '{print $2}')
R1=$(grep fastq_R1 $CFG | awk '{print $2}')
R2=$(grep fastq_R2 $CFG | awk '{print $2}')
REF=$(grep reference_fasta $CFG | awk '{print $2}')
OUT=$(grep outdir $CFG | awk '{print $2}')
threads=$(grep threads $CFG | awk '{print $2}')

mkdir -p ${OUT}/{00_fastp,01_align,02_bqsr,03_gvcf,04_vcf}

fastp -i "$R1" -I "$R2" -o ${OUT}/00_fastp/${sample}_R1.trim.fq.gz -O ${OUT}/00_fastp/${sample}_R2.trim.fq.gz \
      -h ${OUT}/00_fastp/${sample}.fastp.html -j ${OUT}/00_fastp/${sample}.fastp.json -w $threads

bwa-mem2 mem -t $threads -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" "$REF" \
      ${OUT}/00_fastp/${sample}_R1.trim.fq.gz ${OUT}/00_fastp/${sample}_R2.trim.fq.gz | \
      samtools sort -@ $threads -o ${OUT}/01_align/${sample}.sorted.bam
samtools index ${OUT}/01_align/${sample}.sorted.bam

gatk MarkDuplicatesSpark -I ${OUT}/01_align/${sample}.sorted.bam -O ${OUT}/01_align/${sample}.md.bam \
     -M ${OUT}/01_align/${sample}.md.metrics
samtools index ${OUT}/01_align/${sample}.md.bam

gatk BaseRecalibrator -R "$REF" -I ${OUT}/01_align/${sample}.md.bam \
    $(grep -A10 known_sites $CFG | awk '{print $2}' | sed '/^$/d' | sed 's/^/-known-sites /') \
    -O ${OUT}/02_bqsr/${sample}.recal.table
gatk ApplyBQSR -R "$REF" -I ${OUT}/01_align/${sample}.md.bam \
    --bqsr-recal-file ${OUT}/02_bqsr/${sample}.recal.table \
    -O ${OUT}/02_bqsr/${sample}.recal.bam
samtools index ${OUT}/02_bqsr/${sample}.recal.bam

gatk HaplotypeCaller -R "$REF" -I ${OUT}/02_bqsr/${sample}.recal.bam -ERC GVCF \
    -O ${OUT}/03_gvcf/${sample}.g.vcf.gz
gatk GenotypeGVCFs -R "$REF" -V ${OUT}/03_gvcf/${sample}.g.vcf.gz -O ${OUT}/04_vcf/${sample}.raw.vcf.gz

bcftools norm -m -any -f "$REF" ${OUT}/04_vcf/${sample}.raw.vcf.gz -Oz -o ${OUT}/04_vcf/${sample}.norm.vcf.gz
tabix -p vcf ${OUT}/04_vcf/${sample}.norm.vcf.gz

echo "Next: annotation (pipeline/02_annotation/run_vep.sh)"
