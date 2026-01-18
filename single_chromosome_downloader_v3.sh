#!/bin/bash

# ==========================================
# DESCRIPTION
# ==========================================

# This script downloads chromosome 20 from NCBI server of Genome in a Bottle project.
# Streaming performs downsampling and filtration of discordant pairs.
# Parameters from SETTINGS must be double-checked before launching.

# ==========================================
# INSTALLATION (MUST BE DONE ONLY ONCE)
# ==========================================

# set correct priority of conda channels, this must be done only once on a system:
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge
# conda config --set channel_priority strict

# create "samtools_env" environment for picard and samtools:
# conda create -n picard_env picard samtools -y

# ==========================================
# SETTINGS
# ==========================================

# 0. This affects only samtools
THREADS=6

# 1. Genome Settings
CHR="20" # Target Chromosome name for the header

# 2. Downsampling Settings
DOWNSAMPLE_FRACTION=0.1 # Fraction of reads to keep (e.g., 0.1 = 10%)

# 3. Input/Output
BAM_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x_chr20.bam"
OUTPUT_BAM="HG002_Chr20_paired_reheadered.bam"

# address of hg38 chr20 300x:
# https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x_chr20.bam

# address of hg19 chr20 300x:
# https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x_chr20.bam

# ==========================================
# PIPELINE
# ==========================================
echo "--------------------------------------"
echo "Settings:"
echo "  Target URL: $BAM_URL"
echo "  Target Chr: $CHR"
echo "  Fraction:   $DOWNSAMPLE_FRACTION"
echo "--------------------------------------"
echo "Download and processing in progress..."

curl -L --retry 10 --retry-delay 5 --keepalive-time 60 "$BAM_URL" | \
picard DownsampleSam \
    --INPUT /dev/stdin \
    --OUTPUT /dev/stdout \
    --PROBABILITY 0.1 \
    --STRATEGY Chained \
    --VALIDATION_STRINGENCY SILENT \
    --QUIET true | \
samtools view -h -@ "$THREADS" - | \
awk -v c="$CHR" 'BEGIN {OFS="\t"} 
    /^@SQ/ {if ($0 ~ "SN:"c"\t" || $0 ~ "SN:chr"c"\t") print; next} 
    /^@/ {print; next} 
    ($7 == "=" || $7 == $3) {print}' | \
samtools view -b -@ "$THREADS" - > "$OUTPUT_BAM"
echo "Indexing new BAM..."
samtools index "$OUTPUT_BAM"

echo "Success! Created $OUTPUT_BAM"