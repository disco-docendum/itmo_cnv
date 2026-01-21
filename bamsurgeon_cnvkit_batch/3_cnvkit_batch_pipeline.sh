#!/bin/bash

# Official github with cnvkit info:
# https://github.com/etal/cnvkit

# ==========================================
# DESCRIPTION
# ==========================================

# This pipeline scans for .bam files in its directory and performs batch cnv calling,
# based on the parameters in "CONFIGURATION" section.

# ==========================================
# INSTALLATION (MUST BE DONE ONLY ONCE)
# ==========================================

# set correct priority of conda channels, this must be done only once on a system:
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge
# conda config --set channel_priority strict

# create "cnvkit_env" environment for cnvkit and samtools only:
# conda create -n cnvkit_env cnvkit samtools -y

# ==========================================
# CONFIGURATION
# ==========================================

# Reference genome: "hg19", "hs37d5", or "hg38"
REF="hs37d5"

# Set CHR to "all" for WGS, or a specific number inside quotes (e.g., "20") for single chromosome mode
CHR="20"

# Bin size, list of bin sizes for typical .bam depth values:
# 60x:  500
# 30x:  1 000
# 20x:  1 500
# 15x:  2 000
# 10x:  3 000
# 5x:   6 000
# 3x:   10 000
# 1x:   30 000
# 0.5x: 60 000
BIN_SIZE=1000

# Blastlist exclusion of noisy regions, 1 = On, 0 = Off
# Theoretical info: https://github.com/Boyle-Lab/Blacklist
BLACKLIST=1         

# Threads per file
THREADS=4  

# ==========================================
# 0. DEPENDENCY CHECK
# ==========================================

echo "=========================================="
echo "CHECKING DEPENDENCIES"
echo "=========================================="

# Check for CNVkit
if ! command -v cnvkit.py &> /dev/null; then
    echo "❌ ERROR: 'cnvkit.py' not found."
    echo "   Please activate proper conda environment with cnvkit (e.g., 'conda activate cnvkit_env')."
    exit 1
else
    CNV_VER=$(cnvkit.py version 2>&1 | head -n 1) # Capture version safely
    echo "✅ Found CNVkit: $CNV_VER"
fi

# Check for Samtools
if ! command -v samtools &> /dev/null; then
    echo "❌ ERROR: 'samtools' not found."
    echo "   Please install samtools (e.g., 'conda install samtools')"
    exit 1
else
    SAM_VER=$(samtools --version | head -n 1)
    echo "✅ Found Samtools: $SAM_VER"
fi

# ==========================================
# DIRECTORY SETUP
# ==========================================

REF_DIR="./reference"
BLACKLIST_DIR="./blacklist"
VCF_DEFAULT="./vcf_default"
VCF_HMM="./vcf_hmm"
OUTPUT_ROOT="./output"

mkdir -p "$REF_DIR"
mkdir -p "$BLACKLIST_DIR"
mkdir -p "$VCF_DEFAULT"
mkdir -p "$VCF_HMM"
mkdir -p "$OUTPUT_ROOT"

echo "=========================================="
echo " 1. PREPARING REFERENCE (SHARED)"
echo "=========================================="

# Define URLs and Paths based on CHR and REF settings
if [ "$REF" == "hs37d5" ]; then
    # HS37D5 (b37 + decoys) - used by 1000G and Broad
    REF_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
    REF_RAW="${REF_DIR}/hs37d5_raw.fa.gz"
    REF_CLEAN="${REF_DIR}/hs37d5.fa"

    # Warning for Single Chromosome Mode with hs37d5
    if [ "$CHR" != "all" ]; then
        echo "⚠️  WARNING: Single chromosome mode ($CHR) requested for 'hs37d5'."
        echo "    'hs37d5' is typically distributed as a single full-genome file."
        echo "    Downloading the FULL genome reference regardless of CHR setting."
    fi

elif [ "$CHR" == "all" ]; then
    # Standard WGS for hg19/hg38
    REF_URL="https://hgdownload.cse.ucsc.edu/goldenpath/${REF}/bigZips/${REF}.fa.gz"
    REF_RAW="${REF_DIR}/${REF}_raw.fa.gz"
    REF_CLEAN="${REF_DIR}/${REF}.fa"
else
    # Single Chromosome for hg19/hg38
    REF_URL="https://hgdownload.cse.ucsc.edu/goldenpath/${REF}/chromosomes/chr${CHR}.fa.gz"
    REF_RAW="${REF_DIR}/chr${CHR}_raw.fa.gz"
    REF_CLEAN="${REF_DIR}/chr${CHR}.fa"
fi

# Download and format reference if missing
if [ ! -f "$REF_CLEAN" ]; then
    echo "Downloading Reference from $REF_URL..."
    wget -q -O "$REF_RAW" "$REF_URL"

    echo "Decompressing $REF..."
    # Simply decompress; removed sed 's/>chr/>/g' logic
    zcat "$REF_RAW" > "$REF_CLEAN"

    rm "$REF_RAW"
else
    echo "Reference already exists at $REF_CLEAN"
fi

# Index Reference
if [ ! -f "${REF_CLEAN}.fai" ]; then
    echo "Indexing Reference..."
    samtools faidx "$REF_CLEAN"
fi

# Prepare Blacklist if needed
if [ "$BLACKLIST" -eq 1 ]; then
    BLACKLIST_URL="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/${REF}-blacklist.v2.bed.gz"
    BLACKLIST_RAW="${BLACKLIST_DIR}/${REF}_blacklist.bed.gz"
    BLACKLIST_CLEAN="${BLACKLIST_DIR}/${REF}_blacklist_clean.bed"

    if [ ! -f "$BLACKLIST_CLEAN" ]; then
        echo "Downloading Blacklist..."
        wget -q -O "$BLACKLIST_RAW" "$BLACKLIST_URL"
        if [ "$REF" == "hg19" ]; then
            zcat "$BLACKLIST_RAW" | sed 's/^chr//g' > "$BLACKLIST_CLEAN"
        else
            zcat "$BLACKLIST_RAW" > "$BLACKLIST_CLEAN"
        fi
        rm "$BLACKLIST_RAW"
    fi
fi

echo "=========================================="
echo " 2. GENERATING SHARED TARGETS & REFS"
echo "=========================================="

# 1. Access File (Shared)
ACCESS_FILE="${REF_DIR}/access.bed"
if [ ! -f "$ACCESS_FILE" ]; then
    echo "Calculating Access..."
    if [ "$BLACKLIST" -eq 1 ]; then
        cnvkit.py access "$REF_CLEAN" -x "$BLACKLIST_CLEAN" -o "$ACCESS_FILE"
    else
        cnvkit.py access "$REF_CLEAN" -o "$ACCESS_FILE"
    fi
fi

# 2. Target File (Shared)
# We generate one target file for the whole batch
SHARED_TARGETS="${REF_DIR}/shared_targets.bed"
if [ ! -f "$SHARED_TARGETS" ]; then
    echo "Creating Shared Targets (Tiling)..."
    cnvkit.py target "$ACCESS_FILE" --split -a "$BIN_SIZE" -o "$SHARED_TARGETS"
fi

# 3. Flat Reference (Shared)
# We generate one flat reference for the whole batch
SHARED_FLAT_REF="${REF_DIR}/flat_reference.cnn"
if [ ! -f "$SHARED_FLAT_REF" ]; then
    echo "Creating Shared Flat Reference..."
    cnvkit.py reference -o "$SHARED_FLAT_REF" -f "$REF_CLEAN" -t "$SHARED_TARGETS"
fi

echo "=========================================="
echo " 3. BATCH PROCESSING BAM FILES"
echo "=========================================="

# Loop through all BAM files in current directory
for INPUT_BAM in *.bam; do

    # Skip if no bam files found
    [ -e "$INPUT_BAM" ] || continue

    BASENAME=$(basename "$INPUT_BAM" .bam)
    SAMPLE_OUT_DIR="${OUTPUT_ROOT}/${BASENAME}"
    
    echo ">>> Processing Sample: $BASENAME"
    mkdir -p "$SAMPLE_OUT_DIR"

    # A. Run Batch (using shared reference/targets)
    # We pass -r (reference) and -t (targets) so it skips recalculating them
    cnvkit.py batch "$INPUT_BAM" \
        -r "$SHARED_FLAT_REF" \
        -p "$THREADS" \
        --output-dir "$SAMPLE_OUT_DIR" \
        --scatter --diagram

    # B. HMM Segmentation (Extra step not in standard batch)
    cnvkit.py segment "$SAMPLE_OUT_DIR/${BASENAME}.cnr" \
        -m hmm-germline \
        -o "$SAMPLE_OUT_DIR/${BASENAME}_hmm.cns"

    # C. Export Default VCF
    cnvkit.py export vcf "$SAMPLE_OUT_DIR/${BASENAME}.cns" \
        -o "${VCF_DEFAULT}/${BASENAME}.vcf"

    # D. Export HMM VCF
    cnvkit.py export vcf "$SAMPLE_OUT_DIR/${BASENAME}_hmm.cns" \
        -o "${VCF_HMM}/${BASENAME}_hmm.vcf"

    echo "Finished $BASENAME"

done

# ==========================================
# 4. FINAL CLEANUP
# ==========================================

echo "Cleaning up shared temporary files..."
# Deleting flat_reference and targets at the end
if [ -f "$SHARED_FLAT_REF" ]; then
    rm "$SHARED_FLAT_REF"
fi

# We created shared_targets in the reference dir, if you want to delete them:
if [ -f "$SHARED_TARGETS" ]; then
    rm "$SHARED_TARGETS"
fi

echo "Pipeline Complete."