#!/bin/bash

# Official github with theoretical info:
# https://github.com/adamewing/bamsurgeon

# ==========================================
# INSTALLATION (MUST BE DONE ONLY ONCE)
# ==========================================

# set correct priority of conda channels, this must be done only once on a system:
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge
# conda config --set channel_priority strict

# create "bamsurgeon_env" environment for bamsurgeon and its dependancies only:
# conda create -n bamsurgeon_env bamsurgeon -y

# ==========================================
# HOW TO EXECUTE
# ==========================================

# Examples how to run this file:

# bash bamsurgeon.sh -i my_file_1.bam -r hg38 -c 22 -@ 8
# bash bamsurgeon.sh -i my_file_2.bam -r hg19 -c all

# ==========================================
# USAGE & ARGUMENT PARSING
# ==========================================

usage() {
    echo "Usage: $0 -i <input_bam> -r <ref_genome> -c <chromosome> [-@ <threads>]"
    echo ""
    echo "  -i  Input BAM file (Mandatory)"
    echo "  -r  Reference genome (e.g., 'hs37d5', 'hg19', 'hg38') (Mandatory)"
    echo "  -c  Chromosome (e.g., '20' or 'all') (Mandatory)"
    echo "  -@  Number of threads (Default: 4)"
    exit 1
}

# Default value
THREADS=4
CONDA_ENV_NAME="bamsurgeon_env" 

# Parse arguments
while getopts ":i:r:c:@:" opt; do
  case ${opt} in
    i) INPUT_BAM=${OPTARG} ;;
    r) REF=${OPTARG} ;;
    c) CHR=${OPTARG} ;;
    @) THREADS=${OPTARG} ;;
    \?) echo "   Invalid option: -$OPTARG" >&2; usage ;;
    :) echo "   Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

# Check mandatory variables
if [[ -z "$INPUT_BAM" || -z "$REF" || -z "$CHR" ]]; then
    echo "   Error: Missing mandatory arguments."
    usage
fi

echo "=========================================="
echo "CONFIGURATION"
echo "=========================================="
echo "Input BAM: $INPUT_BAM"
echo "Reference: $REF"
echo "Chromosome: $CHR"
echo "Threads:    $THREADS"

# ==========================================
# 0. DEPENDENCY CHECK
# ==========================================

echo "=========================================="
echo "CHECKING DEPENDENCIES"
echo "=========================================="

# Initialize conda
if command -v conda &> /dev/null; then
    CONDA_BASE=$(conda info --base)
    source "$CONDA_BASE/etc/profile.d/conda.sh"
else
    source "$HOME/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || \
    source "$HOME/anaconda3/etc/profile.d/conda.sh" 2>/dev/null
fi

# Activate env
echo "Activating conda environment: $CONDA_ENV_NAME"
conda activate $CONDA_ENV_NAME

# Check for BamSurgeon (addsv.py)
if ! command -v addsv.py &> /dev/null; then
    echo "   ERROR: 'addsv.py' (BamSurgeon) not found."
    echo "   Please set up proper conda environment (e.g., 'conda activate bamsurgeon_env')."
    exit 1
else
    echo "   Found BamSurgeon (addsv.py)"
fi

# Check for BWA
if ! command -v bwa &> /dev/null; then
    echo "   ERROR: 'bwa' not found."
    echo "   Please install BWA (e.g., 'conda install bwa')."
    exit 1
else
    BWA_VER=$(bwa 2>&1 | grep -i "Version" | head -n 1)
    echo "   Found BWA: $BWA_VER"
fi

# Check for Samtools
if ! command -v samtools &> /dev/null; then
    echo "   ERROR: 'samtools' not found."
    echo "   Please install Samtools (e.g., 'conda install samtools')."
    exit 1
else
    SAM_VER=$(samtools --version | head -n 1)
    echo "   Found Samtools: $SAM_VER"
fi

# ==========================================
# 1. PRE-FLIGHT CHECKS
# ==========================================

SIM_DIR="./sim"

if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM '$INPUT_BAM' not found."
    exit 1
fi

if [ ! -d "$SIM_DIR" ]; then
    echo "Error: Directory '$SIM_DIR' does not exist."
    exit 1
fi

shopt -s nullglob
TXT_FILES=("$SIM_DIR"/*.txt)
shopt -u nullglob

if [ ${#TXT_FILES[@]} -eq 0 ]; then
    echo "Error: No .txt files found in '$SIM_DIR'. Exiting."
    exit 1
fi

echo "Found ${#TXT_FILES[@]} variant list(s)."

# ==========================================
# 2. UNIVERSAL REFERENCE SETUP
# ==========================================

REF_DIR="./reference_${REF}"
mkdir -p "$REF_DIR"

echo "------------------------------------------------"
echo "Setting up Reference Genome ($REF)"
echo "------------------------------------------------"

# ------------------------------------
# LOGIC BLOCK: HS37D5 (Special handling)
# ------------------------------------
if [ "$REF" == "hs37d5" ]; then
    
    # 1000 Genomes FTP URL
    REF_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
    REF_RAW="${REF_DIR}/hs37d5_raw.fa.gz"

    # Define Output Clean Path
    if [ "$CHR" == "all" ]; then
        REF_CLEAN="${REF_DIR}/hs37d5.fa"
    else
        REF_CLEAN="${REF_DIR}/hs37d5_chr${CHR}.fa"
    fi

    # Check if final reference exists
    if [ ! -f "$REF_CLEAN" ]; then
        
        # Check if RAW file exists; if not, download it.
        if [ ! -f "$REF_RAW" ]; then
            echo "Downloading HS37D5 Raw Reference (Full)..."
            wget -q -O "$REF_RAW" "$REF_URL"
            if [ $? -ne 0 ]; then
                echo "   Error: Download failed."
                rm -f "$REF_RAW"
                exit 1
            fi
        else
            echo "Found existing HS37D5 Raw Reference ($REF_RAW). Skipping download."
        fi

        echo "Extracting reference sequence..."
        
        if [ "$CHR" == "all" ]; then
            # Extract everything
            gunzip -c "$REF_RAW" > "$REF_CLEAN"
        else
            # Extract SPECIFIC chromosome
            # hs37d5 uses headers like >1, >2 (no 'chr' prefix usually)
            echo "Extracting chromosome '$CHR' from full reference..."
            
            # Using gunzip -c to feed awk
            gunzip -c "$REF_RAW" | awk -v chr="$CHR" 'BEGIN {P=0} /^>/ {if ($1 == ">"chr) P=1; else P=0} {if (P) print $0}' > "$REF_CLEAN"
            
            # Check if extraction resulted in empty file (wrong chr name?)
            if [ ! -s "$REF_CLEAN" ]; then
                echo "   Error: Extraction resulted in an empty file."
                echo "   Did you use the correct chromosome number? (hs37d5 uses '1', not 'chr1')"
                rm "$REF_CLEAN"
                exit 1
            fi
        fi
    else
        echo "   Reference already exists at $REF_CLEAN"
    fi

# ------------------------------------
# LOGIC BLOCK: HG19 / HG38
# ------------------------------------
elif [[ "$REF" == "hg19" || "$REF" == "hg38" ]]; then

    if [ "$CHR" == "all" ]; then
        # WGS Mode
        if [[ "$REF" == "hg19" ]]; then
            REF_URL="https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz"
        else
            REF_URL="https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz"
        fi
        REF_RAW="${REF_DIR}/${REF}_raw.fa.gz"
        REF_CLEAN="${REF_DIR}/${REF}.fa"
    else
        # Single Chromosome Mode
        REF_URL="https://hgdownload.cse.ucsc.edu/goldenpath/${REF}/chromosomes/chr${CHR}.fa.gz"
        REF_RAW="${REF_DIR}/chr${CHR}_raw.fa.gz"
        REF_CLEAN="${REF_DIR}/chr${CHR}.fa"
    fi

    # Standard Download Logic for hg19/hg38
    if [ ! -f "$REF_CLEAN" ]; then
        echo "Reference '$REF_CLEAN' not found."
        echo "Downloading from: $REF_URL"
        
        wget -O "$REF_RAW" "$REF_URL"
        
        if [ $? -ne 0 ]; then
            echo "   Error: Download failed."
            rm -f "$REF_RAW"
            exit 1
        fi

        echo "Decompressing reference..."
        gunzip -c "$REF_RAW" > "$REF_CLEAN"
    else
        echo "   Reference already exists: $REF_CLEAN"
    fi

else
    echo "   Error: Unsupported reference '$REF'. Use 'hs37d5', 'hg19', or 'hg38'."
    exit 1
fi

# ------------------------------------
# INDEXING
# ------------------------------------
# Indexing (Only if missing)
if [ ! -f "${REF_CLEAN}.bwt" ]; then
    echo "Indexing for BWA (this may take a while)..."
    bwa index "$REF_CLEAN"
fi
if [ ! -f "${REF_CLEAN}.fai" ]; then
    echo "Indexing for Samtools..."
    samtools faidx "$REF_CLEAN"
fi

# ==========================================
# 3. BATCH EXECUTION
# ==========================================

for VAR_LIST in "${TXT_FILES[@]}"; do
    
    BASE_INPUT=$(basename "$INPUT_BAM" .bam)
    BASE_VAR=$(basename "$VAR_LIST" .txt)
    OUTPUT_BAM="${BASE_INPUT}_${BASE_VAR}.bam"
    
    # Create isolated temp directory for each run
    TEMP_DIR="bamsurgeon_tmp_${BASE_VAR}"
    rm -rf "$TEMP_DIR" 
    mkdir -p "$TEMP_DIR"
    TEMP_UNSORTED="${TEMP_DIR}/unsorted.bam"

    echo "------------------------------------------------"
    echo "Processing: $VAR_LIST"
    echo "------------------------------------------------"

    # Clean legacy temp file if present
    rm -rf addsv.tmp 

    addsv.py \
      -v "$VAR_LIST" \
      -f "$INPUT_BAM" \
      -r "$REF_CLEAN" \
      -o "$TEMP_UNSORTED" \
      --aligner mem \
      -p "$THREADS" \
      --alignerthreads "$THREADS" \
      --mindepth 1

    if [ ! -f "$TEMP_UNSORTED" ]; then
        echo "Error: BamSurgeon failed for $VAR_LIST"
        continue
    fi

    echo "Sorting and indexing..."
    samtools sort -@ "$THREADS" -o "$OUTPUT_BAM" "$TEMP_UNSORTED"
    samtools index "$OUTPUT_BAM"

    # Cleanup
    rm -rf "$TEMP_DIR"
    rm -rf addsv.tmp
    rm "$TEMP_UNSORTED" 2>/dev/null
    
    echo "Finished: $OUTPUT_BAM"

done

echo "CNV insertion is complete."
