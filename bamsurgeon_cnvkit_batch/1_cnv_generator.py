# ==========================================
#           DESCRIPTION
# ==========================================

# This python script generates txt files with simulated CNVs to generate new BAM files using Bamsurgeon.
# In current version, it needs to be manually configured in "user configuration" section to work.
# Just follow the comments to set up. The resulting txt files will be place into "sim" subdirectory in the location of the script.

# ==========================================
#           USER CONFIGURATION
# ==========================================

# Reference genome version: "hg38" or "hg19" (hs37d5 also should use "hg19", as dictionaries contain only autosomes)
REF_VERSION = "hg19"

# Chromosome to simulate on: "all" or specific name (e.g., "chr20", "chrX")
# Note: If using a specific chromosome, ensure it exists in the dictionary below
CHR_INPUT = "all"

# Blacklist filtering: 0 (off) or 1 (on)
# Theoretical info: https://github.com/Boyle-Lab/Blacklist
BLACKLIST_FLAG = 1

# Number of sample files to generate (sim1.txt, sim2.txt, etc.)
NUM_SAMPLES = 1

# haploid (0.5) or diploid (1) deletion mutations
DELPEPTH = 0.5

# number of extra copies from duplications mutations, integers only
# 1 is equivalent of extra single heterozygous copy
# 2 is equilvalent of homozygous duplication or 2 copies on the same chromosome, and so on
DUPDEPTH = 1

# CNV Matrix: List of [count, length] pairs.
# Example: [[20, 10_000], [10, 50_000]] means:
# - Generate 20 CNVs of length 10,000 bp
# - Generate 10 CNVs of length 50,000 bp
CNV_MATRIX = [
    [100, 5_000],
    [100, 10_000],
    [50, 50_000],
    [50, 100_000],
    [20, 300_000],
    [20, 500_000],
    [10, 1_000_000]
]

# ==========================================
#           SYSTEM SETUP
# ==========================================

import os
import sys
import random
import gzip
import requests

# Get the directory where THIS script is located
try:
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:
    BASE_DIR = os.getcwd()

print(f"Working directory: {BASE_DIR}")

# The hardcoded chromosome length were extracted from UCSC:
# hg19: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
# hg38: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

# Chromosome lengths in base pairs (bp) for UCSC hg19
hg19 = {
    'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276,
    'chr5': 180915260, 'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022,
    'chr9': 141213431, 'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895,
    'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753,
    'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr20': 63025520,
    'chr21': 48129895, 'chr22': 51304566
}

# Chromosome lengths in base pairs (bp) for UCSC hg38
hg38 = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
    'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
    'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
    'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
    'chr21': 46709983, 'chr22': 50818468
}

# Select Dictionary
if REF_VERSION == "hg19":
    genome_db = hg19
    blacklist_name = "hg19-blacklist.v2.bed.gz"
    blacklist_url = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz"
elif REF_VERSION == "hg38":
    genome_db = hg38
    blacklist_name = "hg38-blacklist.v2.bed.gz"
    blacklist_url = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz"
else:
    sys.exit("Error: REF_VERSION must be 'hg19' or 'hg38'.")

# Normalize CHR_INPUT
if CHR_INPUT != "all":
    if str(CHR_INPUT).isdigit(): 
        CHR_INPUT = f"chr{CHR_INPUT}"
    if CHR_INPUT not in genome_db:
        sys.exit(f"Error: Chromosome '{CHR_INPUT}' not found in {REF_VERSION}.")

# Sort Matrix
# This way it is easier to first squeeze large CNVs without overlapping with anything.
CNV_MATRIX.sort(key=lambda x: x[1], reverse=True)
print(f"Sorted CNV Matrix: {CNV_MATRIX}")

# ==========================================
#           BLACKLIST LOGIC
# ==========================================

blocked_regions = {} 

if BLACKLIST_FLAG == 1:
    blacklist_dir = os.path.join(BASE_DIR, "blacklist")
    blacklist_path = os.path.join(blacklist_dir, blacklist_name)

    # 1. Ensure directory exists
    if not os.path.exists(blacklist_dir):
        print(f"Creating blacklist directory: {blacklist_dir}")
        os.makedirs(blacklist_dir)

    # 2. Ensure file exists
    if not os.path.exists(blacklist_path):
        print(f"Downloading {blacklist_name}...")
        try:
            response = requests.get(blacklist_url)
            response.raise_for_status()
            with open(blacklist_path, 'wb') as f:
                f.write(response.content)
            print("Download complete.")
        except Exception as e:
            sys.exit(f"Error downloading blacklist: {e}")
    else:
        print(f"Using blacklist: {blacklist_path}")

    # 3. Parse file
    print("Parsing blacklist regions...")
    try:
        with gzip.open(blacklist_path, 'rt') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    b_chrom = parts[0]
                    b_start = int(parts[1])
                    b_end = int(parts[2])
                    if b_chrom not in blocked_regions:
                        blocked_regions[b_chrom] = []
                    blocked_regions[b_chrom].append((b_start, b_end))
    except Exception as e:
         sys.exit(f"Error reading blacklist file: {e}")

# ==========================================
#           GENERATION LOGIC
# ==========================================

output_dir = os.path.join(BASE_DIR, "sim")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def check_overlap(chrom, start, end, blocked_dict):
    """
    Generic overlap checker.
    Returns True if the interval (start, end) overlaps ANY interval in blocked_dict[chrom].
    blocked_dict format: {'chr1': [(s1, e1), (s2, e2)], ...}
    """
    if chrom not in blocked_dict: return False
    for b_start, b_end in blocked_dict[chrom]:
        if start < b_end and end > b_start:
            return True
    return False

print(f"Generating {NUM_SAMPLES} samples...")

padding_width = len(str(NUM_SAMPLES))

for i in range(1, NUM_SAMPLES + 1):
    file_number = str(i).zfill(padding_width)
    filename = os.path.join(output_dir, f"sim{file_number}.txt")
    # Track placed CNVs for THIS sample to prevent self-overlap
    # Format: {'chr1': [(start, end), ...], 'chr2': ...}
    placed_cnvs = {} 
    
    results = [] # Store rows to write later

    # Iterate through the sorted Matrix
    for group_idx, (count_target, cnv_len) in enumerate(CNV_MATRIX):
        
        generated_count = 0
        safety_break = 0
        
        while generated_count < count_target:
            safety_break += 1
            if safety_break > 200000:
                print(f"Warning: Sample {i}, Group {group_idx} (Len {cnv_len}): Space exhausted. Skipping.")
                break

            # 1. Pick Chromosome
            if CHR_INPUT == "all":
                target_chrom = random.choice(list(genome_db.keys()))
            else:
                target_chrom = CHR_INPUT
            
            chrom_len = genome_db[target_chrom]
            
            if chrom_len <= cnv_len:
                continue 

            # 2. Pick Coordinates
            start_pos = random.randint(0, chrom_len - cnv_len)
            end_pos = start_pos + cnv_len
            
            # 3. Check Blacklist (External noise)
            if BLACKLIST_FLAG == 1:
                if check_overlap(target_chrom, start_pos, end_pos, blocked_regions):
                    continue
            
            # 4. Check Overlap with Previously Placed CNVs (Internal collision)
            if check_overlap(target_chrom, start_pos, end_pos, placed_cnvs):
                continue

            # 5. Success - Register CNV
            if target_chrom not in placed_cnvs:
                placed_cnvs[target_chrom] = []
            placed_cnvs[target_chrom].append((start_pos, end_pos))

            cnv_type = random.choice(["DEL", "DUP"])

            if cnv_type == "DEL":
                depth = DELPEPTH
            else:
                depth = int(DUPDEPTH)
            
            # Prepare Output Name
            output_chrom = target_chrom
            if REF_VERSION == "hg19":
                output_chrom = output_chrom.replace("chr", "")
            
            results.append(f"{output_chrom}\t{start_pos}\t{end_pos}\t{cnv_type}\t{depth}\n")
            generated_count += 1

    # Write file
    with open(filename, 'w') as f:
        f.write("#Chrom\tStart\tEnd\tType\tDepth\n")
        for line in results:
            f.write(line)

print(f"Generation successful. Files located in: {output_dir}")