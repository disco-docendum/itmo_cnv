#!/usr/bin/env python3
# ==========================================
#           DESCRIPTION
# ==========================================

# This python script generates txt files with simulated CNVs specifications to generate new BAM files using Bamsurgeon.
# UPDATED: Cross-platform compatible (Windows/Linux), no external dependencies.
# resulting txt files will be placed into "sim" subdirectory of the script location.

# ==========================================
#           USER CONFIGURATION
# ==========================================

# Reference genome version: "hg38" or "hg19" (hs37d5 also should use "hg19", as dictionaries contain only nuclear autosomes)
REF_VERSION = "hg19"

# Chromosome to simulate on: "all" or specific name (e.g., "chr20", "chrX")
# Note: If using a specific chromosome, ensure it exists in the dictionary below.
CHR_INPUT = "chr1"

# CNV Matrix: List of [count, length] pairs.
CNV_MATRIX = [
    [100, 5_000],
    [100, 10_000],
    [50, 50_000],
    [50, 100_000],
    [30, 300_000],
    [30, 500_000],
    [20, 1_000_000]
]

# Number of sample files/patients to generate (sim1.txt, sim2.txt, etc.)
NUM_SAMPLES = 1

# -------------- TECHNICAL SETTINGS OF FILTRATION --------------

# Blacklist filtering: 0 (off) or 1 (on)
# Theoretical info: https://github.com/Boyle-Lab/Blacklist
BLACKLIST_FLAG = 1

# Simple Repeats filtering: 0 (off) or 1 (on)
# Downloads UCSC simple repeats track to avoid generating CNVs in highly repetitive regions.
# This significantly reduces bamsurgeon assembly errors.
AVOID_REPEATS = 1

# Safety buffer (bp) around simple repeats for endpoints
# Ensures reads overlapping the breakpoint don't map to repeats.
REPEAT_BUFFER = 800

# haploid (0.5) or diploid (1) deletion mutations
DELPEPTH = 0.5

# number of extra copies (1 is heterozygous) from duplications mutations, integers only
DUPDEPTH = 1

# ==========================================
#           SYSTEM SETUP
# ==========================================

import os
import sys
import random
import gzip
import bisect
import shutil
import urllib.request
from pathlib import Path

# Cross-platform clear screen
def clear_screen():
    os.system('cls' if os.name == 'nt' else 'clear')

clear_screen()
print("Initializing CNV Generator...")

# Determine Base Directory safely
try:
    BASE_DIR = Path(__file__).resolve().parent
except NameError:
    BASE_DIR = Path.cwd()

print(f"Working directory: {BASE_DIR}")

# Chromosome lengths
hg19 = {
    'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276,
    'chr5': 180915260, 'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022,
    'chr9': 141213431, 'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895,
    'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753,
    'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr20': 63025520,
    'chr21': 48129895, 'chr22': 51304566
}

hg38 = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
    'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
    'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
    'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
    'chr21': 46709983, 'chr22': 50818468
}

if REF_VERSION == "hg19":
    genome_db = hg19
    blacklist_name = "hg19-blacklist.v2.bed.gz"
    blacklist_url = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz"
    simplerepeat_url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz"
elif REF_VERSION == "hg38":
    genome_db = hg38
    blacklist_name = "hg38-blacklist.v2.bed.gz"
    blacklist_url = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz"
    simplerepeat_url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz"
else:
    sys.exit("Error: REF_VERSION must be 'hg19' or 'hg38'.")

if CHR_INPUT != "all":
    if str(CHR_INPUT).isdigit(): 
        CHR_INPUT = f"chr{CHR_INPUT}"
    if CHR_INPUT not in genome_db:
        sys.exit(f"Error: Chromosome '{CHR_INPUT}' not found in {REF_VERSION}.")

CNV_MATRIX.sort(key=lambda x: x[1], reverse=True)
print(f"Sorted CNV Matrix: {CNV_MATRIX}")

# ==========================================
#           BLOCKLIST MANAGEMENT
# ==========================================

blacklist_dict = {} 
repeats_dict = {}

def download_file(url, filepath):
    """
    Downloads a file using standard library urllib (no pip requests required).
    """
    filepath = Path(filepath)
    if not filepath.exists():
        print(f"Downloading {filepath.name}...")
        try:
            with urllib.request.urlopen(url) as response, open(filepath, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)
            print("Done.")
        except Exception as e:
            sys.exit(f"Error downloading {url}: {e}")
    else:
        print(f"Using existing file: {filepath}")

def load_intervals(filepath, target_dict, file_format):
    """
    Parses files into the target dictionary based on format.
    file_format: 
      'bed'  -> standard BED (chrom, start, end, ...) -> Used for Blacklist
      'ucsc' -> UCSC dump (bin, chrom, start, end, ...) -> Used for Simple Repeats
    """
    print(f"Parsing {Path(filepath).name}...")
    try:
        # 'rt' with encoding='utf-8' ensures cross-platform text reading
        with gzip.open(filepath, 'rt', encoding='utf-8') as f:
            for line in f:
                # Skip comments
                if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                    continue
                
                p = line.strip().split('\t')
                
                try:
                    if file_format == 'bed':
                        # Standard BED: chrom(0), start(1), end(2)
                        if len(p) < 3: continue
                        c, s, e = p[0], int(p[1]), int(p[2])
                        
                    elif file_format == 'ucsc':
                        # UCSC Database dump: bin(0), chrom(1), start(2), end(3)
                        if len(p) < 4: continue
                        c, s, e = p[1], int(p[2]), int(p[3])
                    
                    else:
                        print("Unknown file format specified.")
                        return

                    if c in genome_db:
                        if c not in target_dict: target_dict[c] = []
                        target_dict[c].append((s, e))
                        
                except ValueError:
                    # Skips lines where start/end are not integers (like headers)
                    continue
                    
    except Exception as e:
        sys.exit(f"Error reading file {filepath}: {e}")

# Load Files
if BLACKLIST_FLAG == 1:
    bl_dir = BASE_DIR / "blacklist"
    bl_dir.mkdir(parents=True, exist_ok=True)
    bl_path = bl_dir / blacklist_name
    
    download_file(blacklist_url, bl_path)
    load_intervals(bl_path, blacklist_dict, file_format='bed')

if AVOID_REPEATS == 1:
    sr_dir = BASE_DIR / ("simple_repeat_" + REF_VERSION)
    sr_dir.mkdir(parents=True, exist_ok=True)
    sr_path = sr_dir / "simpleRepeat.txt.gz"
    
    download_file(simplerepeat_url, sr_path)
    load_intervals(sr_path, repeats_dict, file_format='ucsc')

# Optimization: Merge Intervals
def merge_intervals(intervals):
    if not intervals: return []
    intervals.sort()
    merged = []
    curr_s, curr_e = intervals[0]
    for next_s, next_e in intervals[1:]:
        if next_s <= curr_e:
            curr_e = max(curr_e, next_e)
        else:
            merged.append((curr_s, curr_e))
            curr_s, curr_e = next_s, next_e
    merged.append((curr_s, curr_e))
    return merged

print("Optimizing intervals...")
for c in blacklist_dict: blacklist_dict[c] = merge_intervals(blacklist_dict[c])
for c in repeats_dict: repeats_dict[c] = merge_intervals(repeats_dict[c])

# ==========================================
#           GENERATION LOGIC
# ==========================================

output_dir = BASE_DIR / "sim"
output_dir.mkdir(parents=True, exist_ok=True)

# 1. Full Body Check (Blacklist & Internal Overlaps)
def check_overlap_fast(chrom, start, end, interval_dict):
    if chrom not in interval_dict: return False
    intervals = interval_dict[chrom]
    idx = bisect.bisect_right(intervals, (end, -1))
    if idx > 0:
        if intervals[idx - 1][1] > start:
            return True
    return False

# 2. Endpoints Check + Buffer (Simple Repeats)
def check_endpoints_in_repeats(chrom, start, end, repeat_dict, buffer):
    if chrom not in repeat_dict: return False
    intervals = repeat_dict[chrom]
    
    # Check Start Point
    idx_s = bisect.bisect_right(intervals, (start + buffer, float('inf')))
    if idx_s > 0:
        if intervals[idx_s - 1][1] > (start - buffer):
            return True 
        
    # Check End Point
    idx_e = bisect.bisect_right(intervals, (end + buffer, float('inf')))
    if idx_e > 0:
        if intervals[idx_e - 1][1] > (end - buffer):
            return True 
        
    return False

print(f"Generating {NUM_SAMPLES} samples...")
padding_width = len(str(NUM_SAMPLES))

for i in range(1, NUM_SAMPLES + 1):
    file_number = str(i).zfill(padding_width)
    filename = output_dir / f"sim{file_number}.txt"
    placed_cnvs = {} 
    results = [] 

    for group_idx, (count_target, cnv_len) in enumerate(CNV_MATRIX):
        generated_count = 0
        safety_break = 0
        
        while generated_count < count_target:
            safety_break += 1
            if safety_break > 500000:
                print(f"Warning: Sample {i}, Group {group_idx} (Len {cnv_len}): Space exhausted. Skipping.")
                break

            if CHR_INPUT == "all": target_chrom = random.choice(list(genome_db.keys()))
            else: target_chrom = CHR_INPUT
            
            chrom_len = genome_db[target_chrom]
            if chrom_len <= cnv_len: continue 

            start_pos = random.randint(0, chrom_len - cnv_len)
            end_pos = start_pos + cnv_len
            
            # 3. STRICT CHECK: Blacklist
            if BLACKLIST_FLAG == 1:
                if check_overlap_fast(target_chrom, start_pos, end_pos, blacklist_dict):
                    continue
            
            # 4. LOOSE CHECK + BUFFER: Simple Repeats
            if AVOID_REPEATS == 1:
                if check_endpoints_in_repeats(target_chrom, start_pos, end_pos, repeats_dict, REPEAT_BUFFER):
                    continue

            # 5. Internal Overlap
            if check_overlap_fast(target_chrom, start_pos, end_pos, placed_cnvs):
                continue

            if target_chrom not in placed_cnvs: placed_cnvs[target_chrom] = []
            bisect.insort(placed_cnvs[target_chrom], (start_pos, end_pos))

            cnv_type = random.choice(["DEL", "DUP"])
            if cnv_type == "DEL": depth = DELPEPTH
            else: depth = int(DUPDEPTH)
            
            output_chrom = target_chrom
            if REF_VERSION == "hg19": output_chrom = output_chrom.replace("chr", "")
            
            results.append(f"{output_chrom}\t{start_pos}\t{end_pos}\t{cnv_type}\t{depth}\n")
            generated_count += 1

    # Writing with UTF-8 and Unix newlines
    with open(filename, 'w', encoding='utf-8', newline='\n') as f:
        f.write("#Chrom\tStart\tEnd\tType\tDepth\n")
        for line in results:
            f.write(line)

print(f"Generation successful. Files located in: {output_dir}")
