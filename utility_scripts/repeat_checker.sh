#!/bin/bash

# ==========================================
#           USER CONFIGURATION
# ==========================================

# This script checks if bamsurgeon's reported failure to insert synthetic CNV was caused by tandem repeats near CNV endpoints.

# ==========================================
#           USER CONFIGURATION
# ==========================================

# Reminder: UCSC files use "chr22" type of format
CHROM="chr10"
START=102532375
END=102537375

# Put this target file to the working directory manually, if you don't have it, here are the links:
# hg19:
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz
# hg38:
# http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz
TARGET_FILE="simpleRepeat.txt.gz"

# The buffer size determines how far from repeating regions can endpoints of synthetic CNV be placed
BUFFER=300

# ==========================================
#           WORKING DEPARTMENT
# ==========================================

# Check if file exists
if [ ! -f "$TARGET_FILE" ]; then
    echo "Error: File '$TARGET_FILE' not found in current directory."
    exit 1
fi

echo "Searching for simple repeats within ${BUFFER}bp of endpoints for $CHROM:$START-$END..."
echo "---------------------------------------------------------------"

# Search using awk
# UCSC simpleRepeat format: bin(1) chrom(2) start(3) end(4) ...
zcat "$TARGET_FILE" | awk -v c="$CHROM" -v s="$START" -v e="$END" -v buf="$BUFFER" '
    BEGIN { hits=0 }
    
    # Filter for chromosome matches first to save speed
    $2 == c {
        rep_start = $3
        rep_end = $4
        
        # --- Check START Point ---
        # We check if the repeat overlaps the interval [START - BUFFER, START + BUFFER]
        # Overlap Logic: (QueryStart < RepEnd) AND (QueryEnd > RepStart)
        if ((s - buf) < rep_end && (s + buf) > rep_start) {
            print "HIT NEAR START: " $0
            hits++
        }

        # --- Check END Point ---
        # We check if the repeat overlaps the interval [END - BUFFER, END + BUFFER]
        if ((e - buf) < rep_end && (e + buf) > rep_start) {
            print "HIT NEAR END:   " $0
            hits++
        }
    }
    
    END {
        print "---------------------------------------------------------------"
        if (hits > 0) 
            print "FAIL: Found " hits " simple repeats interfering with endpoints."
        else 
            print "PASS: No simple repeats found within " buf "bp of endpoints."
    }
'
