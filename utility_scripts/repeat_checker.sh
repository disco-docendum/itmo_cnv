#!/bin/bash

# Configuration

# Reminder: UCSC files use "chr22" type of format
TARGET_FILE="simpleRepeat.txt.gz"
CHROM="chr22"
START=33547491
END=33552491

# Check if file exists, if you need it, just download manually,
# hg19:
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz
# hg38:
# http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz
if [ ! -f "$TARGET_FILE" ]; then
    echo "Error: File '$TARGET_FILE' not found in current directory."
    exit 1
fi

echo "Searching for overlaps with $CHROM:$START-$END in $TARGET_FILE..."
echo "---------------------------------------------------------------"

# Search using awk
# UCSC simpleRepeat format: bin(1) chrom(2) start(3) end(4) ...
# Overlap logic: (TargetStart < DbEnd) AND (TargetEnd > DbStart)
zcat "$TARGET_FILE" | awk -v c="$CHROM" -v s="$START" -v e="$END" '
    BEGIN { hits=0 }
    
    # Filter for chromosome matches first to save speed
    $2 == c {
        # Check overlap
        if (s < $4 && e > $3) {
            print "HIT FOUND: " $0
            hits++
        }
    }
    
    END {
        print "---------------------------------------------------------------"
        if (hits > 0) 
            print "Total Overlaps Found: " hits
        else 
            print "No overlaps found. This region is NOT in the simple repeats file."
    }
'
