#!/bin/bash

# Check if genomes_dir argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 /path/to/genomes_dir"
    exit 1
fi

GENOMES_DIR=$1
PRE_ALIGN_DIR="${GENOMES_DIR}/records/compleasm/alignments/01_pre_alignments"
LOG_FILE="${GENOMES_DIR}/records/records_kc-align.txt"

# Check if the pre-alignment directory exists
if [ ! -d "$PRE_ALIGN_DIR" ]; then
    echo "Error: Directory $PRE_ALIGN_DIR not found."
    exit 1
fi

# File presence check and initialization
if [ ! -f "$LOG_FILE" ]; then
    echo "Creating new log file: $LOG_FILE"
    echo "KC-Align Pipeline Log Created: $(date)" > "$LOG_FILE"
else
    echo "Appending to existing log file: $LOG_FILE"
    echo -e "\n--- NEW RUN STARTED: $(date) ---" >> "$LOG_FILE"
fi

echo "Starting KC-Align. Logging output to: $LOG_FILE"

# Navigate to the pre-alignments directory
cd "$PRE_ALIGN_DIR" || exit

# Loop through every gene_id directory
for d in */; do
    # Enter the gene directory
    pushd "$d" > /dev/null || continue
    
    GENE_NAME=${d%/}
    
    # Calculate number of species
    # Count '>' in oth_* files, then add 1 for the ref_* file
    OTH_COUNT=$(grep -c "^>" oth_* 2>/dev/null || echo 0)
    TOTAL_SPECIES=$((OTH_COUNT + 1))
    
    # Terminal output
    echo "Processing $GENE_NAME ($TOTAL_SPECIES species total)"
    
    # Log Entry
    echo "DATE/TIME: $(date)" >> "$LOG_FILE"
    echo "GENE ID: $GENE_NAME" >> "$LOG_FILE"
    echo "SPECIES COUNT: $TOTAL_SPECIES" >> "$LOG_FILE"
    echo "COMMAND: kc-align -m gene -r ref_* -S oth_*" >> "$LOG_FILE"
    
    # Run kc-align and redirect BOTH stdout and stderr
    kc-align -m gene -r ref_* -S oth_* >> "$LOG_FILE" 2>&1
    
    echo "----------------------------------------------------------" >> "$LOG_FILE"
    
    # Return to the parent directory
    popd > /dev/null
done

echo "KC-Align processing complete. Check log at: $LOG_FILE"