#!/bin/bash

# Arguments
INPUT_VCF="$1" #Assigns INPUT_VCF to the first argument passed to the script
OUTPUT_TSV="$2" # Assigns OUTPUT_TSV to the second argument passed to the script


# Extract non-header lines (SNP entries)
grep -v "^#" "$INPUT_VCF" > "$OUTPUT_TSV"

# Explanation:
# -v: Inverts the match, so it selects lines that do not match the pattern.
# "^#": This pattern matches lines that start with "#", which are typically header lines in VCF files.
# The grep command is used to filter lines from the input VCF file.
# The output is redirected to the specified output TSV file.



