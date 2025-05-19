#!/bin/bash

# This script takes a VCF line as input and appends the ENF score to the INFO column.
# It uses the fake-enformer tool to predict the score based on the VCF line.
# echo lines are for debugging purposes and can be uncommented if needed.

# Takes one full VCF line as input
line="$1"
#echo "Processing line: $line" >&2


# Use AWK to extract columns (tab-separated)
# cut -fN extracts the Nth column of a tab-separated string.
chrom=$(echo "$line" | cut -f1)
pos=$(echo "$line" | cut -f2)
ref=$(echo "$line" | cut -f4)
alt=$(echo "$line" | cut -f5)
#echo "Chromosome: $chrom" >&2
#echo "Position: $pos" >&2
#echo "Reference: $ref" >&2
#echo "Alternate: $alt" >&2

# Construct fake-enformer input string (hg38:chr9:127749702:C:T)
query="hg38:${chrom}:${pos}:${ref}:${alt}"
#echo "Query: $query" >&2


#Runs the scoring command using the fake-enformer and saves the numeric result (e.g. 0.8238475) in the variable score.
score=$(fake-enformer "$query")
#echo "Score: $score" >&2

# Passes two variables into awk:
# - line: The original VCF line
# - score: The score obtained from the fake-enformer prediction

# BFS=OFS="\t"  → Sets the field separator for input and output to tab
# split(line, fields, FS) → Splits the input line into an array called fields using the tab character as a delimiter
# fields[8] = fields[8] ";ENF=" score → Appends the ENF score to the 8th field (INFO column) of the VCF line
# for (i = 1; i <= length(fields); i++) printf "%s%s", fields[i], (i<length(fields)?OFS:ORS) 
# → Loops through the fields array and prints each field, separated by tabs

echo | awk -v line="$line" -v score="$score" '
BEGIN { 
  FS=OFS="\t"
  split(line, fields, FS)
  fields[8] = fields[8] ";ENF=" score
  for (i = 1; i <= length(fields); i++) 
    printf "%s%s", fields[i], (i<length(fields)?OFS:ORS)
  exit
}
'
