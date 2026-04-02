# Purpose:
# This script filters the WGD gene pairs identified by DupGenFinder into syntenic secondary duplication.
#
# Only gene pairs are retained if:
# 1. Both genes occur on the same chromosome (e.g., If there both on Chromosome 1)
# 2. Each gene matches only one other gene.
#
# This removes cases where one gene matches multiple genes, keeping only
# unique one-to-one WGD pairs.

# Input file:
#    Ltr.wgd.pairs

# Output file:
#    Ltr.wgd.pairs.samechr
import csv
from collections import Counter

# -----------------------------------------------------------------------------
# Function to remove trailing version suffix from a gene ID (e.g., .1, .2)
# -----------------------------------------------------------------------------
def strip_version(gene):
    # Remove trailing .1, .2, etc.
    return gene.split('.')[0]

# -----------------------------------------------------------------------------
# Function to extract chromosome number from location string
# -----------------------------------------------------------------------------
def get_chr(location):
    contig = location.split(':')[0]   # e.g., 'Contig0'
    num = ''.join(ch for ch in contig if ch.isdigit())
    return int(num)

# -----------------------------------------------------------------------------
# Load gene pairs 
# -----------------------------------------------------------------------------
def load_pairs(filename):
    pairs = []
    with open(filename) as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # skip header

        for row in reader:
            g1, loc1, g2, loc2, _ = row
            g1 = strip_version(g1)
            g2 = strip_version(g2)
            pairs.append((g1, loc1, g2, loc2))
    return pairs

# -----------------------------------------------------------------------------
# Filter gene pairs to keep only those on the same chromosome in one-to-one mapping
# -----------------------------------------------------------------------------
def filter_same_chromosome(pairs):
    # Step 1 — keep only pairs where chromosome1 == chromosome2
    valid = []
    for g1, l1, g2, l2 in pairs:
        c1 = get_chr(l1)
        c2 = get_chr(l2)

        if c1 == c2:  # SAME chromosome only
            valid.append((g1, l1, g2, l2))

    # Step 2 — count occurrences of each gene
    c_g1 = Counter(g1 for g1, _, _, _ in valid)
    c_g2 = Counter(g2 for _, _, g2, _ in valid)

    # Step 3 — keep only one-to-one gene pairs
    final = [
        (g1, l1, g2, l2)
        for (g1, l1, g2, l2) in valid
        if c_g1[g1] == 1 and c_g2[g2] == 1
    ]

    return final

# -----------------------------------------------------------------------------
# Main script
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    input_file = "Ltr.wgd.pairs"           # input file
    output_file = "Ltr.wgd.pairs.samechr"  # output file after filtering

    # Load gene pairs
    pairs = load_pairs(input_file)

    # Filter for same-chromosome one-to-one pairs
    filtered = filter_same_chromosome(pairs)

    # Write filtered pairs to output file
    with open(output_file, "w") as out:
        for g1, l1, g2, l2 in filtered:
            out.write(f"{g1}\t{l1}\t{g2}\t{l2}\n")

    # Print completion message
    print(f"Done! Same-chromosome filtered output saved to: {output_file}")