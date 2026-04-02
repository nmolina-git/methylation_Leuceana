# Purpose:
# This script filters the WGD gene pairs identified by DupGenFinder.
#
# Only gene pairs are retained if:
# 1. Both genes occur on their expected duplicated chromosome pair.
# 2. Each gene matches only one other gene.
#
# This removes cases where one gene matches multiple genes, keeping only
# unique one-to-one WGD pairs.

# Input file:
#    Ltr.wgd.pairs

# Output file:
#    Ltr.wgd.pairs.filtered

import csv
from collections import Counter

# WGD duplicate chromosome pairs 
WGD_PAIRS = {
    (0, 21), (1, 25), (10, 23), (11, 27), (12, 14), (13, 24),
    (15, 22), (16, 17), (18, 19), (20, 26), (2, 9), (3, 7),
    (4, 6), (5, 8)
}

# make bidirectional so both (a,b) and (b,a) are valid
WGD_PAIRS.update({(b,a) for (a,b) in WGD_PAIRS})

# -----------------------------------------------------------------------------
# Function to remove version suffix from a gene ID (e.g., .1 or .2)
# -----------------------------------------------------------------------------
def strip_version(gene):
    return gene.split('.')[0]

# -----------------------------------------------------------------------------
# Function to extract chromosome number from gene ID like Leutr000g024000
# -----------------------------------------------------------------------------
def get_chr_from_gene(gene):
    gene = strip_version(gene)
    start = 5  # after 'Leutr'
    end = gene.index('g')
    return int(gene[start:end])

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
            g1_clean = strip_version(g1)
            g2_clean = strip_version(g2)
            pairs.append((g1_clean, loc1, g2_clean, loc2))
    return pairs

# -----------------------------------------------------------------------------
# Filter gene pairs to retain only WGD one-to-one relationships
# Steps:
# 1. Keep only pairs on known WGD chromosomes
# 2. Count occurrences of each gene
# 3. Keep only genes that appear exactly once (one-to-one)
# -----------------------------------------------------------------------------
def filter_wgd_pairs(pairs):
    # Step 1 — keep only WGD chromosome pairs
    valid = []
    for g1, loc1, g2, loc2 in pairs:
        c1 = get_chr_from_gene(g1)
        c2 = get_chr_from_gene(g2)
        if (c1, c2) in WGD_PAIRS:
            valid.append((g1, loc1, g2, loc2))

    # Step 2 — count occurrences for one-to-one filtering
    c_g1 = Counter(g1 for g1, _, _, _ in valid)
    c_g2 = Counter(g2 for _, _, g2, _ in valid)

    # Step 3 — keep only one-to-one pairs
    final = [
        (g1, loc1, g2, loc2)
        for (g1, loc1, g2, loc2) in valid
        if c_g1[g1] == 1 and c_g2[g2] == 1
    ]
    return final

# -----------------------------------------------------------------------------
# Main script
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    input_file = "Ltr.wgd.pairs"            # input file from DupGenFinder
    output_file = "Ltr.wgd.pairs.filtered"  # output filtered file

    # Load all gene pairs
    pairs = load_pairs(input_file)

    # Filter for one-to-one WGD gene pairs
    filtered = filter_wgd_pairs(pairs)

    # Write filtered pairs to output file
    with open(output_file, "w") as out:
        for g1, loc1, g2, loc2 in filtered:
            out.write(f"{g1}\t{loc1}\t{g2}\t{loc2}\n")

    # Print completion message
    print(f"Done! Filtered WGD pairs saved to: {output_file}")