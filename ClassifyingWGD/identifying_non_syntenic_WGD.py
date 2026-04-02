# Purpose:
# This script filters the WGD gene pairs identified by DupGenFinder into a file where it doesnt
# meet WGD expectation or syntenic secondary expectation 
#
# Only gene pairs are retained if:
# 1. They are not WGD expectation or syntenic secondary expectation 
# 2. Each gene matches only one other gene.
#


# Input file:
#    Ltr.wgd.pairs
#    Ltr.wgd.pairs.filtered
#    Ltr.wgd.pairs.samechr

# Output file:
#    Ltr.wgd.pairs.otherchr
from collections import Counter


# Normalize a gene pair so that order doesn't matter

def normalize_pair(g1, g2):
    return tuple(sorted([g1, g2]))

# Load all pairs 
def load_pairs(filename):
    pairs = []
    with open(filename) as f:
        for line in f:
            if not line.strip():
                continue
            g1, l1, g2, l2, *_ = line.strip().split("\t")
            pairs.append((g1, l1, g2, l2))
    return pairs

# Load pairs into a set 
def load_pair_set(filename):
    s = set()
    with open(filename) as f:
        for line in f:
            if not line.strip():
                continue
            g1, _, g2, _, *_ = line.strip().split("\t")
            s.add(normalize_pair(g1, g2))
    return s

# Filter for one-to-one gene pairs only
# Keeps only pairs where each gene appears exactly once
def one_to_one(pairs):
    c1 = Counter(g1 for g1, _, _, _ in pairs)
    c2 = Counter(g2 for _, _, g2, _ in pairs)

    return [
        (g1, l1, g2, l2)
        for g1, l1, g2, l2 in pairs
        if c1[g1] == 1 and c2[g2] == 1
    ]

# -----------------------------------------------------------------------------
# MAIN 
# -----------------------------------------------------------------------------

if __name__ == "__main__":

    # Input files
    all_file = "Ltr.wgd.pairs"          # all pairs from DupGenFinder
    dup_file = "Ltr.wgd.pairs.filtered" # file created from filtering WGD
    same_file = "Ltr.wgd.pairs.samechr" # file created from filtering syntenic secondary duplication

    # Output file for other chromosome pairs
    out_other = "Ltr.wgd.pairs.otherchr"

    # Load all pairs
    all_pairs = load_pairs(all_file)

    # Load exclusion sets 
    dup_pairs = load_pair_set(dup_file)
    same_pairs = load_pair_set(same_file)

    # Keep only pairs not already classified as WGD or syntenic secondary
    other = [
        (g1, l1, g2, l2)
        for g1, l1, g2, l2 in all_pairs
        if normalize_pair(g1, g2) not in dup_pairs
        and normalize_pair(g1, g2) not in same_pairs
    ]

    # Apply one-to-one rule
    other = one_to_one(other)

    # Write filtered "other chromosome" pairs to output file
    with open(out_other, "w") as out:
        for g1, l1, g2, l2 in other:
            out.write(f"{g1}\t{l1}\t{g2}\t{l2}\n")

    # Print summary
    print("Done!")
    print(f"Other-chromosome pairs written to: {out_other}")
    print(f"Count: {len(other)}")