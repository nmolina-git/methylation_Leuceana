## Purpose: Calculate RPK 
# File: sorted_gene_length.txt, contains gene ids and associated length
# File: reversestr_counts.txt, contains gene id and the counts
# Load necessary libraries
library(dplyr)
library(readr)

# Files needs
gene_lengths_file <- "sorted_gene_length.txt"
counts_file <- "reversestr_counts.txt"

gene_lengths <- read_tsv(gene_lengths_file, col_names = c("GeneID", "Length"))

# Read in the date
counts <- read_delim(counts_file, delim = " ", col_names = c("GeneID", "Count"))
gene_lengths <- read_delim(gene_lengths_file, delim = " ", col_names = c("GeneID", "Length"))

# Check files
head(counts)
head(gene_lengths)

# Merge counts with gene lengths
merged_data <- inner_join(counts, gene_lengths, by = "GeneID")

# Calculate RPK
merged_data <- merged_data %>%
  mutate(RPK = Count / (Length / 1000)/1000000)

# Sort by GeneID
sorted_data <- merged_data %>%
  arrange(GeneID)

# Save the sorted results
write_csv(sorted_data, "RPK_scaled_reverse.csv")

