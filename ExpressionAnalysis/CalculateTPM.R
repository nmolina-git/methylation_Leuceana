## Purpose: Calculate TPM for genes in Leuceana trichandra
# File: "RPK_scaled_reverse.csv", has gene ID, and associated RPK value
# Load necessary libraries
library(dplyr)
library(readr)


# Read the RPK data
rpk_data <- read_tsv("RPK_scaled_reverse.csv")

# Calculate the total RPK
total_rpk <- sum(rpk_data$RPK, na.rm = TRUE)

# Calculate TPM
rpk_data <- rpk_data %>%
  mutate(TPM = (RPK / total_rpk) * 1e6)    


# Sort the genes
sorted_tpm <- rpk_data %>%
  arrange(GeneID)   

# Save the sorted TPM results

write_csv(sorted_tpm, "TPM_reverse.csv")



