## Purpose: Perform a statistical test to determine if methylation types differ in expression
# Note: methylation types : gbM, teM, uM
# File needed:geneID_meth_RPK.txt

# Load necessary libraries
library(FSA)
library(ggplot2)

# Load data
rpk_df <- read.csv("geneID_meth_RPK.txt", sep = "\t")
head(rpk_df)

# Add column names
colnames(rpk_df) <- c("GeneName", "MethylType", "RPK")

# Check the dataframe
head(rpk_df)

# Perform the Kruskal-Wallis test
kruskal_test_result <- kruskal.test(RPK ~ rpk_df$MethylType, data = rpk_df)

# Show the result
print(kruskal_test_result)

# Perform Dunn's test with Benjamini-Hochberg correction (adjusts for multiple comparisons)
dunn_result <- dunnTest(RPK ~ rpk_df$MethylType, data = rpk_df, method = "bh")

# View results
print(dunn_result)

# Create Violing plot
ggplot(rpk_df, aes(x = MethylType, y = RPK, fill = MethylType)) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.4, scale = "width") +
  geom_jitter(aes(color = MethylType), width = 0.15, size = 1, alpha = 0.5) +
  
  scale_fill_manual(values = c("gbm" = "#0072B2", "tem" = "#D55E00", "um" = "#009E73")) +
  scale_color_manual(values = c("gbm" = "#0072B2", "tem" = "#D55E00", "um" = "#009E73")) +
  scale_y_log10() +
  
  # Significance bars + asterisks
  geom_segment(aes(x = 1, xend = 2, y = 1300, yend = 1300)) +
  geom_text(aes(x = 1.5, y = 1400, label = "***")) +
  
  geom_segment(aes(x = 1, xend = 3, y = 1450, yend = 1450)) +
  geom_text(aes(x = 2, y = 1550, label = "***")) +
  
  geom_segment(aes(x = 2, xend = 3, y = 1600, yend = 1600)) +
  geom_text(aes(x = 2.5, y = 1700, label = "***")) +
  
  labs(
    title = "RPK by Methylation Type",
    x = "Methylation Type",
    y = "RPK (log10 scale)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

