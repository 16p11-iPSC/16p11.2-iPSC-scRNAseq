# Install necessary packages if not already installed
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("readxl")) install.packages("readxl")
if (!require("dplyr")) install.packages("dplyr")

library(ggplot2)
library(readxl)
library(dplyr)

# 1. Read the data
file_name <- "ASDandBMIbeta_withlog2FC_20012026.xlsx"
data <- read_excel(file_name)

# 2. Define Significance Thresholds
padj_cutoff <- 0.05
log2fc_cutoff <- 1

# ==========================================
# Figure 1: ASD (BETA vs log2FC)
# ==========================================

# Filter data for ASD related rows (ASD specific + Shared)
asd_data <- data %>%
  filter(Condition %in% c("ASD", "Shared")) %>%
  mutate(
    Significance = case_when(
      ASD_padj < padj_cutoff & (ASD_log2FoldChange > log2fc_cutoff | ASD_log2FoldChange < -log2fc_cutoff) ~ "Significant",
      TRUE ~ "Not Significant"
    )
  )

# Plot ASD Figure
asd_plot <- ggplot(asd_data, aes(x = ASD_BETA, y = ASD_log2FoldChange, color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  geom_hline(yintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "black", alpha=0.5) +
  labs(
    title = "ASD: BETA Value vs log2FoldChange",
    subtitle = paste("Colored significant if padj <", padj_cutoff, "& |log2FC| >", log2fc_cutoff),
    x = "ASD BETA",
    y = "ASD log2FoldChange"
  ) +
  theme_minimal()

# ==========================================
# Figure 2: BMI (BETA vs log2FC)
# ==========================================

# Filter data for BMI related rows (BMI specific + Shared)
bmi_data <- data %>%
  filter(Condition %in% c("BMI", "Shared")) %>%
  mutate(
    Significance = case_when(
      BMI_padj < padj_cutoff & (BMI_log2FoldChange > log2fc_cutoff | BMI_log2FoldChange < -log2fc_cutoff) ~ "Significant",
      TRUE ~ "Not Significant"
    )
  )

# Plot BMI Figure
# Note: Using 'BMI_Mean_BETA' based on file header structure
bmi_plot <- ggplot(bmi_data, aes(x = BMI_Mean_BETA, y = BMI_log2FoldChange, color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "blue")) +
  geom_hline(yintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "black", alpha=0.5) +
  labs(
    title = "BMI: BETA Value vs log2FoldChange",
    subtitle = paste("Colored significant if padj <", padj_cutoff, "& |log2FC| >", log2fc_cutoff),
    x = "BMI Mean BETA",
    y = "BMI log2FoldChange"
  ) +
  theme_minimal()

# ==========================================
# Save or Display Plots
# ==========================================

# Display plots
print(asd_plot)
print(bmi_plot)


ggsave("ASD_BETA_vs_Log2FC.png", plot = asd_plot, width = 6, height = 5)
ggsave("BMI_BETA_vs_Log2FC.png", plot = bmi_plot, width = 6, height = 5)