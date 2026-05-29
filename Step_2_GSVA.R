
suppressMessages(library(fgsea))
suppressMessages(library(gage))
suppressMessages(library(stringr))
suppressMessages(library(tibble))

# 1. Prepare the Ranked Gene List from DESeq2 Output (res_df from Step 5)
# Filter out NA log2FoldChanges and establish the ranked list
res_gsea <- res_df[!is.na(res_df$log2FoldChange), ]

# Create the named vector required for GSEA
gene_list <- res_gsea$log2FoldChange
names(gene_list) <- res_gsea$Gene_Symbol

# Sort decreasingly and remove duplicates
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!duplicated(names(gene_list))]

# 2. Load the MSigDB C2 pathways
# Make sure "c2.all.v7.5.1.symbols.gmt" is in your working directory
GO_file <- "./c2.all.v7.5.1.symbols.gmt"
myGO <- fgsea::gmtPathways(GO_file)
pval_threshold <- 0.05

# ---------------------------------------------------------
# A. Run Method 1: FGSEA
# ---------------------------------------------------------
set.seed(123) # For reproducibility
fgRes <- fgsea::fgsea(pathways = myGO, 
                      stats = gene_list,
                      minSize = 15,
                      maxSize = 600) %>% 
  as.data.frame() %>% 
  dplyr::filter(padj < pval_threshold)

# ---------------------------------------------------------
# B. Run Method 2: GAGE
# ---------------------------------------------------------
gaRes <- gage::gage(gene_list, gsets = myGO, same.dir = TRUE, set.size = c(15, 600))

ups <- as.data.frame(gaRes$greater) %>% 
  tibble::rownames_to_column("Pathway") %>% 
  dplyr::filter(!is.na(p.geomean) & q.val < pval_threshold) %>%
  dplyr::select("Pathway")

downs <- as.data.frame(gaRes$less) %>% 
  tibble::rownames_to_column("Pathway") %>% 
  dplyr::filter(!is.na(p.geomean) & q.val < pval_threshold) %>%
  dplyr::select("Pathway")

# ---------------------------------------------------------
# C. Intersect the Results (Keep only consensus pathways)
# ---------------------------------------------------------
keepups <- fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ]
keepdowns <- fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ]

fgRes_filtered <- fgRes[!is.na(match(fgRes$pathway, c(keepups$pathway, keepdowns$pathway))), ] %>% 
  arrange(desc(NES))

# Clean up pathway names and assign direction
fgRes_filtered$pathway <- stringr::str_replace(fgRes_filtered$pathway, "GO_", "")
fgRes_filtered$Enrichment <- ifelse(fgRes_filtered$NES > 0, "Up-regulated", "Down-regulated")

# Save the robust consensus pathways to CSV
write.csv(fgRes_filtered, "Robust_GSEA_Consensus_Pathways_Progenitors.csv", row.names = FALSE)

# ---------------------------------------------------------
# D. Visualization (Top 10 Up and Down)
# ---------------------------------------------------------
filtRes <- rbind(head(fgRes_filtered, n = 10), tail(fgRes_filtered, n = 10))

# Define colors
upcols <- colorRampPalette(colors = c("red4", "red1", "lightpink"))(sum(filtRes$Enrichment == "Up-regulated"))
downcols <- colorRampPalette(colors = c("lightblue", "blue1", "blue4"))(sum(filtRes$Enrichment == "Down-regulated"))
colos <- c(upcols, downcols)
names(colos) <- 1:length(colos)

filtRes$Index <- as.factor(1:nrow(filtRes))

# Generate the plot
gsea_plot <- ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = Index)) +
  scale_fill_manual(values = colos) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score", title = "GSEA - Consensus Pathways (fgsea + gage)") + 
  theme_minimal() + 
  theme(legend.position = "none")

ggsave("GSEA_Top_Pathways_Barplot.pdf", plot = gsea_plot, width = 12, height = 7)

cat("GSEA completed successfully. Consensus plot saved.\n")