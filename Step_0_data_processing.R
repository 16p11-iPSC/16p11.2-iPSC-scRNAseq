# ==============================================================================
# Project: 16p11.2 Microdeletion in Forebrain Interneuron Progenitors
# Description: Complete scRNA-seq analysis pipeline including data QC, 
#              Harmony integration, cell type annotation, in-vivo correlation, 
#              and pseudobulk differential expression analysis (DESeq2).
# ==============================================================================

# --- 0. Load Required Libraries ---
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(gplots)
  library(harmony)
  library(sctransform)
  library(DoubletFinder)
  library(DESeq2)
})

# Define global color palette
color_plane <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",  
                 "maroon", "gold1", "skyblue2", "palegreen2","#FDBF6F", 
                 "#8b8b00", "orchid1","black", "darkturquoise", "darkorange4", 
                 "brown", "gray70")

# ==============================================================================
# STEP 1: Data Cleaning & Quality Control
# ==============================================================================
# Define directories for your biological replicates (CON and DEL)
data_dirs <- c("path/to/CON_1", "path/to/CON_2", "path/to/CON_3", 
               "path/to/DEL_1", "path/to/DEL_2", "path/to/DEL_3")
sample_ids <- c("CON_1", "CON_2", "CON_3", "DEL_1", "DEL_2", "DEL_3")

seurat_list <- list()

for (i in 1:length(data_dirs)) {
  # Read raw 10X data
  counts <- Read10X(data.dir = data_dirs[i])
  obj <- CreateSeuratObject(counts = counts, project = sample_ids[i], min.cells = 60)
  
  # Calculate QC metrics
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  
  # Filter based on strict thresholds
  obj <- subset(obj, subset = nFeature_RNA > 1000 & 
                  nFeature_RNA < 10000 & 
                  nCount_RNA < 40000 & 
                  percent.mt < 10 & 
                  percent.ribo < 40)
  
  # Store condition and replicate metadata
  obj$Conditions <- ifelse(grepl("CON", sample_ids[i]), "CON", "DEL")
  obj$Replicate <- sample_ids[i]
  
  seurat_list[[i]] <- obj
}

# Merge all filtered objects
iPSC_16p11_Obj <- merge(seurat_list[[1]], y = seurat_list[2:length(seurat_list)], 
                        add.cell.ids = sample_ids, project = "16p11_scRNAseq")
rm(seurat_list) # Clear memory

# Remove sex-linked and mitochondrial genes before downstream analysis
sex.genes <- c('DDX3Y','EIF2S3Y','UTY','KDM5D','XIST','TSIX')
mt.genes <- rownames(iPSC_16p11_Obj)[grep("^MT-", rownames(iPSC_16p11_Obj))]
iPSC_16p11_Obj <- subset(iPSC_16p11_Obj, features = setdiff(rownames(iPSC_16p11_Obj), c(mt.genes, sex.genes)))

# ==============================================================================
# STEP 2: Normalization, Integration (Harmony), and Clustering
# ==============================================================================
iPSC_16p11_Obj <- NormalizeData(iPSC_16p11_Obj, scale.factor = 10000)
iPSC_16p11_Obj <- FindVariableFeatures(iPSC_16p11_Obj, selection.method = "vst")

# Cell cycle scoring and regression
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
iPSC_16p11_Obj <- CellCycleScoring(iPSC_16p11_Obj, s.features = s.genes, g2m.features = g2m.genes)
iPSC_16p11_Obj$CC.Difference <- iPSC_16p11_Obj$S.Score - iPSC_16p11_Obj$G2M.Score  
iPSC_16p11_Obj <- ScaleData(iPSC_16p11_Obj, vars.to.regress = 'CC.Difference')

# PCA and Harmony Batch Correction
iPSC_16p11_Obj <- RunPCA(iPSC_16p11_Obj, npcs = 30)
iPSC_16p11_Obj <- RunHarmony(iPSC_16p11_Obj, group.by.vars = "Replicate", dims.use = 1:30)

# Clustering and UMAP on Harmony embeddings
iPSC_16p11_Obj <- FindNeighbors(iPSC_16p11_Obj, reduction = "harmony", dims = 1:30)
iPSC_16p11_Obj <- FindClusters(iPSC_16p11_Obj, resolution = 0.5)
iPSC_16p11_Obj <- RunUMAP(iPSC_16p11_Obj, reduction = "harmony", dims = 1:30)

# Note: Manual cluster annotation is performed here. For the purpose of this script,
# we assume metadata column "broad_Type" is populated with cell type annotations.

# ==============================================================================
# STEP 3: Visualizations (UMAPs, Barplots, DotPlots, Heatmaps)
# ==============================================================================
Idents(iPSC_16p11_Obj) <- "broad_Type"

# 1. UMAPs and Cell Proportions
# ====================================================================
# PART 1: The Overview UMAPs (Side-by-Side)
# ====================================================================
# A. UMAP by Genotype ("Conditions")
plot_geno <- DimPlot(seurat_obj, reduction = "umap", group.by = "Conditions", 
                     cols = c("Controls" = "#4C84B9", "Deletions" = "#D1605E"), pt.size = 0.5) + 
             ggtitle("Genotype")

# B. UMAP by Cell Type ("broad_Type")
plot_cell <- DimPlot(seurat_obj, reduction = "umap", group.by = "broad_Type", 
                     label = TRUE, pt.size = 0.5) + NoLegend() + 
             ggtitle("Cell Types")

# C. UMAP by Replicate ("Replicate")
plot_rep <- DimPlot(seurat_obj, reduction = "umap", group.by = "Replicate", pt.size = 0.5) + 
            ggtitle("Replicates")

# Save the 3-panel plot
overview_plots <- plot_geno + plot_cell + plot_rep
ggsave("UMAP_Overview_Geno_CellType_Rep.pdf", plot = overview_plots, width = 18, height = 5)

# ====================================================================
# PART 2: Split UMAPs (Showing cell types separated by Replicate)
# ====================================================================
# This creates a separate UMAP for each replicate, colored by cell type.
plot_split <- DimPlot(seurat_obj, reduction = "umap", group.by = "broad_Type", 
                      split.by = "Replicate", ncol = 3, pt.size = 0.5) +
              ggtitle("Cell Types across Replicates")

ggsave("UMAP_Split_by_Replicate.pdf", plot = plot_split, width = 12, height = 7)

# ====================================================================
# PART 3: Cell Proportion Stacked Bar Plot 
# ====================================================================
# Calculate the proportions using your specific column names
prop_data <- seurat_obj@meta.data %>%
  group_by(Conditions, Replicate, broad_Type) %>%
  summarise(CellCount = n(), .groups = 'drop') %>%
  group_by(Replicate) %>%
  mutate(Percent = CellCount / sum(CellCount) * 100)

# Create the Stacked Bar Plot
bar_plot <- ggplot(prop_data, aes(x = Replicate, y = Percent, fill = broad_Type)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
  facet_grid(~ Conditions, scales = "free_x", space = "free_x") + # Groups bars by CON/DEL
  theme_classic() +
  labs(title = "Cell Type Composition per Replicate",
       x = "Replicate", y = "Percentage of Cells (%)", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "lightgrey"))

ggsave("Barplot_Cell_Proportions_by_Replicate.pdf", plot = bar_plot, width = 8, height = 6)

# 2. test only: Average Expression Heatmap
iPSC_16p11_Obj$Unique_ID <- paste(iPSC_16p11_Obj$broad_Type, iPSC_16p11_Obj$Replicate, sep = "_")
avg_expr <- as.data.frame(as.matrix(AverageExpression(iPSC_16p11_Obj, group.by = "Unique_ID")$RNA))
target_genes_expr <- avg_expr[as.character(c("VIM","HES1","CLU","LYPD1","FABP7","SFRP1","HMGB2","ASCL1","INSM1","GAD2","DLX2","LHX6","SOX6","SST")), ]

pheatmap(target_genes_expr[rowMeans(target_genes_expr) > 0, ], 
         clustering_distance_rows = "euclidean", clustering_method = "ward.D2",
         scale = "row", cluster_rows = FALSE, cluster_cols = FALSE, col = bluered(100), 
         main = "Figure 1A3: In Vitro Heatmap")

Progenitors <- subset(iPSC_16p11_Obj, idents = "Progenitor")

# Create Pseudobulk matrix by summing counts across replicates
pseudo_counts <- matrix(nrow = nrow(Progenitors), ncol = length(unique(Progenitors$Replicate)))
rownames(pseudo_counts) <- rownames(Progenitors)
colnames(pseudo_counts) <- unique(Progenitors$Replicate)

for (rep in colnames(pseudo_counts)) {
  cells <- colnames(Progenitors)[Progenitors$Replicate == rep]
  pseudo_counts[, rep] <- rowSums(Progenitors@assays$RNA@counts[, cells])
}

# Clean column names (converting hyphens to underscores for clean formatting)
colnames(pseudo_counts) <- gsub("-", "_", colnames(pseudo_counts))

# Setup DESeq2 metadata
colData <- data.frame(
  Replicate = colnames(pseudo_counts),
  Condition = ifelse(grepl("CON", colnames(pseudo_counts)), "CON", "DEL")
)
rownames(colData) <- colData$Replicate

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = pseudo_counts, colData = colData, design = ~ Condition)
dds$Condition <- relevel(dds$Condition, ref = "CON")
dds <- DESeq(dds)
res <- results(dds, contrast = c("Condition", "DEL", "CON"))

# Format Output (Ensuring Gene_Symbol is mapped uniquely to rownames)
res_df <- as.data.frame(res)
res_df$Gene_Symbol <- rownames(res_df)
res_df <- res_df[!duplicated(res_df$Gene_Symbol), ] # Ensure uniqueness
rownames(res_df) <- res_df$Gene_Symbol
res_df <- res_df[, c("Gene_Symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

# Filter significant DEGs (padj < 0.05 and |log2FC| > 0.25)
sig_DEGs <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.25, ]

# Save Results
write.csv(sig_DEGs, "Significant_DEGs_Progenitors_DEL_vs_CON.csv", row.names = FALSE)
write.csv(res_df, "All_Genes_Progenitors_DEL_vs_CON.csv", row.names = FALSE)

cat("Pipeline completed successfully. Results saved to working directory.\n")