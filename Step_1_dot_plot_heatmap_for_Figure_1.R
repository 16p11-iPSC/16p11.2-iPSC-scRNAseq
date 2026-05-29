# ==============================================================================
# Script Name: Figure_1_and_Reviewer_Response.R
# Description: Unified pipeline for in-vitro cell type annotation, reviewer-
#              requested UMAPs/DotPlots, and in-vivo transcriptomic correlation.
# ==============================================================================

# --- 1. Load Required Libraries ---
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(gplots)
  library(harmony)
  library(sctransform)
})

# Define global color palette
color_plane <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",  
                 "maroon", "gold1", "skyblue2", "palegreen2","#FDBF6F", 
                 "#8b8b00", "orchid1","black", "darkturquoise", "darkorange4", 
                 "brown", "gray70")

# ==============================================================================
# PART 1: In-Vitro Characterization (Reviewer Plots & Figure 1A3)
# ==============================================================================

# --- Load the primary in-vitro object ---
# Ensure the file path matches your local directory
iPSC_16p11_Obj <- readRDS("./filtered_object_19th_Jan_2025.rds") 
Idents(iPSC_16p11_Obj) <- "broad_Type"

# ---------------------------------------------------------
# A. Reviewer UMAPs (Genotype, Cell Type, Replicate)
# ---------------------------------------------------------
plot_geno <- DimPlot(iPSC_16p11_Obj, reduction = "umap", group.by = "Conditions", 
                     cols = c("CON" = "#4C84B9", "DEL" = "#D1605E"), pt.size = 0.5) + 
             ggtitle("Genotype")

plot_cell <- DimPlot(iPSC_16p11_Obj, reduction = "umap", group.by = "broad_Type", 
                     label = TRUE, pt.size = 0.5) + NoLegend() + 
             ggtitle("Cell Types")

plot_rep <- DimPlot(iPSC_16p11_Obj, reduction = "umap", group.by = "Replicate", pt.size = 0.5) + 
            ggtitle("Replicates")

overview_plots <- plot_geno + plot_cell + plot_rep
ggsave("UMAP_Overview_Geno_CellType_Rep.pdf", plot = overview_plots, width = 18, height = 5)

plot_split <- DimPlot(iPSC_16p11_Obj, reduction = "umap", group.by = "broad_Type", 
                      split.by = "Replicate", ncol = 3, pt.size = 0.5) +
              ggtitle("Cell Types across Replicates")
ggsave("UMAP_Split_by_Replicate.pdf", plot = plot_split, width = 12, height = 4)

# ---------------------------------------------------------
# B. Reviewer Bar Plot (Cell Proportions)
# ---------------------------------------------------------
prop_data <- iPSC_16p11_Obj@meta.data %>%
  group_by(Conditions, Replicate, broad_Type) %>%
  summarise(CellCount = n(), .groups = 'drop') %>%
  group_by(Replicate) %>%
  mutate(Percent = CellCount / sum(CellCount) * 100)

bar_plot <- ggplot(prop_data, aes(x = Replicate, y = Percent, fill = broad_Type)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
  facet_grid(~ Conditions, scales = "free_x", space = "free_x") +
  theme_classic() +
  labs(title = "Cell Type Composition per Replicate", x = "Replicate", y = "Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_rect(fill = "lightgrey"))
ggsave("Barplot_Cell_Proportions_by_Replicate.pdf", plot = bar_plot, width = 8, height = 6)

# ---------------------------------------------------------
# C. Reviewer Dot Plots (Identity Verification)
# ---------------------------------------------------------
marker_genes <- c("FOXG1", "DLX2", "MKI67", "TUBB3", "VIM", "SOX6", "GAD2") 

dot_plot_standard <- DotPlot(iPSC_16p11_Obj, features = marker_genes, group.by = "broad_Type", dot.scale = 8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  labs(title = "Marker Gene Expression by Cell Type") +
  scale_color_gradient(low = "lightgrey", high = "darkred")
ggsave("New_DotPlot_CellTypes.pdf", plot = dot_plot_standard, width = 8, height = 5)

dot_plot_split <- DotPlot(iPSC_16p11_Obj, features = marker_genes, group.by = "broad_Type", 
                          split.by = "Conditions", cols = c("#4C84B9", "#D1605E"), dot.scale = 8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  labs(title = "Marker Gene Expression (CON vs DEL)")
ggsave("New_DotPlot_Split_by_Genotype.pdf", plot = dot_plot_split, width = 10, height = 6)

# ---------------------------------------------------------
# D. Figure 1A3 (In Vitro Marker Heatmap)
# ---------------------------------------------------------
iPSC_16p11_Obj$Unique_ID <- paste(iPSC_16p11_Obj$broad_Type, iPSC_16p11_Obj$Replicate, sep = "_")
averages_seurat <- AverageExpression(iPSC_16p11_Obj, return.seurat = TRUE, group.by = "Unique_ID")
averages_expr <- as.data.frame(as.matrix(averages_seurat@assays$RNA@data))

cluster_levels <- levels(as.factor(iPSC_16p11_Obj$Unique_ID))
annotation_col <- data.frame(CellType = cluster_levels)
rownames(annotation_col) <- cluster_levels
ann_colors <- list(CellType = setNames(color_plane[1:length(cluster_levels)], cluster_levels))

fig1_markers <- c("FOXG1","VIM","CLU","HES1","SOX2","FABP7","INSM1","ZFXH3","NKX2-1",
              "OLIG2","GSX2","ASCL1","DLX1","DLX2","LHX6","NR2F2","GAD1","GAD2","VGAT","SLC32A1",
			  "SST","PVALB","EMX1","EMX2","NEUROG2","PAX6","TBR2","TBR1")
target_genes_expr <- averages_expr[as.character(fig1_markers), cluster_levels]
target_genes_expr <- target_genes_expr[rowMeans(target_genes_expr) > 0, ]

pheatmap(target_genes_expr, clustering_distance_rows = "euclidean", clustering_method = "ward.D2",
         scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col = annotation_col, annotation_colors = ann_colors,
         col = bluered(100), breaks = seq(-1, 1, length.out = 101), 
         show_rownames = TRUE, show_colnames = FALSE, main = "Figure 1A3: In Vitro Heatmap")


# ==============================================================================
# PART 2: In-Vivo Processing & Comparison (Figure 1A4)
# ==============================================================================

# ---------------------------------------------------------
# A. Process Yu et al. (GSE165388) Reference
# ---------------------------------------------------------
sample_names <- c('./9W','./10W','./11W','./12W')
project.data <- lapply(seq_along(sample_names), function(i) {
  dat <- Read10X(data.dir = sample_names[i])
  colnames(dat) <- gsub("1", i, colnames(dat))
  return(dat)
})

cm <- do.call(cbind, project.data)
cell_used <- read.csv("cellMetaData.csv", row.names = 1, header = TRUE)
cm <- cm[, rownames(cell_used)]

Yu_project <- CreateSeuratObject(counts = cm, meta.data = cell_used, min.cells = 60)
Yu_project <- subset(Yu_project, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & 
                       nCount_RNA < 30000 & percent.mt < 10 & percent.redcell < 10 & percent.ribo < 40)

# Filter genes
sex.genes <- c('DDX3Y','EIF2S3Y','UTY','KDM5D','XIST','TSIX')
mt.genes <- rownames(Yu_project)[grep("^MT-", rownames(Yu_project))]
red.genes <- c("HBA1","HBA2","HBB",'HBD','HBE1','HBG1','HBG2','HBM','HBQ1','HBZ') 
Yu_project <- subset(Yu_project, features = setdiff(rownames(Yu_project), c(red.genes, mt.genes, sex.genes)))

# Normalize and Cluster
Yu_project <- NormalizeData(Yu_project, scale.factor = 10000) %>% FindVariableFeatures()
Yu_project$CC.Difference <- Yu_project$S.Score - Yu_project$G2M.Score  
Yu_project <- ScaleData(Yu_project, vars.to.regress = 'CC.Difference') %>% RunPCA(npcs = 30) %>% RunHarmony("gestational.week")	

Yu_project$clusterIdent <- Yu_project$cell.type
Yu_project$clusterIdent[grep('CGE', Yu_project$cell.type)] <- "CGE"
Yu_project$clusterIdent[grep('LGE', Yu_project$cell.type)] <- "LGE"
Yu_project$clusterIdent[grep('MGE', Yu_project$cell.type)] <- "MGE"

Idents(Yu_project) <- "clusterIdent"
Yu_GE_Prog <- subset(Yu_project, idents = c("MGE","CGE","LGE","P1","P2","P3","P4","P5","P6"))
Yu_GE_Prog <- SCTransform(Yu_GE_Prog, ncells = 5000, conserve.memory = TRUE) %>% FindVariableFeatures()

# ---------------------------------------------------------
# B. Figure 1A4 (In-Vitro vs In-Vivo Correlation)
# ---------------------------------------------------------
avg_expr_Yu <- sapply(levels(Yu_GE_Prog), function(ct) {
  rowMeans(Yu_GE_Prog@assays$RNA@data[, which(Idents(Yu_GE_Prog) == ct)])
})

genes2cor <- intersect(VariableFeatures(Yu_GE_Prog), VariableFeatures(iPSC_16p11_Obj))
corr2ref_cl <- cor(avg_expr_Yu[genes2cor, ], averages_expr[genes2cor, ], method = "pearson")

heatmap.2(corr2ref_cl, scale = "none", trace = "none", key = TRUE, keysize = 0.5, 
          labRow = colnames(avg_expr_Yu), labCol = colnames(averages_expr), 
          cexRow = 0.8, cexCol = 0.8,
          col = colorRampPalette(rev(c("#b2182b","#d6604d","#f4a582","#fddbc7",
                                       "#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")))(30),
          main = "Figure 1A4: In-Vitro vs In-Vivo Correlation")


# ==============================================================================
# PART 3: Process Alternative In-Vivo Reference (Shi et al., GSE135827)
# Test only
# ==============================================================================
raw_count_shi <- read.table("GSE135827_GE_mat_raw_count_with_week_info.txt", header = TRUE, row.names = 1)
metadata_shi <- as.data.frame(readxl::read_excel("science.abj6641_table_s2.xlsx"))
metadata_shi$Week <- sapply(colnames(raw_count_shi), function(x) strsplit(x, "\\.")[[1]][2])
rownames(metadata_shi) <- colnames(raw_count_shi)

Shi <- CreateSeuratObject(counts = raw_count_shi, meta.data = metadata_shi)
Shi <- NormalizeData(Shi) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(dims = 1:3) %>% RunUMAP(dims = 1:3)

Idents(Shi) <- "Major.types"
Shi_Prog <- subset(Shi, idents = "progenitor") %>% SCTransform()

DEGs_RG_IPC <- readxl::read_excel("science.abj6641_table_s5.xlsx")
VariableFeatures(Shi_Prog) <- unique(DEGs_RG_IPC$gene)
Shi_Prog <- RunPCA(Shi_Prog, npcs = 30) %>% FindNeighbors(dims = 1:10) %>% FindClusters(resolution = 0.8)

Shi_Prog$clusterIdent <- ifelse(Shi_Prog$seurat_clusters %in% c("0","4","8","10"), "RGC", "IPC")
Shi_Prog$clusterAssign <- paste(Shi_Prog$clusterIdent, Shi_Prog$Week)

Shi_GE <- subset(Shi, idents = c("MGE","LGE","CGE"))
Shi_GE$clusterIdent <- Shi_GE$Major.types
Shi_GE$clusterAssign <- paste(Shi_GE$clusterIdent, Shi_GE$Week)

Shi_GE_Prog <- merge(Shi_Prog, y = Shi_GE, add.cell.ids = c("Prog", "GE"), project = "Shi")