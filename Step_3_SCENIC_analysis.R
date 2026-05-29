# ==============================================================================
# STEP 1: Data Preparation
# ==============================================================================
suppressMessages({
  library(Seurat)
  library(Matrix)
})

# 1. Load the Master Seurat Object
cat("Loading Seurat Object...\n")
iPSC_16p11_Obj <- readRDS("./filtered_object_19th_Jan.rds") # Use your actual path
Idents(iPSC_16p11_Obj) <- "broad_Type"

# 2. Define the target cell type (e.g., Progenitors)
target_celltype <- "Progenitor"
Prog_Obj <- subset(iPSC_16p11_Obj, idents = target_celltype)

# 3. Split by Condition to run SCENIC separately (CON vs DEL)
# This allows us to compare how the regulatory network breaks/changes in the DEL group
conditions <- unique(Prog_Obj$Conditions) # Should be "CON" and "DEL"

for (cond in conditions) {
  cat(paste0("Processing ", target_celltype, " - ", cond, "...\n"))
  
  # Subset by condition
  subset_obj <- subset(Prog_Obj, subset = Conditions == cond)
  
  # Extract raw counts (SCENIC requires raw integer counts, not normalized data)
  counts_matrix <- subset_obj@assays$RNA@counts
  
  # Filter genes (keep genes expressed in at least 5% of cells in this subset)
  min_cells <- round(0.05 * ncol(counts_matrix))
  counts_matrix <- counts_matrix[rowSums(counts_matrix > 0) >= min_cells, ]
  
  # Convert to dense matrix and transpose (PySCENIC expects Cells as Rows, Genes as Columns)
  dense_mat <- t(as.matrix(counts_matrix))
  
  # Export to CSV (Standard PySCENIC input)
  output_file <- paste0("SCENIC_Input_", target_celltype, "_", cond, ".csv")
  write.csv(dense_mat, file = output_file, row.names = TRUE)
  
  cat(paste0("Exported: ", output_file, "\n"))
}
cat("All SCENIC inputs generated successfully!\n")

#!/bin/bash
# ==============================================================================
# STEP 2: The Core PySCENIC Pipeline
# ==============================================================================

# --- Define Global Variables & Database Paths ---
CELLTYPE="Progenitor"
CONDITIONS=("CON" "DEL")
THREADS=16 # Adjust based on your server capacity

# PySCENIC Reference Databases (Update these paths to your downloaded cisTarget files)
TFS_LIST="hs_hgnc_tfs.txt"
FEATHER_DB="hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
MOTIF_ANNOTATION="motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
MANUAL_REGULONS="Manual_Regulons_Chu.tsv"

echo "Starting Unified PySCENIC Pipeline..."

# Loop through each condition (CON and DEL)
for COND in "${CONDITIONS[@]}"; do
    echo "========================================================="
    echo " processing condition: ${COND}"
    echo "========================================================="
    
    INPUT_MAT="SCENIC_Input_${CELLTYPE}_${COND}.csv"
    ADJ_FILE="adjacencies_${CELLTYPE}_${COND}.tsv"
    REGULON_FILE="regulons_${CELLTYPE}_${COND}.csv"
    OUTPUT_LOOM="pyscenic_output_${CELLTYPE}_${COND}.loom"
    
    # --------------------------------------------------------------------------
    # Step 2.1: GRN Inference (Arboreto / GRNBoost2)
    # Infers co-expression modules between TFs and target genes
    # --------------------------------------------------------------------------
    echo "Running GRNBoost2 for ${COND}..."
    pyscenic grn \
        --num_workers ${THREADS} \
        -o ${ADJ_FILE} \
        --method grnboost2 \
        ${INPUT_MAT} \
        ${TFS_LIST}

    # --------------------------------------------------------------------------
    # Step 2.2: Regulon Prediction (cisTarget)
    # Prunes co-expression modules using motif enrichment (Direct Targets only)
    # --------------------------------------------------------------------------
    echo "Running cisTarget for ${COND}..."
    pyscenic ctx \
        ${ADJ_FILE} \
        ${FEATHER_DB} \
        --annotations_fname ${MOTIF_ANNOTATION} \
        --expression_mtx_fname ${INPUT_MAT} \
        --mode "dask_multiprocessing" \
        --num_workers ${THREADS} \
        --output ${REGULON_FILE}

    # --------------------------------------------------------------------------
    # Step 2.3: Cellular Enrichment (AUCell) + Manual Regulons
    # Scores the activity of each regulon in every single cell.
    # We append the Manual_Regulons_Chu.tsv to be scored simultaneously!
    # --------------------------------------------------------------------------
    echo "Running AUCell (including Manual Chu Regulons) for ${COND}..."
    
    # Optional: If your manual TSV needs to be converted to GMT/CSV for pySCENIC, 
    # it is usually passed alongside the computed regulons.
    pyscenic aucell \
        ${INPUT_MAT} \
        ${REGULON_FILE} \
        --signatures ${MANUAL_REGULONS} \
        --output ${OUTPUT_LOOM} \
        --num_workers ${THREADS}
        
    echo "Finished processing ${COND}. Output saved to ${OUTPUT_LOOM}."
done

echo "All SCENIC Pipeline Complete!"