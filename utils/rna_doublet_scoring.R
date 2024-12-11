#' Run Doublet Detection with scDblFinder
#' 
#' Applies the scDblFinder algorithm to detect doublets in a list of Seurat objects.
#' 
#' @param seurat_split List of Seurat objects split by sample
#' @return List of Seurat objects with scDblFinder results added to metadata
#' @details Processes each sample separately and adds two columns to the metadata:
#'          'barcode' and 'scDbl_class'. The classification is stored in the 
#'          'scDbl_class' column of the metadata.
run_scDblFinder <- function(seurat_split){
  for (sample_index in 1:length(seurat_split)) {
    dbl_out <- scDblFinder(
      GetAssayData(seurat_split[[sample_index]], assay = "RNA", slot = "counts"),
      returnType = "table"
    ) %>% as.data.frame()
    dbl_out$barcode <- rownames(dbl_out)
    dbl_out <- dbl_out[, c("barcode", "class")]
    colnames(dbl_out) <- c("barcode", "scDbl_class")
    
    ## add scoring into seurat metadata
    metadata <- seurat_split[[sample_index]]@meta.data
    metadata <- left_join(metadata, dbl_out, by = "barcode")
    rownames(metadata) <- metadata$barcode
    seurat_split[[sample_index]]@meta.data <- metadata
  }
  return(seurat_split)
}

#' Run SCDS Doublet Detection
#' 
#' Applies the SCDS algorithm suite (cxds, bcds, and hybrid) to detect doublets 
#' in a list of Seurat objects.
#' 
#' @param seurat_split List of Seurat objects split by sample
#' @return List of Seurat objects with SCDS scores added to metadata
#' @details Adds three scoring columns to the metadata:
#'          - cxds_score: Co-expression based doublet scoring
#'          - bcds_score: Binary classification-based doublet scoring
#'          - hybrid_score: Hybrid scoring combining both methods
run_scds <- function(seurat_split){
  for (sample_index in 1:length(seurat_split)) {
    sce <- as.SingleCellExperiment(seurat_split[[sample_index]])
    
    sce <- cxds(sce, retRes = TRUE, verb = TRUE)
    sce <- bcds(sce, retRes = TRUE, verb = TRUE)
    sce <- cxds_bcds_hybrid(sce, verb = TRUE)
    
    seurat_split[[sample_index]]$cxds_score <- sce$cxds_score
    seurat_split[[sample_index]]$bcds_score <- sce$bcds_score
    seurat_split[[sample_index]]$hybrid_score <- sce$hybrid_score
  }
  return(seurat_split)
}

#' Run DoubletFinder Detection
#' 
#' Applies the DoubletFinder algorithm to detect doublets in a list of Seurat objects.
#' 
#' @param seurat_split List of Seurat objects split by sample
#' @return List of Seurat objects with DoubletFinder results added to metadata
#' @details Performs parameter optimization using paramSweep and assumes a 5% doublet 
#'          formation rate. Uses the first 10 principal components for analysis.
#'          Adds classification and scoring columns to the metadata.
run_doubletfinder <- function(seurat_split){
  ## pK Identification (no ground-truth)
  gc()
  sweep_res <- list()
  sweep_stats <- list()
  pk_out <- list()
  for (sample_index in 1:length(seurat_split)) {
    sweep_res[[sample_index]] <- paramSweep(seurat_split[[sample_index]], PCs = 1:10, sct = FALSE)
    sweep_stats[[sample_index]] <- summarizeSweep(sweep_res[[sample_index]], GT = FALSE)
    pk_out[[sample_index]] <- find.pK(sweep_stats[[sample_index]])
  }
  
  ## Homotypic Doublet Proportion Estimate
  nexp_poi <- list()
  for (sample_index in 1:length(seurat_split)) {
    nexp_poi[[sample_index]] <- round(0.05 * ncol(seurat_split[[sample_index]])) ## Assuming 5.0% doublet formation rate
  }
  
  ## Run DoubletFinder
  for (sample_index in 1:length(seurat_split)) {
    seurat_split[[sample_index]] <- doubletFinder(
      seurat_split[[sample_index]],
      PCs = 1:10, pN = 0.25,
      pK = pk_out[[sample_index]]$pK[which.max(pk_out[[sample_index]]$BCmetric)] %>% as.character() %>% as.numeric(),
      nExp = nexp_poi[[sample_index]], reuse.pANN = FALSE, sct = FALSE
    )
  }
  
  return(seurat_split)
}

#' Add Doublet Detection Results to Merged Seurat Object
#' 
#' Combines doublet detection results from multiple algorithms into a single 
#' merged Seurat object.
#' 
#' @param seurat_obj Original merged Seurat object
#' @param seurat_split List of Seurat objects with doublet detection results
#' @return Merged Seurat object with all doublet detection results in metadata
#' @details Adds the following columns to the metadata:
#'          - scDbl_class: Classification from scDblFinder
#'          - cxds_score: Score from SCDS cxds
#'          - bcds_score: Score from SCDS bcds
#'          - hybrid_score: Score from SCDS hybrid method
#'          - doubletfind_score: Score from DoubletFinder
#'          - doubletfind_class: Classification from DoubletFinder
add_doublet_metadata <- function(seurat_obj, seurat_split){
  metalist <- lapply(seurat_split, function(x) {
    metadata <- x@meta.data
    return(metadata)
  })
  metalist <- lapply(metalist, function(x) {
    colnames(x)[(ncol(x) - 1):ncol(x)] <- c("doubletfind_score", "doubletfind_class")
    return(x)
  })
  
  metalist <- do.call("rbind", metalist[1:length(metalist)])
  
  seurat_obj$scDbl_class <- metalist$scDbl_class
  seurat_obj$cxds_score <- metalist$cxds_score
  seurat_obj$bcds_score <- metalist$bcds_score
  seurat_obj$hybrid_score <- metalist$hybrid_score
  seurat_obj$doubletfind_score <- metalist$doubletfind_score
  seurat_obj$doubletfind_class <- metalist$doubletfind_class %>% tolower()
  
  return(seurat_obj)
}

#' Run Complete Doublet Detection Pipeline
#' 
#' Executes a comprehensive doublet detection workflow using multiple algorithms.
#' 
#' @param seurat_obj Seurat object to analyze
#' @param split.by Column name to split the dataset by (default: "sample")
#' @return Seurat object with results from all doublet detection methods added to metadata
#' @details Pipeline steps:
#'          1. Splits Seurat object by sample
#'          2. Performs standard preprocessing (normalization, feature selection, scaling, PCA)
#'          3. Runs scDblFinder
#'          4. Runs SCDS (cxds, bcds, hybrid)
#'          5. Runs DoubletFinder
#'          6. Combines all results in the original Seurat object
run_all_doublet_detection <- function(seurat_obj, split.by = "sample"){
  seurat_split <- SplitObject(seurat_obj, split.by = "sample")
  
  for (i in 1:length(seurat_split)) {
    seurat_split[[i]] <- NormalizeData(seurat_split[[i]])
    seurat_split[[i]] <- FindVariableFeatures(seurat_split[[i]])
    seurat_split[[i]] <- ScaleData(seurat_split[[i]])
    seurat_split[[i]] <- RunPCA(seurat_split[[i]])
  }
  ## functions pulled in from 'rna_doublet_scoring_functions.R'
  seurat_split <- run_scDblFinder(seurat_split)
  seurat_split <- run_scds(seurat_split)
  seurat_split <- run_doubletfinder(seurat_split)
  seurat_obj <- add_doublet_metadata(seurat_obj, seurat_split)
  return(seurat_obj)
}