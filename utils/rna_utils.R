#' Import and Create Merged Seurat Object
#' 
#' Imports single-cell RNA sequencing data from multiple samples and creates a merged 
#' Seurat object. Supports both H5 and MTX file formats from Cell Ranger output.
#' 
#' @param cellranger_folder_path Character string specifying the base path to Cell Ranger output
#' @param file_names_vec Vector of sample names/folders to import
#' @param file_h5_path Character string specifying the path to H5 file relative to sample folder (default: NA)
#' @param import_method Character string specifying import format ("h5" or "mtx", default: "h5")
#' @param subset_cell_names Vector of cell names to subset (default: NA)
#' @return A merged Seurat object containing all samples
#' @details 
#' For H5 files:
#'   - Expects one H5 file per sample in the specified directory
#' For MTX files:
#'   - Expects matrix, barcodes, and features files in each sample directory
#' 
#' The function:
#' 1. Reads each sample's data
#' 2. Creates individual Seurat objects
#' 3. Merges all objects
#' 4. Adds sample information to metadata
#' 5. Adds cell barcodes to metadata
import_seurat <- function(cellranger_folder_path, 
                          file_names_vec, 
                          file_h5_path = NA, 
                          import_method = "h5",
                          subset_cell_names = NA) {
  
  seurat_list <- list()
  
  for (file_index in 1:length(file_names_vec)) {
    if (import_method == "h5") {
      seurat_data <- Read10X_h5(filename = paste0(
        cellranger_folder_path,
        file_names_vec[file_index],
        file_h5_path
      ))
    } else if (import_method == "mtx") {
      seurat_data <- Read10X(data.dir = paste0(
        cellranger_folder_path,
        file_names_vec[file_index]
      ))
    } else {
      stop("import_method must be either 'h5' or 'mtx' ")
    }
    
    if (all(!is.na(subset_cell_names))){
      celltype_names <- paste0(file_names_vec[file_index],"_",colnames(seurat_data))
      seurat_data <- seurat_data[,celltype_names %in% subset_cell_names]
    }
    
    seurat_list[file_index] <- CreateSeuratObject(
      counts = seurat_data,
      project = file_names_vec[file_index]
    )
  }
  
  names(seurat_list) <- file_names_vec
  
  seurat <- merge(
    x = seurat_list[[1]],
    y = seurat_list[2:length(seurat_list)],
    add.cell.id = file_names_vec
  )
  
  seurat <- JoinLayers(seurat) # merge samples in Seurat V5
  
  colnames(seurat@meta.data)[1] <- "sample"
  seurat$barcode <- colnames(seurat)
  
  return(seurat)
}

#' Import CellBender Output to Merged Seurat Object
#' 
#' Imports and merges single-cell RNA sequencing data that has been processed 
#' with CellBender for ambient RNA removal.
#' 
#' @param cellranger_folder_path Character string specifying the base path to CellBender output
#' @param file_names_vec Vector of sample names/folders to import
#' @param file_h5_path Character string specifying the path to H5 file relative to sample folder (default: NA)
#' @return A merged Seurat object containing all CellBender-processed samples
#' @details 
#' Similar to import_seurat but specifically handles CellBender output format.
#' The function:
#' 1. Reads CellBender-processed H5 files
#' 2. Creates individual Seurat objects
#' 3. Merges all objects
#' 4. Adds sample information and cell barcodes to metadata
import_seurat_from_cellbender <- function(cellranger_folder_path, file_names_vec, file_h5_path = NA) {
  
  seurat_list <- list()
  
  for (file_index in 1:length(file_names_vec)) {
    seurat_data <- Read_CellBender_h5_Mat(file_name = paste0(
      cellranger_folder_path,
      file_names_vec[file_index],
      file_h5_path
    ))
    seurat_list[file_index] <- CreateSeuratObject(
      counts = seurat_data,
      project = file_names_vec[file_index]
    )
  }
  
  names(seurat_list) <- file_names_vec
  
  seurat <- merge(
    x = seurat_list[[1]],
    y = seurat_list[2:length(seurat_list)],
    add.cell.id = file_names_vec
  )
  
  seurat <- JoinLayers(seurat) # merge samples in Seurat V5
  
  colnames(seurat@meta.data)[1] <- "sample"
  seurat$barcode <- colnames(seurat)
  
  return(seurat)
}

#' Score Clusters Using PanglaoDB Markers
#' 
#' Calculates cell type scores for each cluster using marker genes from PanglaoDB database.
#' 
#' @param seurat_obj Seurat object with defined clusters
#' @return Matrix of scaled cell type scores for each cluster
#' @details 
#' The function:
#' 1. Loads PanglaoDB marker genes (human and mouse)
#' 2. Filters for genes present in the dataset
#' 3. Calculates module scores for each cell type
#' 4. Aggregates scores by cluster
#' 5. Returns a matrix where:
#'    - Rows are cell types from PanglaoDB
#'    - Columns are cluster numbers
#'    - Values are scaled average module scores
#' 
#' Note: Requires PanglaoDB marker file at "/media/largedata/universal_files/PanglaoDB_markers_27_Mar_2020.tsv"
score_panglao <- function(seurat_obj){
  # returns a matrix of scaled celltype scores for each cluster 
  panglao <- readr::read_tsv("/media/largedata/universal_files/PanglaoDB_markers_27_Mar_2020.tsv")
  panglao <- panglao[panglao$species %in% c("Hs", "Mm Hs"), ]
  panglao <- panglao[panglao$`official gene symbol` %in% rownames(seurat_obj), ]
  panglao <- split(panglao$`official gene symbol`, panglao$`cell type`)
  
  seurat_obj <- AddModuleScore(seurat_obj, panglao)
  
  metadata <- seurat_obj@meta.data
  
  metadata_name_start <- (ncol(metadata) - length(panglao) + 1)
  metadata_name_end <- ncol(metadata)
  
  pang_ann <- metadata[,metadata_name_start:metadata_name_end] 
  colnames(pang_ann)<- names(panglao)
  
  pang_ann <- scale(pang_ann)
  agg_pang <- aggregate(pang_ann, list(seurat_obj@active.ident), mean)
  rownames(agg_pang) <- agg_pang[, 1]
  agg_pang <- agg_pang[, -1]
  agg_pang <- t(agg_pang)
  return(agg_pang)
}