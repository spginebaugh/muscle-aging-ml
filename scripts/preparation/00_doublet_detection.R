# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                                                            ~~
#          Doublet Detection for Perez snRNAseq Data                      ----
#                                                                            ~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Load Libraries                              ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Seurat)
library(scCustomize)
library(qs)
library(ggplot2)
library(magrittr)
library(dplyr)
library(SingleCellExperiment)

library(scds)
library(scater)
library(DoubletFinder)
library(scDblFinder)


source("scripts/utils/rna_doublet_scoring.R")
source("scripts/utils/rna_utils.R")
set.seed(43648)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                          Import Cellranger Outputs                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## note: want to use cellranger output with nuclei selection of cellbender

file_names_vec <- list.files("data/raw/cellranger_outputs/")


## loop over cellbender outputs for each sample and place in seurat object
seurat_bender_01 <- import_seurat_from_cellbender(
  cellranger_folder_path = "data/raw/cellranger_outputs/",
  file_names_vec = file_names_vec,
  file_h5_path = "/outs/cellbender_3k_20k/cellbender_out_FPR_0.01_filtered.h5"
)

seurat_bender_05 <- import_seurat_from_cellbender(
  cellranger_folder_path = "data/raw/cellranger_outputs/",
  file_names_vec = file_names_vec,
  file_h5_path = "/outs/cellbender_3k_20k/cellbender_out_FPR_0.05_filtered.h5"
)

## import uncorrected data from cellranger as alternative and add to seurat object
## this will allow us to use cellbender's empty droplet selection with cellranger outputs
seurat_ranger <- import_seurat(
  cellranger_folder_path = "data/raw/cellranger_outputs/",
  file_names_vec = file_names_vec,
  file_h5_path = "/outs/raw_feature_bc_matrix.h5",
  subset_cell_names = colnames(seurat_bender_01)
)

## remake seurat_ranger so that has the same cell ordering as seurat_bender
sr_counts <- GetAssayData(seurat_ranger, layer = "counts")
sr_counts <- sr_counts[,order(match(colnames(sr_counts), colnames(seurat_bender_01)))]
seurat_ranger <- CreateSeuratObject(counts = sr_counts)
colnames(seurat_ranger@meta.data)[1] <- "sample"
seurat_ranger$barcode <- colnames(seurat_ranger)
rm(sr_counts)
gc()

## check
table(colnames(seurat_bender_01) == colnames(seurat_bender_05))
table(colnames(seurat_bender_01) == colnames(seurat_ranger))
table(rownames(seurat_bender_01) == rownames(seurat_bender_05))
table(rownames(seurat_bender_01) == rownames(seurat_ranger))

get_cellranger_barcodes <- function(cellranger_folder_path, file_names_vec, file_barcodes_path){
  cell_names_vec <- character()
  for (file_index in 1:length(file_names_vec)) {
    cell_names <- read.csv(paste0(
      cellranger_folder_path,
      file_names_vec[file_index],
      file_barcodes_path), 
      header = FALSE)$V1
    
    cell_names <- paste0(file_names_vec[file_index], "_", cell_names)
    cell_names_vec <- c(cell_names_vec, cell_names)
  }
  return(cell_names_vec)
}

seurat_ranger_filtered_cells <- get_cellranger_barcodes(
  cellranger_folder_path = "data/raw/cellranger_outputs/",
  file_names_vec = file_names_vec,
  file_barcodes_path = "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Detect Doublets                             ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat_bender_01 <- run_all_doublet_detection(seurat_bender_01)
seurat_bender_05 <- run_all_doublet_detection(seurat_bender_05)
seurat_ranger <- run_all_doublet_detection(seurat_ranger)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                Save Metadata                             ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
metadata <- seurat_ranger@meta.data
write.csv(metadata, "data/processed/perez_doublet_metadata.csv")
