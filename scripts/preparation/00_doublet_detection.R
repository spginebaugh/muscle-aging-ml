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
library(qs)
library(ggplot2)
library(magrittr)
library(dplyr)
library(SingleCellExperiment)

library(scds)
library(scater)
library(DoubletFinder)
library(scDblFinder)


source("utils/rna_doublet_scoring.R")
source("utils/rna_utils.R")
set.seed(43648)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                          Import Cellranger Outputs                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## note: want to use cellranger output with nuclei selection of cellbender

file_names_vec <- list.files("data/raw/cellranger_outputs/")


## loop over cellbender outputs for each sample and place in seurat object
seurat_bender <- import_seurat(
  cellranger_folder_path = "data/raw/cellranger_outputs/",
  file_names_vec = file_names_vec,
  file_h5_path = "/outs/cellbender/cellbender_out_filtered.h5"
)

## import uncorrected data from cellranger as alternative and add to seurat object
## this will allow us to use cellbender's empty droplet selection with cellranger outputs
seurat_ranger <- import_seurat(
  cellranger_folder_path = "data/raw/cellranger_outputs/",
  file_names_vec = file_names_vec,
  file_h5_path = "/outs/raw_feature_bc_matrix.h5"
)


## cellbender has better emptydrop removal than cellranger
## subset to only cells, as determined by cellbender
seurat_ranger <- seurat_ranger[, colnames(seurat_ranger) %in% colnames(seurat_bender)]

seurat_ranger$barcodes <- colnames(seurat_ranger) # add barcodes for merging metadata later

## split and process
seurat_split <- SplitObject(seurat_ranger, split.by = "sample")

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
seurat_ranger <- add_doublet_metadata(seurat_ranger, seurat_split)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                Save Metadata                             ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
metadata <- seurat_ranger@meta.data
write.csv(metadata, "data/processed/perez_doublet_metadata.csv")
