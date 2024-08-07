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

set.seed(43648)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Functions                                ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
importSeurat <- function(cellranger_folder_path, file_names_vec, file_h5_path) {
  seurat_list <- list()
  for (file_index in 1:length(file_names_vec)) {
    seurat_data <- Read10X_h5(filename = paste0(
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

  return(seurat)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                          Import Cellranger Outputs                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## note: want to use cellranger output with nuclei selection of cellbender

file_names_vec <- list.files("data/raw/cellranger_outputs/")


## loop over cellbender outputs for each sample and place in seurat object
seurat_bender <- importSeurat(
  cellranger_folder_path = "data/raw/cellranger_outputs/",
  file_names_vec = file_names_vec,
  file_h5_path = "/outs/cellbender/cellbender_out_filtered.h5"
)

## import uncorrected data from cellranger as alternative and add to seurat object
## this will allow us to use cellbender's empty droplet selection with cellranger outputs
seurat_ranger <- importSeurat(
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                scDblFinder                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (sample_index in 1:length(seurat_split)) {
  dbl_out <- scDblFinder(
    GetAssayData(seurat_split[[sample_index]], assay = "RNA", slot = "counts"),
    returnType = "table"
  ) %>% as.data.frame()
  dbl_out$barcodes <- rownames(dbl_out)
  dbl_out <- dbl_out[, c("barcodes", "class")]
  colnames(dbl_out) <- c("barcodes", "scDbl_class")

  ## add scoreing into seurat metadata
  metadata <- seurat_split[[sample_index]]@meta.data
  metadata <- left_join(metadata, dbl_out, by = "barcodes")
  rownames(metadata) <- metadata$barcodes
  seurat_split[[sample_index]]@meta.data <- metadata
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                               scds detection                             ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (sample_index in 1:length(seurat_split)) {
  sce <- as.SingleCellExperiment(seurat_split[[sample_index]])

  sce <- cxds(sce, retRes = TRUE, verb = TRUE)
  sce <- bcds(sce, retRes = TRUE, verb = TRUE)
  sce <- cxds_bcds_hybrid(sce, verb = TRUE)

  seurat_split[[sample_index]]$cxds_scores <- sce$cxds_score
  seurat_split[[sample_index]]$bcds_scores <- sce$bcds_score
  seurat_split[[sample_index]]$hybrid_scores <- sce$hybrid_score
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                DoubletFinder                             ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Extract Metadata                            ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
metalist <- lapply(seurat_split, function(x) {
  metadata <- x@meta.data
  return(metadata)
})
metalist <- lapply(metalist, function(x) {
  colnames(x)[(ncol(x) - 1):ncol(x)] <- c("DF_score", "DF_classification")
  return(x)
})

metalist <- do.call("rbind", metalist[1:length(metalist)])

seurat_ranger$scDbl_class <- metalist$scDbl_class
seurat_ranger$cxds_scores <- metalist$cxds_score
seurat_ranger$bcds_scores <- metalist$bcds_score
seurat_ranger$hybrid_scores <- metalist$hybrid_score
seurat_ranger$DF_score <- metalist$DF_score
seurat_ranger$DF_classification <- metalist$DF_classification


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                Save Metadata                             ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
metadata <- seurat_ranger@meta.data
write.csv(metadata, "data/processed/perez_doublet_metadata.csv")
