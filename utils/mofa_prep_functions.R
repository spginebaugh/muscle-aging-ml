#' Create Pseudobulk Count Matrix from Seurat Object
#' 
#' Aggregates single-cell RNA sequencing data from a Seurat object into pseudobulk 
#' counts based on sample and cell type annotations.
#' 
#' @param seurat_obj Seurat object containing single-cell RNA-seq data
#' @param sample_annotation_colname Character string specifying the column name for sample annotations
#' @param cell_annotation_colname Character string specifying the column name for cell type annotations
#' @return Data frame of pseudobulk counts with genes as rows and sample-celltype combinations as columns
prep_pseudobulk_counts <- function(seurat_obj, sample_annotation_colname, cell_annotation_colname){
  pseudobulk_cts <- AggregateExpression(
    object = seurat_obj,
    assay = "RNA",
    slot = "counts",
    group.by = c(sample_annotation_colname, cell_annotation_colname)
  )
  
  pseudobulk_cts <- pseudobulk_cts[[1]] %>% data.frame()
  colnames(pseudobulk_cts) <- str_replace_all(colnames(pseudobulk_cts),fixed("."),"_")
  pseudobulk_cts <- pseudobulk_cts[, order(colnames(pseudobulk_cts))]
  
  return(pseudobulk_cts)
}

#' Create Pseudobulk Metadata from Seurat Object
#' 
#' Generates metadata for pseudobulk samples, including cell counts for each 
#' sample-celltype combination.
#' 
#' @param seurat_obj Seurat object containing single-cell RNA-seq data
#' @param sample_annotation_colname Character string specifying the column name for sample annotations
#' @param cell_annotation_colname Character string specifying the column name for cell type annotations
#' @return Data frame containing metadata for pseudobulk samples with cell counts
prep_pseudobulk_metadata <- function(seurat_obj, sample_annotation_colname, cell_annotation_colname){
  pseudobulk_meta <- table(seurat_obj@meta.data[,sample_annotation_colname], 
                           seurat_obj@meta.data[,cell_annotation_colname]) %>% 
    as.data.frame()
  colnames(pseudobulk_meta) <- c("donor_id","cell_type","cell_counts")
  pseudobulk_meta <- pseudobulk_meta[pseudobulk_meta$cell_counts > 0, ]
  rownames(pseudobulk_meta) <- paste(pseudobulk_meta$donor_id, pseudobulk_meta$cell_type, sep = "_")
  pseudobulk_meta <- pseudobulk_meta[order(rownames(pseudobulk_meta)), ]
  
  return(pseudobulk_meta)
}

#' Perform Consensus Differential Expression Analysis for Single-Cell Data
#' 
#' Identifies marker genes for a focus cell type by comparing it against all other 
#' cell types (excluding specified cell types) using Seurat's FindMarkers function.
#' 
#' @param seurat_obj Seurat object containing single-cell RNA-seq data
#' @param annotation_colname Character string specifying the column name for cell type annotations
#' @param celltype_focus Character string specifying the cell type of interest
#' @param exclude_celltype Character vector of cell types to exclude from the analysis
#' @return List of marker genes for each comparison, where each element contains genes 
#'         upregulated in the focus cell type
single_cell_de_consensus <- function(seurat_obj, annotation_colname, celltype_focus, exclude_celltype){
  celltypes <- unique(seurat_obj@meta.data[,(annotation_colname)])
  celltypes <- celltypes[!(celltypes %in% c(celltype_focus, exclude_celltype))]
  
  markers_list <- list()
  for(celltype in celltypes){
    de_out <- FindMarkers(seurat_obj, 
                          only.pos = TRUE, 
                          ident.1 = celltype_focus,
                          ident.2 = celltype,
                          min.pct = 0.1,
                          min.diff.pct = 0.1,
                          logfc.threshold = 0.5)
    
    markers_list[[celltype]] <- de_out %>% 
      rownames_to_column("gene_name") %>% 
      dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.01) %>% 
      pull(gene_name)
  }
  return(markers_list)
}

#' Identify Marker Genes from Pseudobulk Data
#' 
#' Performs differential expression analysis on pseudobulk data to identify marker 
#' genes for a focus cell type compared to all other cell types.
#' 
#' @param pseudobulk_counts Matrix of pseudobulk count data
#' @param pseudobulk_metadata Data frame containing metadata for pseudobulk samples
#' @param celltype_focus Character string specifying the cell type of interest
#' @param exclude_celltype Character vector of cell types to exclude from the analysis
#' @return List of marker genes for each comparison, where each element contains genes 
#'         upregulated in the focus cell type
pseudobulk_marker_genes <- function(pseudobulk_counts, pseudobulk_metadata, celltype_focus, exclude_celltype){
  pseudobulk_metadata$cell_type <- as.character(pseudobulk_metadata$cell_type)
  celltypes <- unique(pseudobulk_metadata$cell_type)
  celltypes <- celltypes[!(celltypes %in% c(celltype_focus, exclude_celltype))]
  
  markers_list <- list()
  for(celltype in celltypes){
    print(celltype)
    res <- pseudobulk_de(pseudobulk_counts, pseudobulk_metadata, celltype, celltype_focus)
    
    markers_list[[celltype]] <- topTags(res, n = Inf) %>%
      as.data.frame() %>%
      rownames_to_column("gene_name") %>% 
      dplyr::filter(logFC > 0.5 & FDR < 0.01) %>% 
      pull(gene_name)
    
  }
  return(markers_list)
}

#' Perform Differential Expression Analysis on Pseudobulk Data
#' 
#' Internal function that performs differential expression analysis between two cell 
#' types using edgeR on pseudobulk data.
#' 
#' @param pseudobulk_counts Matrix of pseudobulk count data
#' @param pseudobulk_metadata Data frame containing metadata for pseudobulk samples
#' @param celltype Character string specifying the reference cell type
#' @param celltype_focus Character string specifying the cell type of interest
#' @return edgeR results object containing differential expression statistics
#' @details Uses edgeR's quasi-likelihood F-tests with robust estimation. Includes 
#'          filtering for expressed genes and TMM normalization.
pseudobulk_de <- function(pseudobulk_counts, pseudobulk_metadata, celltype, celltype_focus){
  input_counts <- pseudobulk_counts[, pseudobulk_metadata$cell_type %in% c(celltype, celltype_focus)]
  input_metadata <- pseudobulk_metadata[pseudobulk_metadata$cell_type %in% c(celltype, celltype_focus),]
  
  dat <- DGEList(input_counts, samples = DataFrame(input_metadata))
  keep <- filterByExpr(dat, group = input_metadata$cell_type)
  
  dat <- dat[keep,]
  dat <- calcNormFactors(dat)
  design <- model.matrix(~factor(cell_type,
                                 levels = c(celltype,celltype_focus)), dat$samples)
  
  colnames(design) <- c("int", celltype_focus)
  
  dat <- estimateDisp(dat, design)
  
  fit <- glmQLFit(dat, design, robust=TRUE)
  
  res <- glmQLFTest(fit, coef=ncol(design))
  
  return(res)
}