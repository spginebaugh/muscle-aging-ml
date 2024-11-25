# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Load Libraries                              ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Seurat)
library(MOFA2)
library(MOFAcellulaR)
library(qs)
library(dplyr)
library(tidyverse)
library(magrittr)
library(edgeR)
library(DESeq2)
library(nichenetr)

library(ggprism)
library(ggpmisc)

set.seed(43648)

source("utils/mofa_prep_functions.R")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                Import data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data_dir <- file.path("data/mofa_inputs/")

muscle_background_genes <- read.csv(paste0(data_dir, "muscle_background_genes.csv"))$x
pb_counts <- read.csv(paste0(data_dir, "pb_counts.csv"), row.names = 1)
pb_meta <- read.csv(paste0(data_dir, "pb_meta.csv"), row.names = 1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                            filtering parameters                          ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ncells_filter <- 20
ngenes_filter <- 15
nsamples_filter <- 12
min_count_filter <- 10
min_prop_filter <- 0.25
prop_coverage_filter <- 0.9
scale_factor <- 1e6

background_genes <- muscle_background_genes


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                prep for MOFA                             ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## pseudobulk object for MOFA
pb_obj <- MOFAcellulaR::create_init_exp(count = pb_counts, coldata = pb_meta)

## filter out cell types with low cell counts (<20)
expr_list <- MOFAcellulaR::filt_profiles(
  pb_dat = pb_obj,
  cts = unique(pb_meta$cell_type), 
  ncells = ncells_filter, 
  counts_col = "cell_counts", 
  ct_col = "cell_type" 
) 

## filter out samples with low cell counts
expr_list <- MOFAcellulaR::filt_views_bysamples(
  pb_dat_list = expr_list,
  nsamples = nsamples_filter
)

## filter out genes and samples with low expression
expr_list <- MOFAcellulaR::filt_gex_byexpr(
  pb_dat_list = expr_list,
  min.count = min_count_filter,
  min.prop = min_prop_filter
)
expr_list <- MOFAcellulaR::filt_views_bygenes(pb_dat_list = expr_list, ngenes = ngenes_filter)
expr_list <- MOFAcellulaR::filt_samples_bycov(pb_dat_list = expr_list, prop_coverage = prop_coverage_filter)

## TMM normalization of pseudobulk expression
expr_list <- MOFAcellulaR::tmm_trns(pb_dat_list = expr_list, scale_factor = scale_factor)

## filter by highly variable genes
expr_list <- MOFAcellulaR::filt_gex_byhvg(
  pb_dat_list = expr_list,
  prior_hvg = NULL,
  var.threshold = 0
)

# ## remove celltype specific genes
if (!is.null(background_genes)){
  expr_list <- MOFAcellulaR::filt_gex_bybckgrnd(pb_dat_list = expr_list,
                                                prior_mrks = background_genes)
}

## final round of "view" filtering after more genes have been removed
expr_list <- MOFAcellulaR::filt_views_bygenes(pb_dat_list = expr_list, ngenes = ngenes_filter)

multiview_dat <- pb_dat2MOFA(
  pb_dat_list = expr_list,
  sample_column = "donor_id"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 save files                               ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.csv(multiview_dat, paste0(data_dir, "multiview_dat.csv"), row.names = FALSE)

