## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                    Run MOFA on Perez snRNAseq Data                        ----
##                                                                            ~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              Load Libraries                              ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Seurat)
library(MOFA2)
library(MOFAcellulaR)
library(qs)
library(dplyr)
library(tidyverse)
library(magrittr)
library(edgeR)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                 Functions                                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                        Import and prep snRNAseq Data                     ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(43648)
## initialize filtering cutoff variables
ncells_filter <- 20
ngenes_filter <- 15
nsamples_filter <- 2
min_count_filter <- 10
min_prop_filter <- 0.25
prop_coverage_filter <- 0.9


## import seurat object
seurat <- qread("data/processed/annotated_seurat_perez.qs")

## create pseudobulk counts of each celltype for each sample (called "views")
sudo_cts <- AggregateExpression(
  object = seurat,
  assay = "RNA_ranger",
  slot = "counts",
  group.by = c("sample", "annotation_level1")
)
sudo_cts <- sudo_cts[[1]] %>% data.frame()
sudo_cts <- sudo_cts[, order(colnames(sudo_cts))]


## prep pseudobulk metadata
sudo_meta <- table(seurat$sample, seurat$annotation_level1) %>% as.data.frame()
colnames(sudo_meta) <- c("donor_id", "cell_type", "cell_counts")
sudo_meta <- sudo_meta[sudo_meta$cell_counts > 0, ]
rownames(sudo_meta) <- paste(sudo_meta$donor_id, sudo_meta$cell_type, sep = "_")
sudo_meta <- sudo_meta[order(rownames(sudo_meta)), ]

table(colnames(sudo_cts) == rownames(sudo_meta)) ## should be all true

## pseudobulk object for MOFA
sudo_obj <- MOFAcellulaR::create_init_exp(count = sudo_cts, coldata = sudo_meta)


## filter out cell types with low cell counts (<20)
expr_list <- MOFAcellulaR::filt_profiles(
  pb_dat = sudo_obj,
  cts = unique(seurat$annotation_level1), # all celltypes
  ncells = ncells_filter, # celltype cutoff of <20
  counts_col = "cell_counts", # The column name in testcoldata where the number of cells per profile was stored
  ct_col = "cell_type" # The column name in testcoldata where the cell-type label was stored
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
expr_list <- MOFAcellulaR::tmm_trns(pb_dat_list = expr_list, scale_factor = 1000000)

## filter by highly variable genes
expr_list <- MOFAcellulaR::filt_gex_byhvg(
  pb_dat_list = expr_list,
  prior_hvg = NULL,
  var.threshold = 0
)

# ## remove celltype specific genes
# ct_list <- MOFAcellulaR::filt_gex_bybckgrnd(pb_dat_list = ct_list,
#                                             prior_mrks = celltype_specific_genes)

## final round of "view" filtering after more genes have been removed
expr_list <- MOFAcellulaR::filt_views_bygenes(pb_dat_list = expr_list, ngenes = ngenes_filter)

multiview_dat <- pb_dat2MOFA(
  pb_dat_list = expr_list,
  sample_column = "donor_id"
)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                               Fit MOFA model                             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mofa_obj <- MOFA2::create_mofa(multiview_dat)

data_opts <- MOFA2::get_default_data_options(mofa_obj)
train_opts <- MOFA2::get_default_training_options(mofa_obj)
model_opts <- MOFA2::get_default_model_options(mofa_obj)

## FALSE avoids less sparse gene weights
model_opts$spikeslab_weights <- TRUE

## Define the number of factors needed (based on trial and error, unfortunately)
model_opts$num_factors <- 4

## create MOFA model
mofa_obj <- MOFA2::prepare_mofa(
  object = mofa_obj,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

## specify output path
out_file <- file.path("./vignettemodel.hdf5")

## run MOFA2
mofa_model <- MOFA2::run_mofa(mofa_obj, out_file, use_basilisk = TRUE)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                            Explore MOFA results                          ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## patient info metadata
sample_meta <- unique(seurat@meta.data[, c("sample", "age_group", "age_numeric")])

## get all factor expression in each sample
all_factors <- MOFAcellulaR::get_tidy_factors(
  model = mofa_model,
  metadata = sample_meta,
  factor = "all",
  sample_id_column = "sample"
)


categorical_assoc <- MOFAcellulaR::get_associations(
  model = mofa_model,
  metadata = sample_meta,
  sample_id_column = "sample",
  test_variable = "age_group",
  test_type = "categorical",
  group = FALSE
)
categorical_assoc


continuous_assoc <- MOFAcellulaR::get_associations(
  model = mofa_model,
  metadata = sample_meta,
  sample_id_column = "sample",
  test_variable = "age_numeric",
  test_type = "continuous",
  group = FALSE
)
continuous_assoc


mofa_model@cache$variance_explained$r2_total

mofa_model@cache$variance_explained$r2_per_factor$single_group[, , drop = F]



assoc_list <- list("categorical" = categorical_assoc, "continuous" = continuous_assoc)

plot_MOFA_hmap(
  model = mofa_model,
  group = FALSE,
  metadata = sample_meta,
  sample_id_column = "sample",
  sample_anns = c("age_group", "age_numeric"),
  assoc_list = assoc_list
)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              save MOFA model                             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qsave(mofa_model, "data/processed/mofa_model_perez.qs")













