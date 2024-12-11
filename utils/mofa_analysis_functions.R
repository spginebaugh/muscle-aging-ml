#' Select Important Factor from MOFA Model
#' 
#' This function identifies the most significant factor from a MOFA model based on its 
#' association with a specified test variable in the sample metadata.
#' 
#' @param mofa_model A MOFA model object
#' @param sample_metadata Data frame containing sample metadata
#' @param sample_id_column Character string specifying the column name containing sample IDs
#' @param test_variable Character string specifying the variable to test for associations
#' @param test_type Character string indicating the type of test ("categorical" or "continuous")
#' @param categorical_type Character string specifying the type of categorical test ("parametric" or "non-parametric")
#' @return Character string representing the most important factor name
select_important_factor <- function(
    mofa_model,
    sample_metadata,
    sample_id_column,
    test_variable,
    test_type = "categorical",
    categorical_type = "parametric"
) {
  association_scores <- MOFAcellulaR::get_associations(
    model = mofa_model,
    metadata = sample_metadata,
    sample_id_column = sample_id_column,
    test_variable = test_variable,
    test_type = test_type,
    categorical_type = categorical_type
  )
  
  important_factor <- association_scores$Factor[which.min(association_scores$p.value)]
  return(important_factor)
}

#' Get Associations List for Multiple Variables
#' 
#' Calculates associations between MOFA factors and multiple categorical or continuous 
#' variables from the sample metadata.
#' 
#' @param mofa_model A MOFA model object
#' @param sample_meta Data frame containing sample metadata
#' @param sample_id_column Character string specifying the column name containing sample IDs
#' @param categorical_test_var Vector of categorical variables to test
#' @param continuous_test_var Vector of continuous variables to test
#' @param categorical_type Character string specifying the type of categorical test
#' @return List of association results for each tested variable
get_associations_list <- function(
    mofa_model, 
    sample_meta, 
    sample_id_column,
    categorical_test_var = NULL,
    continuous_test_var = NULL,
    categorical_type = "parametric"
) {
  associations_list <- list()
  if (!is.null(categorical_test_var)){
    for (test_var in categorical_test_var){
      associations_list[[test_var]] <- MOFAcellulaR::get_associations(
        model = mofa_model,
        metadata = sample_meta,
        sample_id_column = sample_id_column,
        test_variable = test_var,
        test_type = "categorical",
        categorical_type = categorical_type
      )
    }
  }
  
  if (!is.null(continuous_test_var)){
    for (test_var in continuous_test_var){
      associations_list[[test_var]] <- MOFAcellulaR::get_associations(
        model = mofa_model,
        metadata = sample_meta,
        sample_id_column = sample_id_column,
        test_variable = test_var,
        test_type = "continous",
        categorical_type = categorical_type
      )
    }
  }
  
  return(associations_list)
}

#' Extract Factor Loadings Matrix
#' 
#' Retrieves and formats the gene weights (loadings) for a specified factor from 
#' a MOFA model into a matrix format.
#' 
#' @param mofa_model A MOFA model object
#' @param important_factor Character string specifying the factor name
#' @return Matrix of factor loadings with genes as rows and cell types as columns
get_factor_loadings_matrix <- function(mofa_model, important_factor){
  factor_loadings_mat <- get_geneweights(model = mofa_model, factor = important_factor) %>%
    pivot_wider(names_from = ctype, values_from = value, values_fill = 0) %>%
    column_to_rownames("feature") %>%
    as.matrix()
  return(factor_loadings_mat)
}

#' Filter Factor Loading Matrix
#' 
#' Filters a factor loading matrix to retain only genes with loadings above a 
#' specified cutoff threshold.
#' 
#' @param factor_loadings_mat Matrix of factor loadings
#' @param cutoff Numeric threshold for filtering (default: 0.3)
#' @return Filtered matrix containing only genes meeting the cutoff criterion
filter_factor_loading_matrix <- function(factor_loadings_mat, cutoff = 0.3){
  factor_loadings_mat_filt <- factor_loadings_mat[abs(rowMin(factor_loadings_mat)) >= cutoff,]
  return(factor_loadings_mat_filt)
}

#' Count Gene Occurrences in Factor
#' 
#' Calculates how many cell types each gene appears in with a loading above
#' the specified cutoff.
#' 
#' @param factor_loadings Data frame of factor loadings
#' @param cutoff Numeric threshold for counting (default: 0.3)
#' @return Data frame with gene counts sorted in descending order
get_gene_counts_in_factor <- function(factor_loadings, cutoff = 0.3){
  gene_count_in_factor <- factor_loadings %>%
    dplyr::filter(abs(value) >= 0.3) %>%
    group_by(feature) %>%
    summarize(n_cells = length(feature)) %>%
    arrange(-n_cells)
  
  return(gene_count_in_factor)
}

#' Summarize Shared Genes
#' 
#' Creates a summary table showing how many genes are shared across different 
#' numbers of cell types.
#' 
#' @param factor_loadings Data frame of factor loadings
#' @param cutoff Numeric threshold for counting (default: 0.3)
#' @return Summary table of gene sharing patterns
get_shared_genes_table <- function(factor_loadings, cutoff = 0.3){
  shared_genes_table <- get_gene_counts_in_factor(factor_loadings, cutoff) %>%
    group_by(n_cells) %>%
    summarize(count = length(n_cells))
  
  return(shared_genes_table)
}

#' Project Factor on Bulk Expression Data
#' 
#' Projects a MOFA factor onto bulk expression data using weighted means to calculate
#' factor scores for bulk samples.
#' 
#' @param norm_expr Normalized expression matrix
#' @param factor_loadings Factor loadings data frame
#' @param bulk_meta Metadata for bulk samples
#' @return Matrix of scaled factor scores for bulk samples
project_factor_on_bulk <- function(norm_expr, factor_loadings, bulk_meta) {
  fact_score <- decoupleR::run_wmean(
    norm_expr,
    network = factor_loadings,
    .source = ctype,
    .target = feature,
    .mor = value
  ) %>%
    dplyr::filter(statistic == "norm_wmean") %>%
    left_join(bulk_meta, by = c("condition" = "sample"))
  
  fact_score_mat <- fact_score %>%
    dplyr::select(condition, source, score) %>%
    pivot_wider(names_from = source, values_from = score) %>%
    column_to_rownames("condition") %>%
    as.matrix() %>%
    scale()
  
  return(fact_score_mat)
}
