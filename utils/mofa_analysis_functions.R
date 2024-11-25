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


get_factor_loadings_matrix <- function(mofa_model, important_factor){
  factor_loadings_mat <- get_geneweights(model = mofa_model, factor = important_factor) %>%
    pivot_wider(names_from = ctype, values_from = value, values_fill = 0) %>%
    column_to_rownames("feature") %>%
    as.matrix()
  return(factor_loadings_mat)
}

filter_factor_loading_matrix <- function(factor_loadings_mat, cutoff = 0.3){
  factor_loadings_mat_filt <- factor_loadings_mat[abs(rowMin(factor_loadings_mat)) >= cutoff,]
  return(factor_loadings_mat_filt)
}

get_gene_counts_in_factor <- function(factor_loadings, cutoff = 0.3){
  gene_count_in_factor <- factor_loadings %>%
    dplyr::filter(abs(value) >= 0.3) %>%
    group_by(feature) %>%
    summarize(n_cells = length(feature)) %>%
    arrange(-n_cells)
  
  return(gene_count_in_factor)
}

get_shared_genes_table <- function(factor_loadings, cutoff = 0.3){
  shared_genes_table <- get_gene_counts_in_factor(factor_loadings, cutoff) %>%
    group_by(n_cells) %>%
    summarize(count = length(n_cells))
  
  return(shared_genes_table)
}

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

