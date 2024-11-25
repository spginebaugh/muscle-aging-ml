create_factor_expr_group_boxplot <- function(all_factors, important_factor, xvar = "tech", yvar = "value", fillvar = "age_bin"){
  group_boxplot <- all_factors %>%
    dplyr::filter(Factor %in% c(important_factor)) %>%
    ggplot(aes(x = .data[[xvar]], y = .data[[yvar]], fill = .data[[fillvar]])) +
    geom_boxplot() +
    geom_point(position = position_dodge(width = .75)) +
    theme_prism() +
    scale_fill_manual(values = c("darkred", "green")) +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    ylab(paste0(important_factor, " score")) +
    xlab(fillvar)
}

create_variance_explained_barplot <- function(mofa_model, important_factor){
  variance_explained_barplot <- mofa_model@cache$variance_explained$r2_per_factor$single_group[important_factor, ] %>%
    enframe() %>%
    ggplot(aes(x = name, y = value)) +
    geom_bar(stat = "identity") +
    theme_prism() +
    theme(
      axis.text.x = element_text(angle = 90, size = 12, hjust = 1, vjust = 0.5),
      axis.text = element_text(size = 12)
    ) +
    ylab("R2") +
    xlab("cell type") 
  
  return(variance_explained_barplot)
}


create_shared_genes_plot <- function(shared_genes_table){
  shared_genes_plot <- ggplot(shared_genes_table,  aes(x = (count), y = reorder(paste(n_cells, "types"), n_cells))) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12)
    ) +
    theme(plot.margin = unit(c(0, 1, 0.2, 0.2), "cm")) +
    ylab("number of cell types") +
    xlab("number of genes")
  
  return(shared_genes_plot)
}
