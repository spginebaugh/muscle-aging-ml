#' Create Boxplot of Factor Expression by Groups
#' 
#' Generates a boxplot showing the distribution of factor scores across different 
#' groups, with optional stratification by a fill variable.
#' 
#' @param all_factors Data frame containing factor scores and metadata
#' @param important_factor Character string specifying the factor to plot
#' @param xvar Character string specifying the variable to use for x-axis grouping (default: "tech")
#' @param yvar Character string specifying the variable containing factor scores (default: "value")
#' @param fillvar Character string specifying the variable for fill color grouping (default: "age_bin")
#' @return A ggplot object containing the boxplot with points
#' @details Uses ggplot2 with theme_prism() styling. Points are dodged to align with 
#'          their respective boxes. Default colors are darkred and green for the fill variable levels.
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

#' Create Variance Explained Barplot
#' 
#' Creates a barplot showing the proportion of variance explained by a specific 
#' factor across different cell types.
#' 
#' @param mofa_model A MOFA model object
#' @param important_factor Character string specifying the factor to analyze
#' @return A ggplot object containing the variance explained barplot
#' @details Extracts R2 values from the MOFA model cache for the specified factor 
#'          and creates a barplot using theme_prism() styling. The x-axis shows 
#'          cell types and the y-axis shows R2 values.
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

#' Create Shared Genes Plot
#' 
#' Generates a horizontal barplot showing the distribution of genes shared across 
#' different numbers of cell types.
#' 
#' @param shared_genes_table Data frame containing counts of genes shared across 
#'        different numbers of cell types
#' @return A ggplot object containing the horizontal barplot
#' @details Creates a horizontal barplot using theme_minimal() styling. The y-axis shows 
#'          the number of cell types, and the x-axis shows the count of genes shared 
#'          across that many cell types. The plot includes custom margins and font sizes.
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