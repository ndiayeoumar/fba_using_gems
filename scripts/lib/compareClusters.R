library(plotly)
library(DT)

# source("../lib/my_compassR/my_compass.R")

#' @description uses wilcoxon rank-sum test to see how much different 2 clusters are
#' @param compass_data compass data instance
#' @param cl1 indexes of the first cluster
#' @param cl2 indexes of the second cluster
#' @param significant threshold of significance
#' @param effect threshold of effect size
#' @param significant_field field that determines significance
#' @param effect_field field that determines effect size
#' @param top top x findings with highest significance
#'
#' @return a list containing a table of cluster numerosity, an interactive table with data about
#' Wilcox test and compass data, the top reactions and also a volcano plot

compareClusters <- function(compass_data, cl1, cl2, 
                           significant_field = "adjusted_p_value", significant = 0.05,
                           effect_field = "avgDiff", effect = 0.5, top = 10) {
  # extracting cell indexes for the two clusters
  cl_index_1 <- which(compass_data$cell_metadata$membership %in% cl1)
  cl_index_2 <- which(compass_data$cell_metadata$membership %in% cl2)
  # creating contingency matrix
  t_numerosity <- data.frame("cluster" = c(paste0("cl_", cl1), 
                                           paste0("cl_", paste0(cl2, collapse = "_"))),
                             "number of cells" = c(length(cl_index_1), length(cl_index_2)))
  colnames(t_numerosity) <- c("cluster", "number of cells")
  
  # performing wilcox test
  group_A <- compass_data$cell_metadata[cl_index_1, ] %>%
    pull(cell_id)
  group_B <- compass_data$cell_metadata[cl_index_2, ] %>%
    pull(cell_id)
  
  wilcoxon_results <- my_wilcoxon_test(
    consistencies_matrix = compass_data$reaction_consistencies[, c(cl_index_1, cl_index_2)],
    settings = compass_data$settings,
    group_A_cell_ids = group_A,
    group_B_cell_ids = group_B
  )
  
  wilcoxon_results$reaction_no_direction <- getReactionNoDirection(wilcoxon_results$reaction_id)
  allReactionMetadata <- getAllReactionMetadata(wilcoxon_results, compass_data)
  
  # vulcano plot
  diff_df <- allReactionMetadata$reactionLevelMetadata[c("reaction_no_direction", "reaction_id", "reaction_name", "subsystem", "cohens_d", "avgDiff", "adjusted_p_value", "p_value")]
  # remove NA values
  diff_df <- na.omit(diff_df)
  interactive_table <- DT::datatable(diff_df, caption = "example of DT", filter="top", extensions = 'Buttons',
                                 options = list(dom = 'Blfrtip',
                                                buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                                lengthMenu = list(c(10,25,50,-1), c(10,25,50,"All"))))
  
  # add a grouping column; default value is "not significant"
  diff_df["group"] <- "NotSignificant"
  
  # for our plot, we want to highlight 
  # significant_field < 0.05 and effect_field > 0.5
  
  # change the grouping for the entries with significance but not a large enough effect size change
  diff_df[which(diff_df[significant_field] < significant & abs(diff_df[effect_field]) < effect ),"group"] <- paste0("Significant (", significant_field, "<", significant, ")")
  
  # change the grouping for the entries a large enough effect size change but not a low enough significance
  diff_df[which(diff_df[significant_field] > significant & abs(diff_df[effect_field]) > effect ),"group"] <- paste0("Effect size (|", effect_field, "|>", effect, ")")
  
  # change the grouping for the entries with both significance and large enough effect size change
  diff_df[which(diff_df[significant_field] < significant & abs(diff_df[effect_field]) > effect ),"group"] <- "Significant & Effect size"
  
  # Find the top peaks..
  half_top <- top / 2
  top_peaks <- diff_df[with(diff_df, order(get(effect_field), get(significant_field))),][1:half_top, ]
  top_peaks <- rbind(top_peaks, diff_df[with(diff_df, order(-get(effect_field), get(significant_field))),][1:half_top, ])
  
  # make the Plot.ly plot
  p_volcano <- plot_ly(data = diff_df,
          x = diff_df[, effect_field],
          y = -log10(diff_df[, significant_field]),
          text = diff_df$reaction_name,
          mode = "markers",
          color = diff_df$group,
          type = "scatter") %>% 
    layout(title ="Volcano Plot")
  
  return(list("cluster_numerosity"=t_numerosity, "stat_table"=interactive_table, "metadata" =diff_df, "top_findings"=top_peaks, "plot"=p_volcano))
}
