library(compassR)
library(ggrepel)
library(tidyverse)

source("../utility.R")

# params: 
  # filePath: path to where all the needed files are
  # gene_identifiers: gene id system
  # phenotype: characteristic used to identify groups to which cells belongs to
  # phenoA: phenotype/characteristic of group A
  # phenoB: phenotype/characteristic of group B
# returns a list with:
  # wilcoxon_data: data frame containing all the fields that may be needed

getCompassData <- function(filePath, gene_identifiers, phenotype, phenoA, phenoB) {
  ###########################################
  ########    LOADING COMPASS     ###########
  ###########################################
  # settings of compass
  compass_settings <- CompassSettings$new(
    user_data_directory = directoryPath,
    cell_id_col_name = "cell_id",
    gene_id_col_name = gene_identifiers
  )
  
  # create instance of compassData
  compass_data <- CompassData$new(compass_settings)
  # instance of compassAnalyzer
  compass_analyzer <- CompassAnalyzer$new(compass_settings)
  
  
  ###########################################
  ########    WILCOXON TEST       ###########
  ###########################################
  
  # defining cell groups based on phenotype
  # GROUP A is for cells with chronic stage
  group_A_cell_ids <-
    compass_data$cell_metadata %>%
    filter(get(phenotype) == phenoA) %>%
    pull(cell_id)
  # GROUP B is for cells with onset stage
  group_B_cell_ids <-
    compass_data$cell_metadata %>%
    filter(get(phenotype) == phenoB) %>%
    pull(cell_id)
  
  # get wilcoxon test
  wilcoxon_results <- compass_analyzer$conduct_wilcoxon_test(
    # using reaction_consistencies
    # alternative: metareaction_consistencies in case of with metacells
    compass_data$reaction_consistencies,
    group_A_cell_ids,
    group_B_cell_ids,
    # put to true when analysing data from metacells
    for_metareactions = FALSE
  )
  
  # add reactions ids without the direction
  wilcoxon_results$reaction_no_direction <- getReactionNoDirection(wilcoxon_results$reaction_id)
  
  ###########################################
  ###########   UMAP COMPONENTS   ###########
  ###########################################
  
  cell_info_with_umap_components <-
    compass_analyzer$get_umap_components(
      compass_data$reaction_consistencies
    ) %>%
    inner_join(
      compass_data$cell_metadata,
      by = "cell_id"
    ) %>%
    left_join(
      compass_data$gene_expression_statistics,
      by = "cell_id"
    )
  # components are vectors of chars and we need to transform them into numerics
  cell_info_with_umap_components$component_1 <- as.numeric(cell_info_with_umap_components$component_1)
  cell_info_with_umap_components$component_2 <- as.numeric(cell_info_with_umap_components$component_2)
  
  
  ###########################################
  #########   MERGING THE OUTPUTS   #########
  ###########################################
  # list containing all needed metadata of the reactions
  allReactionMetadata <- getAllReactionMetadata(wilcoxon_results, compass_data)
  # list containing all needed metadata on the cells
  allCellMetadata <- getAllCellMetadata(compass_data)
  
  globalData <- list(wilcoxon_results,  cell_info_with_umap_components, allReactionMetadata, allCellMetadata)
  names(globalData) <- c("wilcoxon_results",  "cell_info_with_umap_components", "allReactionMetadata", "allCellMetadata")
  
  return(globalData)
}
