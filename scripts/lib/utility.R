library(httr)
library(jsonlite)
library(scuttle)
library("ClusterR")
library("cluster")

##########################################
############   pre-processing  ###########
##########################################

get_cpm_counts <- function(sce, n = dim(sce)[2], colnames_metadata) {
  sample_nozero <- sce[rowSums(counts(sce)) > 0, 1:n]
  cpmcounts <- as.matrix(cpm(sample_nozero))
  metadata <- data.frame(colData(sample_nozero)[, colnames_metadata])
  return(list("cpm" = cpmcounts, "metadata" = metadata))
}

get_tpm_counts <- function(sce, n = dim(sce)[2], colnames_metadata) {
  sample_nozero <- sce[rowSums(counts(sce)) > 0, 1:n]
  tpmcounts <- as.matrix(tpm(sample_nozero))
  metadata <- data.frame(colData(sample_nozero)[, colnames_metadata])
  return(list("tpm" = tpmcounts, "metadata" = metadata))
}

get_counts <- function(sce, n = dim(sce)[2], colnames_metadata) {
  sample_nozero <- sce[rowSums(counts(sce)) > 0, 1:n]
  counts <- as.matrix(counts(sample_nozero))
  metadata <- data.frame(colData(sample_nozero)[, colnames_metadata])
  return(list("counts" = counts, "metadata" = metadata))
}

writeCompassInput <- function(sample, countType, output_dir) {
  write.table(x = data.frame("Symbol" = rownames(sample[[countType]]), 
                             sample[[countType]]), 
              file = paste0(output_dir, "linear_gene_expression_matrix.tsv"), 
              sep = "\t", row.names = FALSE)
# write.table(data.frame("H"=rownames(a),a),"a.txt", row.names=FALSE)
  write.csv(data.frame("cell_id"=rownames(sample$metadata), sample$metadata),
            file = paste0(output_dir, "cell_metadata.csv"),
            row.names=FALSE)
}

###########################################
############  post-processing  ############
###########################################

## params:
  # reactionIds: vector of all reaction ids in wilcoxon_resuls object
## returns:
  # clean_ids: vector of the ids in input without the direction ("_neg"/"_pos")
getReactionNoDirection <- function(reactionIds) {
  clean_ids  <- gsub(x = reactionIds, pattern = "_pos|_neg", replacement = "")
  return(clean_ids)
}


## params: 
  # wilcoxon_results: result of the wilcoxon rank sum test
  # compass_data: metadata from our input file
  # columnTargets: vector of the columns names we want to extract
## return a list:
  # reactionLevelMetadata: data frame containing all info about the reactions
  # reactionConsistencies: data frame where each row is a reaction and each column is a cell
getAllReactionMetadata <- function(wilcoxon_results, compass_data) {
  reactionLevelMetadata <- merge(
    x = wilcoxon_results,
    y = compass_data$reaction_metadata,
    by.x = "reaction_no_direction",
    by.y = 'reaction_no_direction'
  )

  reactionConsistencies <- compass_data$reaction_consistencies
  results <- list(reactionLevelMetadata, reactionConsistencies)
  names(results) <- c("reactionLevelMetadata", "reactionConsistencies")
  return(results)
}


#' @param compass_data metadata from our input filet
#' @teturn data frame containing all info about the cells
getAllCellMetadata <- function(compass_data) {
  allCellMetadata <- merge(
    x = compass_data$cell_metadata,
    y = compass_data$gene_expression_statistics,
    by.x = "cell_id",
    by.y = 'cell_id'
  )
  return(allCellMetadata)
}


#' @param compass_data metadata from our input filet
#' @param gene_identifiers gene id nomenclature system
#' @return data frame containing all info about the genes
getAllGeneMetadata <- function(compass_data, gene_identifiers) {
  allGeneMetadata <- merge(
    x = compass_data$gene_metadata,
    y = compass_data$metabolic_genes,
    by.x = gene_identifiers,
    by.y = "gene"
  )
  return(allGeneMetadata)
}

#' @param matrix the matrix we want to clusterize
#' @param comps a vector of which components to focus on
#' @param ncluster indicates how many clusters we want to make
#' @return a vector where each position indicates the group number to which each cell belongs to
get_cluster <- function(matrix, comps = c(1,2), ncluster){
  d <- dist(matrix[,comps])
  hc <- hclust(d, method = "complete")
  membs <- cutree(hc, k = ncluster)
  return(membs)
}


# function to compute average silhouette for k clusters
avg_sil <- function(df, k) {
  # km.res <- kmeans(df, centers = k, nstart = 25)
  d <- dist(as.matrix(df))
  hc.res <- hclust(d, method = "average")
  ss <- silhouette(cutree(hc.res, k), d)
  mean(ss[, 3])
}

#' @param umap
#' @param

numClustersK <- function (umap, k.values = 2:15, comps = c(1, 2), plot = TRUE) {
  # remove possible NA
  umap <- na.omit(umap[, comps])
  # extract avg silhouette for 2-15 clusters
  avg_sil_values <- map_dbl(k.values, avg_sil, df = umap)
  df_to_plot <- data.frame("n_cluster" = k.values, 
                           "avg_silhouette" = avg_sil_values)
  if (plot) {
      p <- ggplot(df_to_plot, aes(x = n_cluster, y = avg_silhouette)) + 
          geom_line() + 
          geom_point() + 
          ggtitle("Average Silhouette indexes", 
                  subtitle = "Hierarchical clustering on euclidean distances") +
          xlab("Number of clusters") + ylab("Average Silhouette")
      print(p)
  }

  optimalClusterNum <- which(avg_sil_values == max(avg_sil_values))
  
  return (k.values[optimalClusterNum])
}
