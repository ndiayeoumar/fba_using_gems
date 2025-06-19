library(phyloseq)
library(HotLoadings)

source("../lib/fast_mixomics.R")

plotMixOmics <- function (compass_data, clusterA, clusterB, phen_to_map, samples_name = "cell_id") {
  # get counts and metadata
  clA_index <- which(compass_data$cell_metadata$membership == clusterA)
  clB_index <- which(compass_data$cell_metadata$membership == clusterB)
  cluster_metadata <- as.data.frame(compass_data$cell_metadata[c(clA_index, clB_index), ])
  cluster_counts <- compass_data$reaction_consistencies[, c(clA_index, clB_index)]
  
  # compute operational taxanomic unit abundance
  otu <- otu_table(cluster_counts, taxa_are_rows = TRUE)
  # compute sample data
  sd <- sample_data(cluster_metadata)
  # get sample names
  sample_names(otu) <- sample_names(sd)
  # build an experiment-level of phyloseq
  ps <- phyloseq(otu_table = otu, sample_data = sd)
  
  tt <- matrix(rownames(cluster_counts), ncol = 1)
  rownames(tt) <- tt
  taxa_names(ps) <- tt
  tax_table(ps) <- tt
  
  
  # run fast micOmics
  run_cluster <- fast_mixomics(counts = cluster_counts,
                               metadata = cluster_metadata,
                               Y_name = phen_to_map,
                               sample_name = samples_name)
  
  # plot 
  p_grid <- cowplot::plot_grid(
    HotLoadings.plot_loadings_simple(PSOBJ = ps, 
                                     format = "none",
                                     data.splsda = run_cluster$splsda_data,
                                     Y_name = phen_to_map,
                                     comp = 1, 
                                     ndisplay = run_cluster$nfeatures[1]) + 
      ggtitle(paste0("sPLS-DA - EAE Disease, clusters 1 vs 3 (Accuracy = ", 
                     round((1 - run_cluster$error$Max_dist.[1])*100, 2), "%)")),
    HotLoadings.heat_map(PSOBJ = ps,
                         format = "none",
                         data.splsda = run_cluster$splsda_data,
                         Y_name = phen_to_map,
                         component = 1,
                         n_top = run_cluster$nfeatures[1],
                         sample_name = samples_name,
                         facet_formula = ~ get(phen_to_map)), 
    ncol = 2, align = "h", axis = "tb"
  )
  
  res <- list("cluster_counts" = cluster_counts,
              "cluster_metadata" = cluster_metadata,
              "phylo_seq" = ps,
              "fast_mixomics" = run_cluster,
              "plot_grid" = p_grid)
}




