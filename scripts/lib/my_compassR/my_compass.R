# install and upload libraries
library(magrittr)
library(R6)
library(tidyverse)
library(uwot)
library(dplyr)

#' @description
#' define the settings parameters used by compassR
#'
#' @param input_path path to the folder containing expression metadata and reaction scored matrices
#' @param metadata_path compassR's rebuilt metadata
#' @param gene_id_col_name "MGI.symbol" for mice and "HGNC.symbol" for human
#'
#' @return a list with all the settings parameters

setCompassRSettings <- function (input_path, metadata_path, gene_id_col_name){
  if (!is.element(gene_id_col_name, c("MGI.symbol", "HGNC.symbol"))) {
    stop('name of the gene id  must be  either "MGI.symbol" (for mice metadata) or "HGNC.symbol" (for human)')
  }
  settings <- list()
  settings$gene_metadata_path <- paste0(metadata_path, "gene_metadata.csv")
  settings$metabolite_metadata_path <-  paste0(metadata_path, "metabolite_metadata.csv")
  settings$reaction_metadata_path <-  paste0(metadata_path, "reaction_metadata.csv")
  settings$cell_metadata_path <- paste0(input_path, "cell_metadata.csv")
  settings$reaction_scores_path <- paste0(input_path, "reactions.tsv")
  settings$gene_expression_path <-paste0(input_path, "linear_gene_expression_matrix.tsv")
  
  settings$min_reaction_consistency <- 1e-10
  settings$min_reaction_range <- 1e-3
  settings$reaction_direction_separator <- "_"
  settings$reaction_directions <- c("neg", "pos") 
  settings$cluster_strength <- 0.1
  settings$gene_id_col_name <- "MGI.symbol"    # for human cells: HGNC.symbol
  settings$cell_id_col_name = "cell_id"
  
  return(settings)
}


#' @description
#' drops the reactions that that are inconsistent in at least one cell
#'
#' @param reaction_consistencies A param.
#' @param min_consistency A param.
#'
#' @return An output.
#'
#' @noRd
#' 
drop_inconsistent_reactions <- function(reaction_consistencies, min_consistency) {
  reactions_inconsistent_for_all_cells <-  apply(reaction_consistencies, 1, function(x) ifelse(max(x) - min(x) < min_consistency, TRUE, FALSE))
  alert_of_drop(
    reactions_inconsistent_for_all_cells,
    "inconsistent for all cells"
  )
  reaction_consistencies <- reaction_consistencies[!reactions_inconsistent_for_all_cells,]
  return(reaction_consistencies)
}


#' @description
#' get the reactions consistencies which are calculated from the compass scores and indicate
#' the likely hood of a each reaction to occur in a certain cell
#'
#' @param reaction_scores A param.
#' @param min_consistency A param.
#' @param min_range A param.
#'
#' @return An output.
#'
#' @importFrom magrittr %>% %<>%
#'
#' @noRd
get_reaction_consistencies <- function(reaction_scores, min_range) {
  raw_consistencies <- reaction_scores %>% 
      (function(x) -log1p(x))
  all_consistencies <- drop_inconsistent_reactions(
                                reaction_consistencies = raw_consistencies,
                                min_consistency = min_range)
  reaction_consistencies <- all_consistencies - min(all_consistencies)
  return(reaction_consistencies)
}


#' @description
#' calculates the effect size used to indicate the standardised difference between two means.
#'
#' @param group_A_values A param.
#' @param group_B_values A param.
#'
#' @return An output.
#'
#' @noRd
cohens_d <- function(group_A_values, group_B_values) {
  n_A <- length(group_A_values)
  n_B <- length(group_B_values)
  mu_A <- mean(group_A_values)
  mu_B <- mean(group_B_values)
  var_A <- var(group_A_values)
  var_B <- var(group_B_values)
  pooled_sd <- sqrt(((n_A - 1) * var_A + (n_B - 1) * var_B) / (n_A + n_B - 2))
  cohen_d <- (mu_A - mu_B) / pooled_sd
  return(cohen_d)
}



#' @description
#' calculates the difference between two means.
#'
#' @param group_A_values A param.
#' @param group_B_values A param.
#'
#' @return An output.
#'
#' @noRd
avgDiff <- function(group_A_values, group_B_values) {
  mu_A <- mean(group_A_values)
  mu_B <- mean(group_B_values)
  logFC <- (mu_A - mu_B)
  return(logFC)
}




#' @description
#' drops reactions with very small range of values through the different cells
#'
#' @param reaction_consistencies A param.
#' @param min_range A param.
#'
#' @return An output.
#'
#' @noRd
drop_constant_reactions <- function(reaction_consistencies, min_range) {
  constant_reactions <- apply(reaction_consistencies, 1, function(x) { max(x) - min(x) < min_range })
  # alert_of_drop(
  #   constant_reactions,
  #   "have too small a range of values"
  # )
  reaction_consistencies <- reaction_consistencies[!constant_reactions,]
  reaction_consistencies
}


#' @description
#' Description.
#'
#' @param ids A param.
#' @param separator A param.
#' @param annotations A param.
#' @param id_col_name A param.
#' @param unannotated_col_name A param.
#' @param annotation_col_name A param.
#'
#' @return An output.
#'
#' @importFrom magrittr %>% %<>%
#'
#' @noRd
get_annotations <- function(ids, separator, annotations, id_col_name, unannotated_col_name, annotation_col_name) {
  if (is.null(separator) | is.null(annotations)) {
    annotations <- cbind(ids, NA, NA)
  } else {
    patterns <- stringr::str_glue(
      "^(.+){separator}({annotation})$",
      separator = separator,
      annotation = annotations
    )
    annotations <- purrr::map_dfr(
      ids,
      function(id) {
        stringr::str_match(id, patterns) %>%
          na.omit() %>%
          as.data.frame(stringsAsFactors = FALSE)
      }
    )
  }
  colnames(annotations) <- c(
    id_col_name,
    unannotated_col_name,
    annotation_col_name
  )
  annotations <- tibble::as_tibble(annotations)
  return(annotations)
}




#' @description
#' Description.
#'
#' @param linear_gene_expression matrix of gene (rows) expression in each cell (column)
#' @param metabolic_genes Each row describes a gene in terms of its ID and whether it's a metabolic gene.
#'
#' @return Each row describes a cell in terms of its ID, total expression, metabolic expression, and metabolic activity.
#'
#' @importFrom magrittr %>% %<>%
#'
#' @noRd
get_gene_expression_statistics <- function(linear_gene_expression, metabolic_genes) {
  total_expressions <- colSums(linear_gene_expression)
  metabolic_expressions <- colSums(linear_gene_expression[
    metabolic_genes %>%
      dplyr::filter(is_metabolic == TRUE) %>%
      dplyr::pull(gene),
  ])
  gene_expression_statistics <-
    rbind(
      total_expression = total_expressions,
      metabolic_expression = metabolic_expressions,
      metabolic_activity = metabolic_expressions / total_expressions
    ) %>%
    t() %>%
    tibble::as_tibble(rownames = "cell_id")
  gene_expression_statistics
}


#' @description
#' Description.
#'
#' @param reaction_consistencies matrix of compass relaborated compass scrores
#' @param metareactions A param.
#'
#' @return Each row is a metareaction and each column is a cell. metareaction_consistencies[i, j] is the consistency (or "compatibility") between metareaction i and cell j
#'
#' @importFrom magrittr %>% %<>%
#'
#' @noRd
  get_metareaction_consistencies <- function(reaction_consistencies, metareactions) {
  metareaction_consistencies <-
    reaction_consistencies %>%
    tibble::as_tibble(rownames = "reaction_id") %>%
    dplyr::inner_join(metareactions, by = "reaction_id") %>%
    dplyr::select(-reaction_id) %>%
    dplyr::group_by(metareaction_id) %>%
    dplyr::summarize_all(mean) %>%
    dplyr::ungroup() %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("metareaction_id") %>%
    data.matrix()
  
  return(metareaction_consistencies)
}

#' @description
#' Description.
#'
#' @param reaction_consistencies A param.
#' @param cluster_strength A param.
#'
#' @return An output.
#'
#' @importFrom magrittr %>% %<>%
#'
#' @noRd
get_metareactions <- function(reaction_consistencies, cluster_strength) {
  pairwise_reaction_correlations <- cor(t(reaction_consistencies), method = "spearman")
  pairwise_reaction_distances <- as.dist(1 - pairwise_reaction_correlations)
  reaction_hierarchy <- hclust(pairwise_reaction_distances, method = "complete")
  metareactions <-
    cutree(reaction_hierarchy, h = cluster_strength) %>%
    tibble::enframe(name = "reaction_id", value = "metareaction_code") %>%
    dplyr::transmute(
      reaction_id = reaction_id,
      metareaction_id = paste("group", metareaction_code, sep = "_")
    )
  return(metareactions)
}


#' @description
#' Description.
#'
#' @param file_path A param.
#' @param column_specification A param.
#'
#' @return An output.
#'
#' @noRd
read_table <- function(file_path, column_specification) {
  file_reader <- get_file_reader(file_path)
  data <- file_reader(
    file_path,
    col_types = column_specification,
    na = c("", "NA", "N/A", "N.A.", "na", "n/a", "n.a.")
  )
  data
}

#' @description
#' Description.
#'
#' @param file_path A param.
#'
#' @return An output.
#'
#' @noRd
get_file_reader <- function(file_path) {
  if (endsWith(file_path, ".csv") | endsWith(file_path, ".csv.gz")) {
    file_reader <- readr::read_csv
  } else if (endsWith(file_path, ".tsv") | endsWith(file_path, ".tsv.gz")) {
    file_reader <- readr::read_tsv
  } else {
    stop(
      stringr::str_glue("File \"{file_path}\" has an unsupported file extension."),
      call. = FALSE
    )
  }
  file_reader
}

#' @description
#' Description.
#'
#' @param file_path A param.
#'
#' @return An output.
read_compass_metadata <- function(file_path) {
  data <- read_table(file_path, readr::cols(.default = "c"))
  return (data)
}

#' @description
#' Description.
#'
#' @param file_path A param.
#' @param index A param.
#' @param suppress_warnings A param.
#'
#' @return An output.
read_compass_matrix <- function(file_path, index, suppress_warnings = FALSE) {
  warning_handler <- if (suppress_warnings) { suppressWarnings } else { function(expr) { expr } }
  data <-
    warning_handler(read_table(file_path, readr::cols())) %>%
    dplyr::rename(!!index := 1) %>%
    tibble::column_to_rownames(index) %>%
    as.data.frame()
  return(data)
}

#' @description
#' Description.
#'
#' @param package_name A param.
#'
#' @return An output.
#'
#' @noRd
require_suggested_package <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    stop(
      stringr::str_glue(
        "Package \"{package_name}\" is required for this function. Please install it."
      ),
      call. = FALSE
    )
  }
}

#' @description
#' Description.
#'
#' @param reactions_to_drop A param.
#' @param description A param.
#' @param is_warning A param.
#'
#' @return An output.
#'
#' @noRd
alert_of_drop <- function(reactions_to_drop, description, is_warning = FALSE) {
  if (any(reactions_to_drop)) {
    num_reactions_to_drop <- sum(reactions_to_drop)
    alert <- ifelse(is_warning, warning, message)
    alert(stringr::str_glue(
      "Dropping {num_reactions_to_drop} reactions that are {description} ..."
    ))
  }
}

#' @description
#' Description.
#'
#' @param string A param.
#' @param indentation_level A param.
#' @param indentation_style A param.
#'
#' @return An output.
#'
#' @noRd
indent <- function(string, indentation_level = 1, indentation_style = "  ") {
  indented_string <- paste(
    strrep("  ", indentation_level),
    string,
    sep = ""
  )
  indented_string
}

#' @description
#' Description.
#'
#' @param binding_name A param.
#' @param binding_value A param.
#' @param separator A param.
#'
#' @return An output.
#'
#' @importFrom magrittr %>% %<>%
#'
#' @noRd
get_binding_representation <- function(binding_name, binding_value, separator = ", ") {
  if (length(binding_value) > 1) {
    binding_value %<>% paste(collapse = separator)
  }
  binding_representation <- stringr::str_glue(
    "{binding_name}: {binding_value}"
  )
  binding_representation
}

#' @description
#' Description.
#'
#' @param table A param.
#' @param table_name A param.
#' @param table_class A param.
#' @param rows A param.
#' @param cols A param.
#'
#' @return An output.
#'
#' @noRd
get_tabular_data_representation <- function(table, table_name, table_class, rows, cols) {
  tabular_data_representation <- stringr::str_glue(
    "{table_name} {table_class} ({dim(table)[1]} {rows} x {dim(table)[2]} {cols})"
  )
  tabular_data_representation
}

# actual compass object build
build_my_compass_data <- function(settings) {
  gene_metadata <- read_compass_metadata(settings$gene_metadata_path)
  metabolite_metadata <- read_compass_metadata(settings$metabolite_metadata_path)
  reaction_metadata <- read_compass_metadata(settings$reaction_metadata_path)
  cell_metadata <- read_compass_metadata(settings$cell_metadata_path)
  reaction_scores <- read_compass_matrix(settings$reaction_scores_path, "reaction_id", suppress_warnings = TRUE)
  linear_gene_expression <- read_compass_matrix(settings$gene_expression_path, "gene")
  reaction_consistencies <- get_reaction_consistencies(
    reaction_scores,
    min_range = settings$min_reaction_range
  )
  annotated_reactions <- get_annotations(
    rownames(reaction_consistencies),
    separator = settings$reaction_direction_separator,
    annotations = settings$reaction_directions,
    id_col_name = "reaction_id",
    unannotated_col_name = "reaction_no_direction",
    annotation_col_name = "direction"
  )
  metareactions <- get_metareactions(
    reaction_consistencies,
    cluster_strength = settings$cluster_strength
  )
  metareaction_consistencies <- get_metareaction_consistencies(
    reaction_consistencies,
    metareactions
  )
  reaction_partitions <-
    annotated_reactions %>%
    dplyr::left_join(
      metareactions,
      by = "reaction_id"
    )
  metabolic_genes <-
    tibble::tibble(gene = rownames(linear_gene_expression)) %>%
    dplyr::left_join(
      tibble::tibble(
        gene = gene_metadata[[settings$gene_id_col_name]],
        is_metabolic = TRUE
      ),
      by = "gene"
    ) %>%
    tidyr::replace_na(list(
      is_metabolic = FALSE
    ))
  gene_expression_statistics <- get_gene_expression_statistics(
    linear_gene_expression,
    metabolic_genes
  )
  compass_data <- list()
  compass_data$settings <- settings
  compass_data$reaction_consistencies <- reaction_consistencies
  compass_data$metareaction_consistencies <- metareaction_consistencies
  compass_data$metabolic_genes <- metabolic_genes
  compass_data$gene_expression_statistics <- gene_expression_statistics
  compass_data$cell_metadata <- cell_metadata
  compass_data$gene_metadata <- gene_metadata
  compass_data$metabolite_metadata <- metabolite_metadata
  compass_data$reaction_metadata <- reaction_metadata
  compass_data$reaction_partitions <- reaction_partitions
  return(compass_data)
}

#'
#'document later
#' @param 
#' @param 
#' 
#' @return 
build_umap_from_compass <- function(compass_data, num_umap_comps = 2) {
  # buil umap components
  umap_components <- uwot::umap(t(compass_data$reaction_consistencies), n_components = num_umap_comps)
  umap_components <- cbind(colnames(compass_data$reaction_consistencies), umap_components)
  colnames(umap_components) <- append(
    compass_data$settings$cell_id_col_name,
    paste("component", 1:num_umap_comps, sep = "_")
  )
  umap_components <- tibble::as_tibble(umap_components)
  
  umap_components <- umap_components %>%
    inner_join(
      compass_data$cell_metadata,
      by = "cell_id"
    ) %>%
    left_join(
      compass_data$gene_expression,
      by = "cell_id"
    )
  # convert components value back to number
  for (i in 1:num_umap_comps) {
    umap_components[[paste0('component_', i)]] <- as.numeric(umap_components[[paste0('component_', i)]])
  }
  return(umap_components)
}


my_wilcoxon_test <- function(consistencies_matrix, settings, group_A_cell_ids, group_B_cell_ids, for_metareactions = FALSE) {
  # if (0 < length(intersect(group_A_cell_ids, group_B_cell_ids))) {
  #   message("Groups A and B are not mutually exclusive. Continuing anyways ...")
  # }
  group_A_values_per_metareaction <-
    consistencies_matrix %>%
    t() %>%
    tibble::as_tibble(rownames = settings$cell_id_col_name) %>%
    dplyr::right_join(
      tibble::tibble(!!settings$cell_id_col_name := group_A_cell_ids),
      by = settings$cell_id_col_name
    )
  group_A_values_per_metareaction <- as.data.frame(group_A_values_per_metareaction)

  group_B_values_per_metareaction <-
    consistencies_matrix %>%
    t() %>%
    tibble::as_tibble(rownames = settings$cell_id_col_name) %>%
    dplyr::right_join(
      tibble::tibble(!!settings$cell_id_col_name := group_B_cell_ids),
      by = settings$cell_id_col_name
    )
  group_B_values_per_metareaction <- as.data.frame(group_B_values_per_metareaction)

  metareaction_ids <- rownames(consistencies_matrix)
  wilcoxon_results <-
    purrr::map_dfr(
      metareaction_ids,
      function(metareaction_id) {
        group_A_values <- group_A_values_per_metareaction[,metareaction_id]
        group_B_values <- group_B_values_per_metareaction[,metareaction_id]
        wilcoxon_result_obj <- wilcox.test(group_A_values, group_B_values)
        wilcoxon_result_tbl <- data.frame(
          metareaction_id = metareaction_id,
          wilcoxon_statistic = wilcoxon_result_obj$statistic,
          cohens_d = cohens_d(group_A_values, group_B_values),
          avgDiff = avgDiff(group_A_values, group_B_values),
          p_value = wilcoxon_result_obj$p.value,
          stringsAsFactors = FALSE
        )
        wilcoxon_result_tbl
      }
    )
  wilcoxon_results <- tibble::as_tibble(wilcoxon_results) %>%
    dplyr::mutate(adjusted_p_value = p.adjust(dplyr::pull(., p_value), method = "BH"))
  if (!for_metareactions) {
    wilcoxon_results %<>% dplyr::rename(reaction_id = metareaction_id)
  }
  return(wilcoxon_results)
}



#' @param tg vector of our targeted groups
#' @param cell_membership vector where each position indicates the group number to which each cell belongs to
#' @param compass_data compassData istance
pairwise_comparisons <- function(tg, cell_membership, compass_data) {
  tg_list <- list()
  for (i in 1:(length(tg) - 1)) {
    for (j in (i + 1):length(tg)) {
      tg_list <- append(tg_list, list(c(tg[i], tg[j])))
      names(tg_list)[length(tg_list)] <- paste0("Comparison", tg[i], "vs", tg[j])
    }
  }
  compass_data$cell_metadata$cell_membership <- cell_membership
  l_out <- lapply(tg_list, function(x){
    # GROUP A is for cells with phenotype A
    group_A_cell_ids <-
      compass_data$cell_metadata %>%
      filter(cell_membership == x[1]) %>%
      pull(cell_id)
    # GROUP B is for cells with phenotype B
    group_B_cell_ids <-
      compass_data$cell_metadata %>%
      filter(cell_membership == x[2]) %>%
      pull(cell_id)
    
    # get wilcoxon test
    wilcoxon_results <- my_wilcoxon_test(
      # using reaction_consistencies (alt: metareaction consistencies in case of with metacells)
      compass_data$reaction_consistencies,
      compass_data$settings,
      group_A_cell_ids,
      group_B_cell_ids,
      # put to true when analysing data from metacells
      for_metareactions = FALSE
    )
    return(wilcoxon_results)
  })
  return(l_out)
}



#'
#' @param 
#' 
#' may be double but better keep it here too
#' 
get_differential_scores <- function(wilcoxon_results, compass_data, facets) {
  compass_scores_by_cell_type <- wilcoxon_results %>%
    left_join(
      select(compass_data$reaction_partitions, "reaction_id", "reaction_no_direction"),
      by = "reaction_id"
    ) %>%
    left_join(
      compass_data$reaction_metadata,
      by = "reaction_no_direction"
    ) %>%
    # Keep only "confident reactions", as defined in our paper.
    filter(!is.na(EC_number)) %>%
    filter(confidence == "0" | confidence == "4") %>%
    # Exclude non-mitochondrially localized reactions from TCA.
    mutate(subsystem = case_when(
      reaction_id == "SPMDOX_pos" ~ "Arginine and Proline Metabolism",
      subsystem == "Citric acid cycle" & !grepl("[m]", formula, fixed = TRUE) ~ "Other",
      TRUE ~ subsystem
    )) %>%
    # Assign reactions to the appropriate subsystem.
    mutate(
      subsystem_priority = factor(subsystem) %>%
        fct_recode(
          "Glycolysis" = "Glycolysis/gluconeogenesis",
          "TCA cycle" = "Citric acid cycle"
        ) %>%
        fct_collapse("Amino acid metabolism" = c(
          "Alanine and aspartate metabolism",
          "Arginine and Proline Metabolism",
          "beta-Alanine metabolism",
          "Cysteine Metabolism",
          "D-alanine metabolism",
          "Folate metabolism",
          "Glutamate metabolism",
          "Glycine, serine, alanine and threonine metabolism",
          "Histidine metabolism",
          "Lysine metabolism",
          "Methionine and cysteine metabolism",
          "Taurine and hypotaurine metabolism",
          "Tryptophan metabolism",
          "Tyrosine metabolism",
          "Urea cycle",
          "Valine, leucine, and isoleucine metabolism"
        )) %>%
        fct_other(keep = facets) %>%
        fct_relevel(facets)
    ) %>%
    # Keep only the subsystems for which we want to plot a facet.
    filter(subsystem_priority != "Other") %>%
    # Lower-bound the adjusted p-value.
    mutate(adjusted_p_value = if_else(
      subsystem_priority == "Amino acid metabolism" & adjusted_p_value <= 1e-12,
      1e-12,
      adjusted_p_value
    )) %>%
    # Assign descriptive labels to various reactions.
    mutate(label = case_when(
      reaction_id == "PGM_neg" ~ "phosphoglycerate mutase (PGAM)",
      reaction_id == "LDH_L_neg" ~ "lactate dehydrogenase",
      reaction_id == "PDHm_pos" ~ "pyruvate dehydrogenase (PDH)",
      reaction_id == "TPI_neg" ~ "triosephosphate isomerase (DHAP forming)",
      reaction_id == "FACOAL1821_neg" ~ "long-chain fatty-acid-CoA ligase",
      reaction_id == "r1257_pos" ~ "long-chain fatty-acid-CoA ligase",
      reaction_id == "FACOAL1831_neg" ~ "long-chain fatty-acid-CoA ligase",
      reaction_id == "CSNATr_neg" ~ "carnitine O-acetyltransferase",
      reaction_id == "C160CPT1_pos" ~ "carnitine O-palmitoyltransferase",
      reaction_id == "ACONTm_pos" ~ "aconitate hydratase",
      reaction_id == "SUCOASm_pos" ~ "succinate-CoA ligase",
      reaction_id == "AKGDm_pos" ~ "alpha-ketoglutarate dehydrogenase",
      reaction_id == "SUCD1m_pos" ~ "succinate dehydrogenase",
      reaction_id == "ICDHyrm_pos" ~ "isocitrate dehydrogenase",
      reaction_id == "CK_pos" ~ "creatine\nkinase",
      reaction_id == "PGCD_pos" ~ "phosphoglycerate dehydrogenase",
      reaction_id == "ARGSS_pos" ~ "arginosuccinate synthase",
      reaction_id == "r0281_neg" ~ "putrescine diamine oxidase",
      reaction_id == "SPMDOX_pos" ~ "spermidine dehydrogenase (spermidine -> GABA)",
      reaction_id == "ARGDCm_pos" ~ "arginine decarboxylase",
      reaction_id == "AGMTm_pos" ~ "agmatinase",
      reaction_id == "GHMT2r_pos" ~ "serine hydroxymethyltransferase",
      reaction_id == "AHC_pos" ~ "adenosylhomocysteinase",
      reaction_id == "METAT_pos" ~ "methionine adenosyltransferase",
      reaction_id == "METS_pos" ~ "methionine\nsynthase",
      reaction_id == "ARGN_pos" ~ "arginase",
      TRUE ~ ""
    ))
  return(compass_scores_by_cell_type)
}
