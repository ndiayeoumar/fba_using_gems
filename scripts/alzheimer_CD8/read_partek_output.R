library(rhdf5)
library(Matrix)
library(glue)
library(scales)
library(tidyverse)
library(dplyr)
library(scuttle)
# Open the connection
my_h5_file <- "./data/alzheimer/cd8/Matrix_exported.h5"

# Each cell gets a barcode nucleotide sequence.
# e.g. "AAACCTGAGACCGGAT-1"
my_barcodes <- as.character(h5read(my_h5_file, "matrix/barcodes"))
# Read the data into a dgCMatrix (sparse matrix).
h5 <- h5read(my_h5_file, "matrix")
counts <- sparseMatrix(
    dims = h5$shape,
    i = as.numeric(h5$indices),
    p = as.numeric(h5$indptr),
    x = as.numeric(h5$data),
    index1 = FALSE
)
colnames(counts) <- my_barcodes
rownames(counts) <- as.data.frame(h5[["features"]])$id

# Compute the sparsity of the matrix.
sparsity <- 1 - length(counts@x) / (as.numeric(nrow(counts)) * ncol(counts))
n_genes <- nrow(counts)
n_cells <- ncol(counts)
print(glue("Counts matrix has {comma(n_genes)} genes and {comma(n_cells)} cells"))
print(glue("Counts matrix is {signif(100 * sparsity, 3)}% sparse"))

# We'll also want to get the corresponding gene information from the .h5 file.
genes <- h5read(my_h5_file, "matrix/features")
genes <- tibble(ensembl_id = genes$id, symbol = genes$name) %>%
    filter(ensembl_id != ".")

# We add the MGI.symbol to the count matrix
counts_matrix <- data.frame(2^as.matrix(counts) - 1) %>%
    cbind("ensembl_id" = rownames(.), .) %>% 
    left_join(x = genes, y = ., by = "ensembl_id")
colnames(counts_matrix)[2] <- "MGI.symbol"

counts_matrix_noNA <- na.omit(counts_matrix)

to_keep <- rowSums(counts_matrix_noNA[, -c(1:2)]) > 0

counts_TPM <- calculateTPM(x = counts_matrix_noNA[to_keep, -c(1:2)])
write.table(x = cbind(
    "MGI.symbol" = toupper(counts_matrix_noNA$MGI.symbol[to_keep]), 
    counts_TPM), 
    file =  "./data/alzheimer/cd8/tpm_counts_normalized/linear_gene_expression_matrix2.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#  retrieve annotations for cell metadata file
annotations <- read.csv(file = "./data/alzheimer/cd8/exported Annotations.txt", sep = "\t")
cell_metadata <- setNames(
    annotations[, 
                c("Cell.name", 
                  "Tissue", 
                  'Condition', 
                  'Arms.trajectory.for.violins',
                  'CD8..T.CELLS.SUBSETS.CLASSIFICATION',
                  'Expressed.genes',
                  'Ribosomal.reads.percent',
                  'Total.count')], 
    nm = c("cell_id", 
           "tissue", 
           'condition', 
           'trajectories',
           'cell_type',
           'expressed_genes',
           'ribosomal_reads_percent',
           'total_count'))
# Two cells are duplicated: "TAACACGGTGTCTTGA-1" "GCCAGGTTCTCTCTAA-1"
# rownames(cell_metadata) <- cell_metadata$cell_id

cell_metadata$cell_id[which(duplicated(cell_metadata$cell_id))]
# save cell metadata file
write_delim(x = cell_metadata, 
            file = "./data/alzheimer/cd8/tpm_counts_normalized/cell_metadata.csv", 
            delim = ",", quote = "all")

# Retrieve coordinates
reductions_PCA <- read.csv(file = "./data/alzheimer/cd8/exported PCA Coordinates.txt", sep = "\t")
reductions_UMAP <- read.csv(file = "./data/alzheimer/cd8/exported UMAP All samples together.txt", sep = "\t")
reductions <- cbind(
    setNames(data.frame(t(reductions_PCA[1:5,-1])), reductions_PCA[1:5,1]), 
    setNames(data.frame(t(reductions_UMAP[1:2,-1])), reductions_UMAP[1:2,1]))
colnames(reductions) <- c(paste0(rep("PCA",5), 1:5), "UMAP1", "UMAP2")
reductions$cell_id <- gsub(rownames(reductions), pattern = "[.]", replacement = "-")

# save cell metadata file
write_delim(x = reductions, 
            file = "./data/alzheimer/cd8/PCA_UMAP_coordinates.csv", 
            delim = ",", quote = "all")
