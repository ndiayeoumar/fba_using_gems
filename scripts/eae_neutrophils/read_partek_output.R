library(rhdf5)
library(Matrix)
library(glue)
library(scales)
library(tidyverse)
library(dplyr)
library(scuttle)
# Open the connection
my_h5_file <- "./data/eae/neutrophils/Neutrophils_EAE_batched.h5"

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

counts_TPM <- calculateTPM(x = counts_matrix_noNA[, -c(1:2)])
write.table(x = cbind(
    "MGI.symbol" = toupper(counts_matrix_noNA$MGI.symbol), 
    counts_TPM), 
    file =  "./data/eae/neutrophils/tpm_counts_normalized/linear_gene_expression_matrix1.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
