library(rhdf5)
library(Matrix)
library(glue)
library(scales)
library(tidyverse)
library(dplyr)
library(scuttle)
# Open the connection
my_h5_file <- "./data/eae/cd4/CD4_EAE_ONLY.h5"

# Each cell gets a barcode nucleotide sequence.
# e.g. "AAACCTGAGACCGGAT-1"
my_barcodes <- as.character(rhdf5::h5read(my_h5_file, "matrix/barcodes"))
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
genes <- rhdf5::h5read(my_h5_file, "matrix/features")
genes <- tibble(ensembl_id = genes$id, symbol = genes$name) %>%
    filter(ensembl_id != ".")

# We add the MGI.symbol to the count matrix
counts_matrix <- data.frame(2^as.matrix(counts) - 1) %>%
    cbind("ensembl_id" = rownames(.), .) %>% 
    left_join(x = genes, y = ., by = "ensembl_id")
colnames(counts_matrix)[2] <- "MGI.symbol"

counts_matrix_noNA <- na.omit(counts_matrix)

to_keep <- rowSums(counts_matrix_noNA[, -c(1:2)]) > 0
# calculate TPM counts to feed to compass
counts_TPM <- scuttle::calculateTPM(x = counts_matrix_noNA[to_keep, -c(1:2)])
# save linear gene expression mtrix
write.table(x = cbind(
    "MGI.symbol" = toupper(counts_matrix_noNA$MGI.symbol[to_keep]), 
    counts_TPM), 
    file =  "./data/eae/cd4/tpm_counts_normalized/linear_gene_expression_matrix.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
#  retrieve annotations for cell metadata file
annotations <- read.csv(file = "./data/eae/cd4/CD4_EAE_ONLY_annotations.csv", sep = "\t")
annotations$Cell.name <- paste0("X", annotations$Cell.name)
paste0(colnames(annotations), collapse = "','")
cell_metadata <- setNames(
    annotations[, 
        c("Cell.name", 
        "Stage", 
        'CD4.T.cell.subsets', 
        'Expressed.genes',
        'Mitochondrial.reads.percent',
        'Ribosomal.reads.percent',
        'Total.count')], 
    nm = c("cell_id", 
        "stage", 
        'cd4_cell_subset', 
        'expressed_genes',
        'mitochondrial_reads_percent',
        'ribosomal_reads_percent',
        'total_count'))

cell_metadata$cd4_cell_subset <- ifelse(
    cell_metadata$cd4_cell_subset == "N/A", 
    NA, 
    cell_metadata$cd4_cell_subset)

rownames(cell_metadata) <- cell_metadata$cell_id
# save cell metadata file
write_delim(x = cell_metadata, 
    file = "./data/eae/cd4/tpm_counts_normalized/cell_metadata.csv", 
    delim = ",", quote = "all")

reductions_PCA <- read.csv(file = "./data/eae/cd4/CD4_EAE_ONLY_PCA_coordinates.csv", sep = "\t")
reductions_UMAP <- read.csv(file = "./data/eae/cd4/CD4_EAE_ONLY_UMAP_coordinates.csv", sep = "\t")
reductions <- cbind(
    setNames(data.frame(t(reductions_PCA[1:5,-1])), reductions_PCA[1:5,1]), 
    setNames(data.frame(t(reductions_UMAP[1:2,-1])), reductions_UMAP[1:2,1]))
colnames(reductions) <- c(paste0(rep("PCA",5), 1:5), "UMAP1", "UMAP2")
reductions$cell_id <- rownames(reductions)

# save cell metadata file
write_delim(x = reductions, 
            file = "./data/eae/cd4/CD4_EAE_ONLY_PCA_UMAP_coordinates.csv", 
            delim = ",", quote = "all")
