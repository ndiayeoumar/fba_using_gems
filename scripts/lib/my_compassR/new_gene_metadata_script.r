library(plyr)
# Per trovare ortologi
# install.packages("babelgene")
library(babelgene)
# Scarico uomo
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("org.Mm.eg.db")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
x <- org.Hs.eg.db
keytypes(x)
columns(x)

HGNC.ENTREZID <- keys(x, "ENTREZID")
human_genes <- select(x, keys = HGNC.ENTREZID, columns = c("SYMBOL") , keytype = "ENTREZID")
mouse_orthologs <- orthologs(genes = human_genes$SYMBOL, species = "mouse")

y <- org.Mm.eg.db
mouse_genes <- select(y, keys = mouse_orthologs$symbol, columns = c("MGI") , keytype = "SYMBOL")
mouse_genes_unique <- ddply(.data = mouse_genes, .variables = ~ SYMBOL, .fun = function(x){
    gsub(x[1, 2], pattern = "MGI:MGI:", replacement = "MGI:")
})

ord <- match(mouse_orthologs$symbol, mouse_genes_unique$SYMBOL)

gene_metadata <- mouse_orthologs[, c("human_entrez", "human_symbol", "symbol")]
gene_metadata <- cbind(gene_metadata, mouse_genes_unique$V1[ord])
colnames(gene_metadata) <- c("EntrezGene.ID", "HGNC.symbol", "MGI.symbol", "MGI.ID")
write.csv(x = gene_metadata, file = "../gene_metadata.csv", quote = FALSE, row.names = FALSE)
