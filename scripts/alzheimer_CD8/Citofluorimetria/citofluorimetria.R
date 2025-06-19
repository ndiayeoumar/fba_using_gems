# Leggo i dati
annotations <- read.delim("../../../Downloads/Annotations.csv", sep = ";")
data <- t(read.delim("../../../Downloads/Data.csv", 
    sep = ";", row.names = 1, check.names = FALSE, dec = ","))
# Accosto annotazioni e dati di citofluorimetria
df <- cbind(data[annotations$X,], annotations[, -1])

# Passo da un formato wide ad un formato long
library(dplyr)
library(reshape)
df_melt <- melt(df, id.vars = c("Condition", "Population"), 
                variable_name = "Gene")
# Calcolo la media del valore di citofluorimetria per ogni tripletta
df2plot <- df_melt %>% 
    dplyr::group_by(Condition, Population, Gene) %>%
    dplyr::summarise(Value = mean(value)) %>%
    ungroup() %>%
# Per ogni gene calcolo un valore riscalato "normalizzato" rispetto 
# a tutti i valori ottenuti per quel gene
    group_by(Gene) %>%
    mutate(Scaledvalue = scale(Value))

# Creo il grafico
library(ggplot2)
library(viridis)
library(ggh4x)
# Ordino i geni, da quello con valore normalizzato più basso al più alto
gene_order <- df2plot %>% group_by(Gene) %>%
    summarise(mag = sum(Scaledvalue)) %>%
    arrange(mag) %>%
    dplyr::select(Gene) %>% unlist()
# Stampo i valori
ggplot(df2plot, aes(x = Gene, y = Condition)) + 
    geom_point(aes(color = Scaledvalue, size = Scaledvalue)) +
    scale_x_discrete(limits = gene_order) +
    scale_color_viridis(option = "inferno") +
    facet_nested(Population ~ ., scales = "free") + 
    theme(strip.text.y = element_text(angle = 0),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")
