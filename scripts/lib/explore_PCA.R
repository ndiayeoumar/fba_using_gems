# PCA instead of UMAP

```{r} 
library(FactoMineR)
library(factoextra)


res.pca <- stats::prcomp(t(compass_data$reaction_consistencies))
explained_variance_ratio <- 100 * summary(res.pca)[["importance"]]['Proportion of Variance',]
PCA_coord <- data.frame(t(t(res.pca[["x"]])/res.pca[["sdev"]]))
PCA_loads <- mapply(
    data.frame(res.pca[["rotation"]]), 
    summary(res.pca)[["sdev"]], 
    FUN = function(comp, variance) {
        loadings <- comp * variance
        return(loadings)
    }, SIMPLIFY = TRUE, USE.NAMES = TRUE)


plot(res.pca)
biplot(res.pca)
# Results for individuals
res.ind <- get_pca_ind(res.pca)
PCA_coord <- data.frame(res.ind$coord)[, 1:2] # Coordinates
PCA_coord$cell_id <- rownames(PCA_coord)

# We perform the UMAP and we add metadata informations
cell_info_with_PCA_components <- 
    cell_info_with_umap_components %>%
    left_join(
        PCA_coord,
        by = "cell_id"
    )
```

```{r rPCAbystage, fig.cap="Each dot represents a single cell coloured based on their stage (Chronic or Onset)."}
p_m_stage <- ggplot(cell_info_with_PCA_components, 
                    aes(x = Dim.1, y = Dim.2, color = stage)) +
    geom_point(alpha = 0.8) +
    theme(legend.position = "bottom") +
    # ggtitle("PCA colored by stage (Chronic or Onset)") +
    xlab(label = "PCA 1") +
    ylab(label = "PCA 2") +
    theme_void() +
    theme(legend.position = "bottom")
p_m_stage

t_p <- ggplot(cell_info_with_PCA_components, 
              aes(x = Dim.1, fill = stage)) +
    geom_density(alpha = 0.5) +
    theme_void() +
    theme(legend.position = "none")

r_p <- ggplot(cell_info_with_PCA_components, 
              aes(x = Dim.2, fill = stage)) +
    geom_density(alpha = 0.5) +
    coord_flip() +
    theme_void() +
    theme(legend.position = "none")

cowplot::plot_grid(plotlist = list(t_p, NULL, p_m_stage, r_p),
                   nrow = 2, rel_widths = c(5,2), rel_heights = c(2,5), 
                   align = "hv", axis = "lrtb")
```

```{r rPCAbycelltype, fig.cap="Each dot represents a single cell coloured based on their cell type."}
p_m_ct <- ggplot(cell_info_with_PCA_components,
                 aes(x = Dim.1, y = Dim.2, color = pruned_fine)) +
    geom_point(alpha = 0.8) +
    ggtitle("PCA colored by cell type") +
    xlab(label = "PCA 1") +
    ylab(label = "PCA 2") +
    scale_color_manual(values = RColorBrewer::brewer.pal(12, "Set3")) + 
    theme(legend.position = "bottom") +
    guides(color = guide_legend(title = "Cell type", ncol = 2))
p_m_ct
```

```{r rPCAbymetabolicactivity, fig.cap="Each dot represents a single cell coloured based on their metabolic activity."}
# umap according to metabolic activity
p_m_ma <- ggplot(cell_info_with_PCA_components,
                 aes(x = Dim.1, y = Dim.2, color = metabolic_activity)) +
    geom_point(alpha = 0.8) +
    theme(legend.position = "bottom") +
    ggtitle("PCA colored by metabolic activity") +
    xlab(label ="PCA 1") +
    ylab(label = "PCA 2") +
    scale_color_viridis_c() + 
    theme(legend.position = "bottom")

p_m_ma
```

```{r geneUMAPmetabolic, fig.cap="Each dot represents a single cell coloured by cluster membership computed from the reaction consistencies (a), stage (b), cell type (c), and metabolic activity (d). The UMAPs are based on the gene expression values.", fig.width=14, fig.height=16}
cowplot::plot_grid(plotlist = list(p_m_stage, p_m_ct, p_m_ma, NULL), nrow = 2, align = "hv", axis = "tblr", labels = "auto")
```