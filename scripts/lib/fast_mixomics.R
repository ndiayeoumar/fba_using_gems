library(mixOmics)

fast_mixomics <- function(counts, metadata, Y_name, sample_name, multilevel = NULL){
  # In Y avrò i nomi dei miei gruppi (Diverse etichette)
  Y = unlist(metadata[, Y_name])
  # Casto i gruppi come fattori
  Y = as.factor(Y)

  # Faccio la trasposta
  psT = t(counts)

  pca = mixOmics::pca(psT, ncomp = 10)

  # PCA
  #plot(pca, main = paste0("Principal components - ",Y_name))
  plot_pca <- plotIndiv(pca,
            ellipse = TRUE,
            comp = c(1,2), # the components to plot
            ind.names = NULL,
            group = Y,
            legend = TRUE,
            title = paste0("PCA - ",Y_name))

  # PLS-DA

  # Vado a calcolare un modello PLS-DA da usare per stimare il tasso di errore per ogni componente
  plsda = plsda(X = psT, Y, ncomp = nlevels(Y), logratio = 'none', multilevel = multilevel)
  perf = perf(plsda, validation = 'loo', progressBar = FALSE)

  # Tabella errore Overall
  error_TOT_pslda <- perf
  
  plot_plsda <- plotIndiv(plsda, 
            abline = TRUE, centroid = TRUE,
            style = "ggplot2",
            comp = c(1:2),
            ind.names = NULL,
            ellipse = TRUE,
            legend = TRUE,
            title = paste0("PLS-DA - ",Y_name))

  # sPLS-DA
  set.seed(33)

  tune = tune.splsda(psT, Y, ncomp = 2,
                     logratio = 'none',
                     test.keepX = c(seq(10, 100, 10), seq(150, 500, 50)),
                     validation = 'loo', dist = 'max.dist', progressBar = FALSE)

  # plot(tune)
  # optimal number of variables to select on 2 comps:
  keepx = tune$choice.keepX
  print(keepx)

  splsda = splsda(psT, Y, ncomp = 2,
                  logratio = 'none',
                  keepX = keepx, multilevel = multilevel)

  plot_splsda <- plotIndiv(splsda, comp = c(1,2), abline = TRUE, centroid = TRUE,
            style = "ggplot2",
            ind.names = NULL,
            ellipse = TRUE, legend = TRUE,
            title = paste0("sPLS-DA - ",Y_name))

  # Tasso di errore della SPLS-DA
  set.seed(34)  # for reproducible results for this code
  perf.splsda = perf(splsda, validation = 'loo', progressBar = FALSE)
  # Tabella errore Overall
  max.dist.O <- data.frame(Max_dist. = c((perf.splsda$error.rate$overall)[,1]))
  cent.dist.O <- data.frame(Centroid_dist. = c((perf.splsda$error.rate$overall)[,2]))
  mahala.dist.O <- data.frame(Mahalanobis_dist. = c((perf.splsda$error.rate$overall)[,3]))
  error_TOT_O <- data.frame(Max_dist. = max.dist.O, Centroid_dist. = cent.dist.O, Mahalanobis_dist. = mahala.dist.O)

  # Error-rate della sPLS-DA
  error <- error_TOT_O

  # ASV selezionate dalla sPLS-DA

  #Nomi delle ASV più discriminanti
  ASVs_comp1 = selectVar(splsda, comp = 1)$name
  ASVs_comp2 = selectVar(splsda, comp = 2)$name

  # Comp 1 (valore tra le [[]])
  comp1tab = data.frame(perf.splsda$features$stable[[1]][ASVs_comp1])
  colnames(comp1tab)[1] = "Reaction"
  colnames(comp1tab)[2] = "Frequency"
  # Ordino la tabella per le frequenze in valore discendente
  comp1tab = comp1tab[order(-comp1tab$Frequency),]

  # Comp 2
  comp2tab = data.frame(perf.splsda$features$stable[[2]][ASVs_comp2])
  colnames(comp2tab)[1] = "Reaction"
  colnames(comp2tab)[2] = "Frequency"
  # Ordino la tabella per le frequenze in valore discendente
  comp2tab = comp2tab[order(-comp2tab$Frequency),]

  results <- list("pca" = plot_pca,
                  "plsda" = plot_plsda,
                  "plsda_data" = plsda,
                  "splsda" = plot_splsda,
                  "splsda_data" = splsda,
                  "ASVs_1" = comp1tab,
                  "ASVs_2" = comp2tab,
                  "error_plsda" = error_TOT_pslda,
                  "error" = error,
                  "nfeatures" = keepx)
return (results)
}
