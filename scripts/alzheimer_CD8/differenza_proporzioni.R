# Librerie necessarie
library(dplyr)
library(tidyr)

# Funzione principale
# Funzione principale
analyze_cell_proportions <- function(data, n_permutations = 10000) {
    
    # Calcola totali per ciascun gruppo
    total_wt <- sum(data$WT)
    total_tg <- sum(data$Tg)
    
    # Calcola proporzioni e Odds Ratio osservati
    results <- data %>%
        mutate(
            Proportion_WT = WT / total_wt,
            Proportion_Tg = Tg / total_tg,
            Odds_WT = WT / (total_wt - WT),
            Odds_Tg = Tg / (total_tg - Tg),
            Odds_Ratio = Odds_Tg / Odds_WT  # Inverte l'OR: Tg come riferimento
        )
    
    # Funzione per calcolare IC per l'OR osservato
    calc_or_ci <- function(wt_yes, tg_yes, total_wt, total_tg) {
        wt_no <- total_wt - wt_yes
        tg_no <- total_tg - tg_yes
        
        # Calcola OR osservato
        or_obs <- (tg_yes / tg_no) / (wt_yes / wt_no)
        
        # Calcola SE del log(OR)
        se_log_or <- sqrt(1 / wt_yes + 1 / wt_no + 1 / tg_yes + 1 / tg_no)
        
        # Calcola IC per l'OR
        ci_low <- exp(log(or_obs) - 1.96 * se_log_or)
        ci_high <- exp(log(or_obs) + 1.96 * se_log_or)
        
        list(or_obs = or_obs, ci_low = ci_low, ci_high = ci_high)
    }
    
    # Funzione per il test permutazionale
    perm_test <- function(observed_wt, observed_tg, total_wt, total_tg, n_permutations) {
        observed_odds <- (observed_tg / (total_tg - observed_tg)) / (observed_wt / (total_wt - observed_wt))
        
        # Permutazione casuale
        perm_odds <- replicate(n_permutations, {
            perm_tg <- rbinom(1, size = observed_wt + observed_tg, prob = total_tg / (total_wt + total_tg))
            perm_wt <- observed_wt + observed_tg - perm_tg
            perm_odds <- (perm_tg / (total_tg - perm_tg)) / (perm_wt / (total_wt - perm_wt))
            perm_odds
        })
        
        # Calcolo p-value
        p_value <- mean(abs(perm_odds - 1) >= abs(observed_odds - 1))
        list(p_value = p_value)
    }
    
    # Applica il calcolo per ogni tipo cellulare
    results <- results %>%
        rowwise() %>%
        mutate(
            CI = list(calc_or_ci(WT, Tg, total_wt, total_tg)),
            OR_CI_Low = CI$ci_low,
            OR_CI_High = CI$ci_high,
            Permutation_Test = list(perm_test(WT, Tg, total_wt, total_tg, n_permutations)),
            p_value = Permutation_Test$p_value
        ) %>%
        ungroup() %>%
        select(-CI, -Permutation_Test)
    
    # Correzione per confronti multipli (FDR)
    results <- results %>%
        mutate(
            Adjusted_p_value = p.adjust(p_value, method = "fdr")
        )
    
    return(results)
}

# Esegui l'analisi
# Specifica il path al file di input e output
input_file <- "./data/alzheimer/cd8/Number of cells.csv"

data <- read.csv(input_file, row.names = 1)

data_brain <- data[-nrow(data),] %>% 
    select(Brain_3xTg, Brain_WT) %>%
    rename("Tg" = "Brain_3xTg", 
           "WT" = "Brain_WT")

data_meninges <- data[-nrow(data),] %>% 
    select(Meninges_3xTg, Meninges_WT) %>%
    rename("Tg" = "Meninges_3xTg", 
           "WT" = "Meninges_WT")

data_all <- data[-nrow(data),] %>% 
    mutate(Tg = Brain_3xTg + Meninges_3xTg, 
           WT = Brain_WT + Meninges_WT) %>%
    select(Tg, WT)

results_brain <- analyze_cell_proportions(data_brain)
rownames(results_brain) <- rownames(data_brain)

results_meninges <- analyze_cell_proportions(data_meninges)
rownames(results_meninges) <- rownames(data_meninges)

results_all <- analyze_cell_proportions(data_all)
rownames(results_all) <- rownames(data_all)

write.csv(results_meninges, 
    "./data/alzheimer/cd8/meninges_results.csv", row.names = TRUE)
write.csv(results_brain, 
    "./data/alzheimer/cd8/brain_results.csv", row.names = TRUE)
write.csv(results_all, 
    "./data/alzheimer/cd8/all_results.csv", row.names = TRUE)

