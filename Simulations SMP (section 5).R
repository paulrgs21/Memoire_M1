# Code effectuant la section 5.1 et 5.2.1 de l'article
# donc on est que sous H_0 pour l'instant
# pas d'état absorbant
# paramètres des lois aléatoires pour simus
# hypothèse w_lj = w_l

library(ggplot2)

### Simus théoriques ----
#### Paramétrage ----
n1 <- 30  # Nombre de trajectoires simulées
n2 <- 30
n <- n1+n2
n_states <- 4  # Nombre d'états
max_transitions <- 5
dist_type <- "weibull"
niveau_test <- 0.01
R <- 500
n_repetitions <- 100 # pour graphique

# Génération de P et alpha
# pas d'état absorbant
P <- matrix(runif(n_states^2), nrow = n_states)
diag(P) <- 0
P <- P / rowSums(P)

alpha <- runif(n_states,0,1)
alpha <- alpha/sum(alpha)

# Paramètres Gamma/Weibull/Exp pour chaque état (on considère l'hyp w_lj = w_l)
# pris aléatoirement pour le moment
if (dist_type == "gamma") {
  dist_params <- matrix(round(runif(n_states * 2, min = 1, max = 5), 2), ncol = 2, byrow = TRUE)
} else if (dist_type == "weibull") {
  dist_params <- matrix(round(runif(n_states * 2, min = 1, max = 5), 2), ncol = 2, byrow = TRUE)
} else {
  dist_params <- matrix(round(runif(n_states, min = 0.5, max = 2), 2), ncol = 1)
}


#### Graphiques ----
### chi-2
chi2_df <- data.frame(p_value = numeric(n_repetitions))
for (i in 1:n_repetitions) {
  smp_trajectories1 <- simulate_SMP(n1, n_states, P, alpha, dist_type, dist_params, max_transitions)
  smp_trajectories2 <- simulate_SMP(n2, n_states, P, alpha, dist_type, dist_params, max_transitions)
  likelihood_ratio <- compute_LR(smp_trajectories1, smp_trajectories2, n_states, dist_type)
  chi2_df$p_value[i] <- chi2(likelihood_ratio, n_states, dist_type)
}


### permutation
permutation_df <- data.frame(p_value = numeric(n_repetitions))

nb_cores <- parallel::detectCores() - 1
cl <- makeCluster(nb_cores)
registerDoParallel(cl)

for (i in 1:n_repetitions) {
  trajectories1 <- simulate_SMP(n1, n_states, P, alpha, dist_type, dist_params, max_transitions)
  trajectories2 <- simulate_SMP(n2, n_states, P, alpha, dist_type, dist_params, max_transitions)
  likelihood_ratio <- compute_LR(trajectories1, trajectories2, n_states, dist_type)
  
  permutation_df$p_value[i] <- permutation(likelihood_ratio, R, n_states, n1, n2,
                                           trajectories1, trajectories2, dist_type)
}

stopCluster(cl)

# Graphique à la Fig. 4
line_df <- data.frame(x = c(0, 1), y = c(0, 1))
chi2_df$source <- "Test du Chi2"
permutation_df$source <- "Test de permutation"
combined_df <- rbind(chi2_df, permutation_df)
ggplot() +
  stat_ecdf(data = chi2_df, aes(x = p_value, color = "Test du Chi2"), 
            geom = "step", size = 1) +
  stat_ecdf(data = permutation_df, aes(x = p_value, color = "Test de permutation"), 
            geom = "step", size = 1) +
  geom_line(data = line_df, aes(x = x, y = y), color = "red", size = 1, linetype = "dashed") +
  labs(title = "Fonctions de répartition empiriques des p-valeurs",
       x = "P-valeurs",
       y = "F(x)",
       color = "Méthode") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom") +
  scale_color_manual(values = c("Test du Chi2" = "blue", "Test de permutation" = "darkgreen"))




#### Niveaux empiriques ----
# niveau empirique du test à la Table 1
compute_empirical_level <- function(p_values_df, niveau_test) {
  # Calcul du niveau empirique : proportion de p-values < niveau_test
  empirical_level <- mean(p_values_df$p_value < niveau_test)*100
  
  return(empirical_level)
}


chi2_level <- compute_empirical_level(chi2_df, niveau_test)
# 8
permutation_level <- compute_empirical_level(permutation_df, niveau_test)
# 4




#### Puissances ----
# SOUS H1
P1 <- matrix(runif(n_states^2), nrow = n_states)
diag(P1) <- 0
P1 <- P1 / rowSums(P1)

P2 <- matrix(runif(n_states^2), nrow = n_states)
diag(P2) <- 0
P2 <- P2 / rowSums(P2)

alpha <- runif(n_states,0,1)
alpha <- alpha/sum(alpha)

# Paramètres Gamma/Weibull/Exp pour chaque état (on considère l'hyp w_lj = w_l)
# pris aléatoirement pour le moment
if (dist_type == "gamma") {
  dist_params1 <- matrix(round(runif(n_states * 2, min = 1, max = 5), 2), ncol = 2, byrow = TRUE)
  dist_params2 <- matrix(round(runif(n_states * 2, min = 1, max = 5), 2), ncol = 2, byrow = TRUE)
} else if (dist_type == "weibull") {
  dist_params1 <- matrix(round(runif(n_states * 2, min = 1, max = 5), 2), ncol = 2, byrow = TRUE)
  dist_params2 <- matrix(round(runif(n_states * 2, min = 1, max = 5), 2), ncol = 2, byrow = TRUE)
} else {
  dist_params1 <- matrix(round(runif(n_states, min = 0.5, max = 2), 2), ncol = 1)
  dist_params2 <- matrix(round(runif(n_states, min = 0.5, max = 2), 2), ncol = 1)
}

### chi-2
chi2_df <- data.frame(p_value = numeric(n_repetitions))
for (i in 1:n_repetitions) {
  smp_trajectories1 <- simulate_SMP(n1, n_states, P1, alpha, dist_type, dist_params1, max_transitions)
  smp_trajectories2 <- simulate_SMP(n2, n_states, P2, alpha, dist_type, dist_params2, max_transitions)
  likelihood_ratio <- compute_LR(smp_trajectories1, smp_trajectories2, n_states, dist_type)
  chi2_df$p_value[i] <- chi2(likelihood_ratio, n_states, dist_type)
}

mean(chi2_df$p_value < niveau_test)


### permutation
permutation_df <- data.frame(p_value = numeric(n_repetitions))

nb_cores <- parallel::detectCores() - 1
cl <- makeCluster(nb_cores)
registerDoParallel(cl)

for (i in 1:n_repetitions) {
  trajectories1 <- simulate_SMP(n1, n_states, P1, alpha, dist_type, dist_params1, max_transitions)
  trajectories2 <- simulate_SMP(n2, n_states, P2, alpha, dist_type, dist_params2, max_transitions)
  likelihood_ratio <- compute_LR(trajectories1, trajectories2, n_states, dist_type)
  
  permutation_df$p_value[i] <- permutation(likelihood_ratio, R, n_states, n1, n2,
                                           trajectories1, trajectories2, dist_type)
}

stopCluster(cl)

mean(permutation_df$p_value < niveau_test)



### Test sur données réelles ----

#### Data management ----
data <- read.table("C:/Users/paulr/Documents/M1 Dauphine/S1/mémoire/données/donnees.txt", 
                   header = TRUE,      # La première ligne contient les noms de colonnes
                   sep = "\t",         # Délimiteur tabulation
                   quote = "",         # Pas de guillemets pour délimiter les chaînes
                   na.strings = "",    # Chaînes vides comme NA
                   fill = TRUE,        # Compléter les lignes trop courtes
                   stringsAsFactors = FALSE)

# à faire :
# 1) gestion des NA
# 2) formatage des trajectoires pour correspondre à nos fonctions
# 3) test existence état absorbant
# 4) identification des paramètres et de la loi des temps de séjour
# 5) nombre de transitions maximal à estimer



## 1)
library(skimr)
skim(data)
# 0 NA pour les variables qui nous intéressent pour le moment, ie genre (Q1) et c1-94



## 2)
# Appliquer la transformation à chaque individu et créer la liste finale
trajectoires_smp <- list()
for(i in 1:nrow(data)) {
  id <- data$IDENT[i]
  trajectoires_smp[[as.character(id)]] <- formatage_smp(data[i,])
}



## 3) test si état 1 est absorbant
# Fonction pour tester si l'état 1 est absorbant
test_etat_absorbant <- function(trajectoires_smp) {
  # Initialiser les compteurs
  n_individus_total <- length(trajectoires_smp)
  n_individus_avec_etat_1 <- 0
  n_individus_absorbants <- 0
  n_individus_non_absorbants <- 0
  
  # Liste pour stocker les cas non-absorbants (pour inspection détaillée)
  cas_non_absorbants <- list()
  
  # Parcourir tous les individus
  for (id in names(trajectoires_smp)) {
    states <- trajectoires_smp[[id]]$states
    
    # Vérifier si l'individu a l'état 1 dans sa trajectoire
    if (1 %in% states) {
      n_individus_avec_etat_1 <- n_individus_avec_etat_1 + 1
      
      # Trouver la première occurrence de l'état 1
      premiere_pos_etat_1 <- which(states == 1)[1]
      
      # Vérifier si tous les états après le premier état 1 sont aussi 1
      if (premiere_pos_etat_1 < length(states)) {
        etats_apres_premier_1 <- states[(premiere_pos_etat_1 + 1):length(states)]
        
        if (all(etats_apres_premier_1 == 1)) {
          # Cas absorbant
          n_individus_absorbants <- n_individus_absorbants + 1
        } else {
          # Cas non absorbant
          n_individus_non_absorbants <- n_individus_non_absorbants + 1
          
          # Stocker ce cas pour analyse plus détaillée
          cas_non_absorbants[[id]] <- list(
            states = states,
            premiere_pos_etat_1 = premiere_pos_etat_1,
            etats_apres = etats_apres_premier_1
          )
        }
      } else {
        # Si l'état 1 est le dernier état, on le considère comme absorbant
        n_individus_absorbants <- n_individus_absorbants + 1
      }
    }
  }
  
  # Calculer le pourcentage d'absorption
  pourcentage_absorption <- (n_individus_absorbants / n_individus_avec_etat_1) * 100
  
  # Résumé des résultats
  resultats <- list(
    n_individus_total = n_individus_total,
    n_individus_avec_etat_1 = n_individus_avec_etat_1,
    n_individus_absorbants = n_individus_absorbants,
    n_individus_non_absorbants = n_individus_non_absorbants,
    pourcentage_absorption = pourcentage_absorption,
    cas_non_absorbants = cas_non_absorbants
  )
  
  # Afficher un résumé
  cat("Résultats du test d'absorption pour l'état 1:\n")
  cat("Nombre total d'individus:", n_individus_total, "\n")
  cat("Nombre d'individus passant par l'état 1:", n_individus_avec_etat_1, "\n")
  cat("Nombre d'individus pour lesquels l'état 1 est absorbant:", n_individus_absorbants, "\n")
  cat("Nombre d'individus pour lesquels l'état 1 n'est PAS absorbant:", n_individus_non_absorbants, "\n")
  cat("Pourcentage d'absorption:", sprintf("%.2f%%", pourcentage_absorption), "\n")
  
  if (n_individus_non_absorbants > 0) {
    cat("\nExemples de cas où l'état 1 n'est pas absorbant (jusqu'à 5 cas):\n")
    id_exemples <- names(cas_non_absorbants)[1:min(5, length(cas_non_absorbants))]
    for (id in id_exemples) {
      cat("ID:", id, "\n")
      cat("  Séquence d'états:", cas_non_absorbants[[id]]$states, "\n")
      cat("  Position du premier état 1:", cas_non_absorbants[[id]]$premiere_pos_etat_1, "\n\n")
    }
  }
  
  return(resultats)
}

resultats_test <- test_etat_absorbant(trajectoires_smp)



## 4) temps de séjour : échelle et loi
# Fonction pour extraire les temps de séjour par état
extraire_temps_sejour <- function(trajectoires_smp) {
  # Dictionnaire pour stocker les temps de séjour par état
  temps_sejour_par_etat <- list()
  
  # Parcourir toutes les trajectoires
  for (id in names(trajectoires_smp)) {
    states <- trajectoires_smp[[id]]$states
    times <- trajectoires_smp[[id]]$times
    
    # Associer chaque durée avec son état correspondant
    for (i in 1:length(states)) {
      state <- states[i]
      time <- times[i]
      
      # Initialiser la liste si nécessaire
      if (is.null(temps_sejour_par_etat[[as.character(state)]])) {
        temps_sejour_par_etat[[as.character(state)]] <- numeric(0)
      }
      
      # Ajouter la durée à la liste des durées pour cet état
      temps_sejour_par_etat[[as.character(state)]] <- c(temps_sejour_par_etat[[as.character(state)]], time)
    }
  }
  
  return(temps_sejour_par_etat)
}

# Fonction pour ajuster les distributions avec MASS::fitdistr
test_distrib <- function(temps_sejour_par_etat) {
  resultats <- list()
  
  for (state in names(temps_sejour_par_etat)) {
    # Skip if no data for this state
    if (length(temps_sejour_par_etat[[state]]) == 0) {
      next
    }
    
    times <- temps_sejour_par_etat[[state]]
    
    cat("\n--- Analyse de l'état", state, "---\n")
    cat("Nombre d'observations:", length(times), "\n")
    cat("Durée moyenne:", mean(times), "\n")
    cat("Durée min:", min(times), "\n")
    cat("Durée max:", max(times), "\n")
    
    # Ajuster les distributions avec MASS::fitdistr
    tryCatch({
      # Nombre d'observations pour le calcul du BIC
      n <- length(times)
      
      # Distribution exponentielle
      fit_exp <- fitdistr(times, "exponential")
      # Distribution gamma
      fit_gamma <- fitdistr(times, "gamma")
      # Distribution Weibull
      fit_weibull <- fitdistr(times, "weibull")
      
      # Calculer les valeurs de log-vraisemblance 
      loglik_exp <- fit_exp$loglik
      loglik_gamma <- fit_gamma$loglik
      loglik_weibull <- fit_weibull$loglik
      
      # Calculer AIC: -2*logLik + 2*k (k = nombre de paramètres)
      aic_exp <- -2 * loglik_exp + 2 * 1  # Exponentielle a 1 paramètre
      aic_gamma <- -2 * loglik_gamma + 2 * 2  # Gamma a 2 paramètres
      aic_weibull <- -2 * loglik_weibull + 2 * 2  # Weibull a 2 paramètres
      
      # Calculer BIC: -2*logLik + k*ln(n) (k = nombre de paramètres, n = nombre d'observations)
      bic_exp <- -2 * loglik_exp + 1 * log(n)  # Exponentielle a 1 paramètre
      bic_gamma <- -2 * loglik_gamma + 2 * log(n)  # Gamma a 2 paramètres
      bic_weibull <- -2 * loglik_weibull + 2 * log(n)  # Weibull a 2 paramètres
      
      # Afficher les résultats
      cat("\nParamètres estimés:\n")
      cat("Exponentielle: rate =", fit_exp$estimate, "\n")
      cat("Gamma: shape =", fit_gamma$estimate["shape"], ", rate =", fit_gamma$estimate["rate"], "\n")
      cat("Weibull: shape =", fit_weibull$estimate["shape"], ", scale =", fit_weibull$estimate["scale"], "\n")
      
      cat("\nCritères d'ajustement (AIC - plus petit = meilleur):\n")
      aic_df <- data.frame(
        Distribution = c("Exponentielle", "Gamma", "Weibull"),
        AIC = c(aic_exp, aic_gamma, aic_weibull)
      )
      print(aic_df)
      
      cat("\nCritères d'ajustement (BIC - plus petit = meilleur):\n")
      bic_df <- data.frame(
        Distribution = c("Exponentielle", "Gamma", "Weibull"),
        BIC = c(bic_exp, bic_gamma, bic_weibull)
      )
      print(bic_df)
      
      # Trouver la meilleure distribution selon AIC
      best_aic_idx <- which.min(c(aic_exp, aic_gamma, aic_weibull))
      best_dist_aic <- c("Exponentielle", "Gamma", "Weibull")[best_aic_idx]
      cat("\nMeilleure distribution selon AIC:", best_dist_aic, "\n")
      
      # Trouver la meilleure distribution selon BIC
      best_bic_idx <- which.min(c(bic_exp, bic_gamma, bic_weibull))
      best_dist_bic <- c("Exponentielle", "Gamma", "Weibull")[best_bic_idx]
      cat("Meilleure distribution selon BIC:", best_dist_bic, "\n")
      
      # Générer un histogramme des données avec les courbes ajustées
      x_range <- seq(min(times), max(times), length.out = 100)
      
      # Densités théoriques
      dens_exp <- dexp(x_range, rate = fit_exp$estimate)
      dens_gamma <- dgamma(x_range, shape = fit_gamma$estimate["shape"], rate = fit_gamma$estimate["rate"])
      dens_weibull <- dweibull(x_range, shape = fit_weibull$estimate["shape"], scale = fit_weibull$estimate["scale"])
      
      # Créer un dataframe pour ggplot (si besoin de visualisation)
      plot_data <- data.frame(
        x = x_range,
        Exponentielle = dens_exp,
        Gamma = dens_gamma,
        Weibull = dens_weibull
      )
      
      # Enregistrer les résultats
      resultats[[state]] <- list(
        exp = fit_exp,
        gamma = fit_gamma,
        weibull = fit_weibull,
        aic = aic_df,
        bic = bic_df,
        best_dist_aic = best_dist_aic,
        best_dist_bic = best_dist_bic
      )
      
    }, error = function(e) {
      cat("Erreur lors de l'ajustement pour l'état", state, ":", e$message, "\n")
    })
  }
  
  return(resultats)
}



temps_sejour_par_etat <- extraire_temps_sejour(trajectoires_smp)

# Afficher des statistiques de base pour tous les états
for (state in names(temps_sejour_par_etat)) {
  times <- temps_sejour_par_etat[[state]]
  if (length(times) > 0) {
    cat("État", state, ": n =", length(times), ", moyenne =", mean(times), 
        ", médiane =", median(times), ", max =", max(times), "\n")
  }
}

# Tester l'ajustement des distributions pour tous les états
resultats_ajustement <- test_distrib(temps_sejour_par_etat)

# Résumé des meilleures distributions pour chaque état
cat("État\tAIC\tBIC\n")
for (state in names(resultats_ajustement)) {
  cat(state, "\t", resultats_ajustement[[state]]$best_dist_aic, 
      "\t", resultats_ajustement[[state]]$best_dist_bic, "\n")
}



# Graphiques de comparaison avec lois théoriques
visualiser_ajustement <- function(state, temps_sejour_par_etat, resultats_ajustement) {
  if (!(state %in% names(resultats_ajustement))) {
    cat("Pas de résultats disponibles pour l'état", state, "\n")
    return(NULL)
  }
  
  times <- temps_sejour_par_etat[[state]]
  fit_exp <- resultats_ajustement[[state]]$exp
  fit_gamma <- resultats_ajustement[[state]]$gamma
  fit_weibull <- resultats_ajustement[[state]]$weibull
  
  x_range <- seq(min(times), max(times), length.out = 100)
  
  # Calculer les densités théoriques
  dens_exp <- dexp(x_range, rate = fit_exp$estimate) * length(times)
  dens_gamma <- dgamma(x_range, shape = fit_gamma$estimate["shape"], 
                       rate = fit_gamma$estimate["rate"]) * length(times)
  dens_weibull <- dweibull(x_range, shape = fit_weibull$estimate["shape"], 
                           scale = fit_weibull$estimate["scale"]) * length(times)
  
  # Créer le graphique avec ggplot2
  hist_data <- data.frame(time = times)
  
  p <- ggplot(hist_data, aes(x = time)) +
    geom_histogram(aes(y = ..count..), binwidth = max(1, (max(times) - min(times))/30),
                   fill = "lightblue", color = "black", alpha = 0.7) +
    geom_line(data = data.frame(x = x_range, y = dens_exp), 
              aes(x = x, y = y, color = "Exponentielle"), size = 1) +
    geom_line(data = data.frame(x = x_range, y = dens_gamma), 
              aes(x = x, y = y, color = "Gamma"), size = 1) +
    geom_line(data = data.frame(x = x_range, y = dens_weibull), 
              aes(x = x, y = y, color = "Weibull"), size = 1) +
    scale_color_manual(name = "Distribution", 
                       values = c("Exponentielle" = "red", "Gamma" = "blue", "Weibull" = "green")) +
    labs(title = paste("Ajustement des distributions pour l'état", state),
         subtitle = paste("Meilleure distribution:", resultats_ajustement[[state]]$best_dist),
         x = "Durée de séjour", y = "Fréquence") +
    theme_minimal()
  
  print(p)
  return(p)
}

# Pour visualiser l'ajustement d'un état spécifique:
visualiser_ajustement("1", temps_sejour_par_etat, resultats_ajustement)
visualiser_ajustement("2", temps_sejour_par_etat, resultats_ajustement)
visualiser_ajustement("3", temps_sejour_par_etat, resultats_ajustement)
visualiser_ajustement("4", temps_sejour_par_etat, resultats_ajustement)
visualiser_ajustement("5", temps_sejour_par_etat, resultats_ajustement)
visualiser_ajustement("6", temps_sejour_par_etat, resultats_ajustement)
visualiser_ajustement("7", temps_sejour_par_etat, resultats_ajustement)
visualiser_ajustement("8", temps_sejour_par_etat, resultats_ajustement)
visualiser_ajustement("9", temps_sejour_par_etat, resultats_ajustement)

# On a 6 Weibull contre 3 Gamma, donc on choisit la loi Weibull



# 5) nombre de transitions maximal à estimer
nombre_transitions <- sapply(trajectoires_smp, function(trajectoire) {
  length(trajectoire$states) - 1
})
summary(nombre_transitions)




#### 1er test H vs F ----
# Leurs trajectoires professionnelles suivent-ils la même loi SMP ?

# paramètres
n_states <- 9
dist_type <- "weibull"

# Création des sous-bases
hommes <- data[data$Q1==1,]
femmes <- data[data$Q1==2,]

test(hommes, femmes, n_states, dist_type)
# LR de 0, p value de 0, on rejette tout le temps H_0
# différences notables dans les matrices de transition, par exemple au niveau du
# service militaire. Les alpha sont identiques (on débute en études). Les
# paramètres des lois des temps de séjour sont sensiblement similaires.

test(hommes, hommes, n_states, dist_type)
test(femmes, femmes, n_states, dist_type)
# LR de 1, p value de 1, on ne rejette JAMAIS H_0 (test trivial)




#### autres tests ----
# tentative de test significatif :
# test entre individus similaires, n'étant différents que par le fait d'avoir
# eu une mobilité de commune durant le parcours scolaire ou non (Q31A)

base_1 <- data[data$Q1==1&data$Q31==11&data$perefr==1&data$merefr==1&data$Q53==3&data$Q52==3&data$Q31A==1,]
base_2 <- data[data$Q1==1&data$Q31==11&data$perefr==1&data$merefr==1&data$Q53==3&data$Q52==3&data$Q31A==2,]

test2 <- test(base_1, base_2, n_states, dist_type)

chi2(test2$likelihood_ratio, n_states, dist_type)

nb_cores <- parallel::detectCores() - 1
cl <- makeCluster(nb_cores)
registerDoParallel(cl)

permutation(test2$likelihood_ratio, 1000, n_states, nrow(base_1), nrow(base_2),
            test2$trajectoires1, test2$trajectoires2, dist_type)

stopCluster(cl)
# on est sous H0, p valeur de 1 avec chi2, p valeur de 0.471 avec permut


base_1 <- data[data$Q1==1&data$Q31==11&data$perefr==1&data$merefr==1&data$Q53==3&data$Q52==3&data$nivdip7>4&data$Q31A==1,]
base_2 <- data[data$Q1==1&data$Q31==11&data$perefr==1&data$merefr==1&data$Q53==3&data$Q52==3&data$nivdip7>4&data$Q31A==2,]
test3 <- test(base_1, base_2, n_states, dist_type)

chi2(test3$likelihood_ratio, n_states, dist_type)

nb_cores <- parallel::detectCores() - 1
cl <- makeCluster(nb_cores)
registerDoParallel(cl)

permutation(test3$likelihood_ratio, 1000, n_states, nrow(base_1), nrow(base_2),
            test3$trajectoires1, test3$trajectoires2, dist_type)

stopCluster(cl)
# on est sous H0 ? p-valeur de 0.0005 avec chi-2, 0.155 puis 0.16 avec permut
# peut-être sous H1 car chi2 rejette trop souvent H0




### Simus avec param réels ----
### (setup : à changer (permutation marche pas ou prends
### trop de temps pour un n aussi grand) !!!)
### on prend le setup du dernier test
n1 <- 339
n2 <- 192
n <- n1+n2
n_states <- 9
max_transitions <- 5
dist_type <- "weibull"
R <- 1000

test3$likelihood_ratio



### SOUS H0

# Génération de P et alpha
P <- matrix(c(0.00, 0.20, 0.01, 0.04, 0.08, 0.42, 0.17, 0.08, 0.00,
              0.54, 0.00, 0.01, 0.01, 0.04, 0.26, 0.08, 0.06, 0.01,
              0.67, 0.05, 0.00, 0.00, 0.00, 0.29, 0.00, 0.00, 0.00,
              0.46, 0.19, 0.00, 0.00, 0.00, 0.31, 0.04, 0.00, 0.00,
              0.41, 0.19, 0.01, 0.00, 0.00, 0.26, 0.06, 0.08, 0.00,
              0.39, 0.31, 0.02, 0.05, 0.10, 0.00, 0.07, 0.06, 0.02,
              0.26, 0.22, 0.02, 0.01, 0.05, 0.18, 0.00, 0.23, 0.02,
              0.27, 0.19, 0.01, 0.00, 0.03, 0.34, 0.16, 0.00, 0.01,
              0.23, 0.17, 0.00, 0.01, 0.02, 0.17, 0.32, 0.08, 0.00), 
            nrow=9, ncol=9, byrow=TRUE)

alpha <- c(0,0,0,0,0,0,0,0,1)

# Paramètres Weibull (on considère l'hyp w_lj = w_l)
dist_params <- matrix(c(1.15, 22.86,
                        1.07, 10.69,
                        1.54, 17.30,
                        1.43, 36.10,
                        1.08, 10.18,
                        1.06, 7.14,
                        0.98, 4.33,
                        3.15, 12.15,
                        2.14, 8.23), 
                      nrow=9, ncol=2, byrow=TRUE)

# Echantillon 1
smp_trajectories1 <- simulate_SMP(n1, n_states, P, alpha, dist_type, dist_params, max_transitions)
# Echantillon 2 (on est sous H_0 donc même setup que échantillon 1)
smp_trajectories2 <- simulate_SMP(n2, n_states, P, alpha, dist_type, dist_params, max_transitions)

compute_LR(smp_trajectories1, smp_trajectories2, n_states, dist_type)
# paramètres sensiblement identiques aux vrais paramètres
# on les stocke pour utilisation :
likelihood_ratio <- compute_LR(smp_trajectories1, smp_trajectories2, n_states, dist_type)


# on n'utilise que permutation à partir de maintenant

nb_cores <- parallel::detectCores() - 1
cl <- makeCluster(nb_cores)
registerDoParallel(cl)

permutation(likelihood_ratio, R, n_states, n1, n2, smp_trajectories1, smp_trajectories2, dist_type)

stopCluster(cl)
# p-valeur de 0.7, on est bien sous H0



### SOUS H1

# Génération de P et alpha
P1 <- matrix(c(0.00, 0.20, 0.02, 0.04, 0.08, 0.42, 0.14, 0.10, 0.00,
               0.51, 0.00, 0.01, 0.01, 0.04, 0.26, 0.08, 0.07, 0.01,
               0.69, 0.06, 0.00, 0.00, 0.00, 0.25, 0.00, 0.00, 0.00,
               0.44, 0.22, 0.00, 0.00, 0.00, 0.28, 0.06, 0.00, 0.00,
               0.42, 0.20, 0.02, 0.00, 0.00, 0.27, 0.05, 0.05, 0.00,
               0.33, 0.31, 0.02, 0.05, 0.12, 0.00, 0.07, 0.08, 0.03,
               0.27, 0.22, 0.01, 0.02, 0.05, 0.20, 0.00, 0.22, 0.01,
               0.25, 0.19, 0.02, 0.00, 0.04, 0.33, 0.16, 0.00, 0.01,
               0.20, 0.17, 0.01, 0.01, 0.02, 0.18, 0.34, 0.08, 0.00), 
             nrow=9, ncol=9, byrow=TRUE)

P2 <- matrix(c(0.00, 0.20, 0.00, 0.05, 0.07, 0.43, 0.21, 0.04, 0.00,
               0.61, 0.00, 0.00, 0.00, 0.03, 0.25, 0.07, 0.03, 0.01,
               0.60, 0.00, 0.00, 0.00, 0.00, 0.40, 0.00, 0.00, 0.00,
               0.50, 0.12, 0.00, 0.00, 0.00, 0.38, 0.00, 0.00, 0.00,
               0.38, 0.17, 0.00, 0.00, 0.00, 0.21, 0.08, 0.17, 0.00,
               0.54, 0.29, 0.02, 0.03, 0.06, 0.00, 0.05, 0.01, 0.00,
               0.26, 0.22, 0.03, 0.01, 0.06, 0.14, 0.00, 0.25, 0.03,
               0.31, 0.20, 0.00, 0.00, 0.00, 0.35, 0.14, 0.00, 0.00,
               0.27, 0.17, 0.00, 0.01, 0.04, 0.14, 0.29, 0.09, 0.00), 
             nrow=9, ncol=9, byrow=TRUE)

alpha <- c(0,0,0,0,0,0,0,0,1)

# Paramètres Weibull (on considère l'hyp w_lj = w_l)
dist_params1 <- matrix(c(1.09, 22.11,
                         1.04, 10.07,
                         1.72, 17.98,
                         1.32, 30.68,
                         1.07, 11.00,
                         1.03, 7.30,
                         1.02, 4.38,
                         3.22, 12.23,
                         2.07, 8.18), 
                       nrow=9, ncol=2, byrow=TRUE)


dist_params2 <- matrix(c(1.30, 24.25,
                         1.16, 12.15,
                         1.18, 14.86,
                         2.02, 48.61,
                         1.16, 8.01,
                         1.17, 6.74,
                         0.92, 4.23,
                         2.97, 11.93,
                         2.40, 8.35), 
                       nrow=9, ncol=2, byrow=TRUE)


# Echantillon 1
smp_trajectories1 <- simulate_SMP(n1, n_states, P1, alpha, dist_type, dist_params1, max_transitions)
smp_trajectories2 <- simulate_SMP(n2, n_states, P2, alpha, dist_type, dist_params2, max_transitions)

compute_LR(smp_trajectories1, smp_trajectories2, n_states, dist_type)
# paramètres sensiblement identiques aux vrais paramètres
# on les stocke pour utilisation :
likelihood_ratio <- compute_LR(smp_trajectories1, smp_trajectories2, n_states, dist_type)

nb_cores <- parallel::detectCores() - 1
cl <- makeCluster(nb_cores)
registerDoParallel(cl)

permutation(likelihood_ratio, R, n_states, n1, n2, smp_trajectories1, smp_trajectories2, dist_type)

stopCluster(cl)
# p-valeur de 0.011, on rejette bien H0