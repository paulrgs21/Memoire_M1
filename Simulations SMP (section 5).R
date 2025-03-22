# Code effectuant la section 5.1 et 5.2.1 de l'article
# donc on est que sous H_0 pour l'instant
# pas d'état absorbant
# paramètres des lois aléatoires pour simus
# hypothèse w_lj = w_l

library(MASS)
library(ggplot2)

### Simulation de trajectoires SMP ----
simulate_SMP <- function(n, n_states, P, alpha, dist_type, dist_params, max_transitions) {
  trajectories <- vector("list", n)
  
  for (i in 1:n) {
    states <- numeric(max_transitions)
    times <- numeric(max_transitions)
    states[1] <- sample(1:n_states, 1, prob = alpha) # État initial
    
    for (t in 1:(max_transitions - 1)) {
      states[t + 1] <- sample(1:n_states, 1, prob = P[states[t], ]) # Transition
      
      # Tirage du temps de séjour selon la distribution spécifiée
      if (dist_type == "gamma") {
        times[t] <- rgamma(1, shape = dist_params[states[t], 1],
                           rate = dist_params[states[t], 2])
      } else if (dist_type == "weibull") {
        times[t] <- rweibull(1, shape = dist_params[states[t], 1],
                             scale = dist_params[states[t], 2])
      } else {
        times[t] <- rexp(1, rate = dist_params[states[t], 1])
      }
    
    }
    trajectories[[i]] <- list(states = states, times = times)
  }
  return(trajectories)
}



### Estimation theta ----
estimate_SMP_params <- function(trajectories, n_states, dist_type) {
  # Estimation α (fréquences initiales)
  initial_states <- sapply(trajectories, function(x) x$states[1])
  alpha <- table(factor(initial_states, levels=1:n_states))
  alpha <- alpha/sum(alpha)
  
  # Estimation P (matrice de transitions)
  transition_counts <- matrix(0, n_states, n_states)
  for (traj in trajectories) {
    states <- traj$states
    for (i in 1:(length(states)-1)) {
      transition_counts[states[i], states[i+1]] <- transition_counts[states[i], states[i+1]] + 1
    }
  }
  P <- transition_counts/rowSums(transition_counts)
  
  # Estimation ω (MLE par état)
  if(dist_type == "gamma" || dist_type == "weibull") {
    dist_params <- matrix(0, nrow = n_states, ncol = 2)
  } else {
    dist_params <- matrix(0, nrow = n_states, ncol = 1)
  }
  
  for (state in 1:n_states) {
    durations <- unlist(lapply(trajectories, function(traj) {
      idx <- which(traj$states[-length(traj$states)] == state)
      if(length(idx) > 0) return(traj$times[idx]) 
      else return(NULL)
    }))
    
    if(dist_type == "gamma") {
      fit <- fitdistr(durations, "gamma", start=list(shape=1, rate=1))
      dist_params[state, ] <- c(fit$estimate["shape"], fit$estimate["rate"])
    } else if(dist_type == "weibull") {
      fit <- fitdistr(durations, "weibull", start=list(shape=1, scale=1))
      dist_params[state, ] <- c(fit$estimate["shape"], fit$estimate["scale"])
    } else {
      fit <- fitdistr(durations, "exponential", start=list(shape=1))
      dist_params[state, ] <- fit$estimate["rate"]
    }
  }
  return(list(alpha=alpha, P=P, dist_params=dist_params))
}


### Log-vraisemblance ----
log_likelihood <- function(params, trajectories, dist_type) {
  alpha <- params$alpha
  P <- params$P
  dist_params <- params$dist_params
  logL <- 0
  epsilon <- 1e-10
  
  for (traj in trajectories) {
    states <- traj$states
    times <- traj$times
    
    # Terme initial alpha_j1
    logL <- logL + log(max(alpha[states[1]], epsilon))
    
    m <- length(states) - 1
    
    for (i in 1:m) {
      from <- states[i]
      to <- states[i+1]
      
      # Terme de transition P_ij
      logL <- logL + log(max(P[from, to], epsilon))
      
      # Terme de durée
      if (dist_type == "gamma") {
        logL <- logL + dgamma(times[i], shape = dist_params[from, 1], rate = dist_params[from, 2], log = TRUE)
      } else if (dist_type == "weibull") {
        logL <- logL + dweibull(times[i], shape = dist_params[from, 1], scale = dist_params[from, 2], log = TRUE)
      } else if (dist_type == "exponential") {
        logL <- logL + dexp(times[i], rate = dist_params[from, 1], log = TRUE)
      }
    }
  }
  
  return(logL)
}



### Ratio de vraisemblance ----
compute_LR <- function(trajectories1, trajectories2, n_states, dist_type) {
  # Estimation sous H0 (données regroupées)
  mle_H0 <- estimate_SMP_params(c(trajectories1, trajectories2), n_states, dist_type)
  
  # Estimation sous H1 (données séparées)
  mle_H1_t1 <- estimate_SMP_params(trajectories1, n_states, dist_type)
  mle_H1_t2 <- estimate_SMP_params(trajectories2, n_states, dist_type)
  
  # Calcul des log-vraisemblances
  logL_H0 <- log_likelihood(mle_H0, c(trajectories1, trajectories2), dist_type)
  
  logL_H1 <- log_likelihood(mle_H1_t1, trajectories1, dist_type) + log_likelihood(mle_H1_t2, trajectories2, dist_type)
  
  LR <- exp(logL_H0 - logL_H1)
  return(list(LR = LR, log_LR = log(LR)))
}



### 1) Chi-2 ----
compute_asymptotic_pvalue <- function(LR_val, n_states, dist_type) {
  # Calcul des degrés de liberté
  k <- ifelse(dist_type %in% c("gamma", "weibull"), 2, 1) #calcul du nombre de paramètres 
  #if (is.null(absorbing_state)) { #le degré de liberté depend de la présence ou non d'un état absorbant.
  d <- n_states^2 - n_states - 1 + k*n_states # car hyp w_lj=w_l
  #} else { pour l'instant, pas d'état absorbant
    #d <- n_states^2 - 2*n_states +k*(n_states - 1)
  #}
  
  # Calcul de la p-value
  test_stat <- -2*LR_val[[2]]
  pval <- 1 - pchisq(test_stat, df = d)
  return(pval)
}



### 2) Parametric bootstrap ----
# on vient générer un échantillon de trajectoires suivant les paramètres du MLE sous H0 :
simulate_SMP_bootstrap <- function(n, n_states, params, dist_type, max_transitions) {
  alpha <- params$alpha
  P <- params$P
  dist_params <- params$dist_params
  
  trajectories <- vector("list", n)
  
  for (i in 1:n) {
    states <- numeric(max_transitions)
    times <- numeric(max_transitions)
    states[1] <- sample(1:n_states, 1, prob = alpha) # État initial
    
    for (t in 1:(max_transitions - 1)) {
      states[t + 1] <- sample(1:n_states, 1, prob = P[states[t], ]) # Transition
      
      # Tirage du temps de séjour selon la distribution spécifiée
      if (dist_type == "gamma") {
        times[t] <- rgamma(1, shape = dist_params[states[t], 1],
                           rate = dist_params[states[t], 2])
      } else if (dist_type == "weibull") {
        times[t] <- rweibull(1, shape = dist_params[states[t], 1],
                             scale = dist_params[states[t], 2])
      } else {
        times[t] <- rexp(1, rate = dist_params[states[t], 1])
      }
      
    }
    trajectories[[i]] <- list(states = states, times = times)
  }
  return(trajectories)
}



#On vient ensuite effectuer R fois ce processus en calculant la statistique de test à chaque fois

parametric_bootstrap <- function(trajectories1, trajectories2, n1, n2, n_states, max_transitions, dist_type, R) {
  
  # On calcule les paramètres du MLE (pour générer les échantillons)
  mle_params <- estimate_SMP_params(c(trajectories1, trajectories2), n_states, dist_type)
  
  # Calcul de la statistique de test T_l
  T_l <- compute_LR(trajectories1, trajectories2, n_states, dist_type)
  
  # Compteur pour la p-value
  count <- 0
  for (r in 1:R) {
    
    # Génére n trajectoires indépendantes sous H0 avec les paramètres du MLE
    bootstrap_trajectories1 <- simulate_SMP_bootstrap(n1, n_states, mle_params, dist_type, max_transitions)
    bootstrap_trajectories2 <- simulate_SMP_bootstrap(n2, n_states, mle_params, dist_type, max_transitions)
    
    # Calcule la stat de test T* pour les trajectoires générées
    T_star <- compute_LR(bootstrap_trajectories1, bootstrap_trajectories2, n_states, dist_type)
    
    # Itérations pour le calcul de la p-value
    if (T_star[[2]] <= T_l[[2]]) {
      count <- count + 1
    }
  }
  
  p_boot <- count / R
  
  return(p_boot)
}



### 3) Permutation ----
permutation_test <- function(trajectories1, trajectories2, n1, n2, R) {
  
  T_l <- compute_LR(trajectories1, trajectories2, n_states, dist_type)
  
  count <- 0  
  
  for (r in 1:R) {
    
    permutation <- sample(1:(n1+n2), n1+n2, replace = FALSE) # Permutation aléatoire
    
    perm_traj1 <- c(trajectories1, trajectories2)[permutation[1:n1]]
    perm_traj2 <- c(trajectories1, trajectories2)[permutation[(n1+1):(n1+n2)]]
    
    T_star <- compute_LR(perm_traj1, perm_traj2, n_states, dist_type)
    
    if (T_star[[2]] <= T_l[[2]]) {
      count <- count + 1
    }
  }
  
  p_perm <- count / R
  
  return(p_perm)
}



### Exemple d'utilisation ----
n1 <- 60  # Nombre de trajectoires simulées
n2 <- 60
n <- n1+n2
n_states <- 7  # Nombre d'états
max_transitions <- 5
dist_type <- "gamma"
niveau_test <- 0.05
R <- 1000
n_repetitions <- 500 # pour graphique

# Génération de P et alpha
# pour l'instant, pas d'état absorbant
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


# Echantillon 1 (100 trajectoires)
smp_trajectories1 <- simulate_SMP(n1, n_states, P, alpha, dist_type, dist_params, max_transitions)
# Echantillon 2 (100 trajectoires, on est sous H_0 donc même setup que échantillon 1)
smp_trajectories2 <- simulate_SMP(n2, n_states, P, alpha, dist_type, dist_params, max_transitions)

likelihood_ratio <- compute_LR(smp_trajectories1, smp_trajectories2, n_states, dist_type)


### chi-2

# Initialisation du data frame pour stocker les p-valeurs
p_values_df <- data.frame(p_value = numeric(n_repetitions))
# Boucle pour répéter l'expérience 500 fois
for (i in 1:n_repetitions) {
  # Génération des échantillons sous H0
  smp_trajectories1 <- simulate_SMP(n1, n_states, P, alpha, dist_type, dist_params, max_transitions)
  smp_trajectories2 <- simulate_SMP(n2, n_states, P, alpha, dist_type, dist_params, max_transitions)
  
  # Calcul du ratio de vraisemblance
  likelihood_ratio <- compute_LR(smp_trajectories1, smp_trajectories2, n_states, dist_type)
  
  # Calcul de la p-valeur
  p_values_df$p_value[i] <- compute_asymptotic_pvalue(likelihood_ratio, n_states, dist_type)
}

# Graphique à la Fig. 4
line_df <- data.frame(x = c(0, 1), y = c(0, 1))
ggplot(p_values_df, aes(x = p_value)) +
  stat_ecdf(geom = "step", color = "blue", size = 1) +  # Fonction de répartition empirique
  geom_line(data = line_df, aes(x = x, y = y), color = "red", size = 1) +  # Droite y = x
  labs(title = "Fonction de répartition empirique des p-valeurs",
       x = "P-valeurs",
       y = "F(x)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# niveau empirique du test à la Table 1
compute_empirical_level <- function(p_values_df, niveau_test) {
  # Calcul du niveau empirique : proportion de p-values < niveau_test
  empirical_level <- mean(p_values_df$p_value < niveau_test)*100
  
  return(empirical_level)
}

empirical_level <- compute_empirical_level(p_values_df, niveau_test)


### parametric bootstrap
p_value_2 <- parametric_bootstrap(smp_trajectories1, smp_trajectories2, n1, n2, n_states, max_transitions, dist_type, R)

### permutation
p_value_3 <- permutation_test(smp_trajectories1, smp_trajectories2, n1, n2, R)

# très (trop) long à tourner : à optimiser pour reproduire Fig. 4




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
# n_states = 9
# Fonction pour transformer les données d'un individu en format semi-markovien
formatage_smp <- function(row) {
  # Extraire les colonnes c1 à c94 qui représentent les états mensuels
  colonnes_etats <- paste0("c", 1:94)
  states <- as.numeric(row[colonnes_etats])
  
  # Identifier les changements d'état (TRUE quand il y a changement)
  changements <- c(TRUE, states[-1] != states[-length(states)])
  
  # Extraire la séquence d'états uniques (sans répétition d'états consécutifs identiques)
  etats_uniques <- states[changements]
  
  # Calculer les temps de séjour dans chaque état
  times <- as.double(rle(states)$lengths)
  
  # Retourner une liste avec les états et les durées
  return(list(
    states = etats_uniques,
    times = times
  ))
}

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
hommes <- data[data$Q1==1,]
femmes <- data[data$Q1==2,]

trajectoires_h <- list()
for(i in 1:nrow(hommes)) {
  id <- hommes$IDENT[i]
  trajectoires_h[[as.character(id)]] <- formatage_smp(hommes[i,])
}

trajectoires_f <- list()
for(i in 1:nrow(femmes)) {
  id <- femmes$IDENT[i]
  trajectoires_f[[as.character(id)]] <- formatage_smp(femmes[i,])
}

# paramètres
n1 <- nrow(hommes)
n2 <- nrow(femmes)
n <- n1+n2
n_states <- 9
dist_type <- "weibull"
niveau_test <- 0.05
R <- 1000
n_repetitions <- 500 # pour graphique
max_transitions <- 5 # à estimer (pour bootstrap ou graphique à la Fig.4)

likelihood_ratio <- compute_LR(trajectoires_h, trajectoires_f, n_states, dist_type)

# p-valeur du Chi-2
compute_asymptotic_pvalue(likelihood_ratio, n_states, dist_type)
