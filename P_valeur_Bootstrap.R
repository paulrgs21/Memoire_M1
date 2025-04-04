rm(list=ls())
library(MASS)
library(ggplot2)

### Simulation de trajectoires SMP 
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




### Test sur données réelles ----

#### Data management ----
data <- read.table("C:/Users/rapha/Desktop/Mémoire/donnees.txt", 
                   header = TRUE,      # La première ligne contient les noms de colonnes
                   sep = "\t",         # Délimiteur tabulation
                   quote = "",         # Pas de guillemets pour délimiter les chaînes
                   na.strings = "",    # Chaînes vides comme NA
                   fill = TRUE,        # Compléter les lignes trop courtes
                   stringsAsFactors = FALSE)



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



#### Fonction de test global ----
test_boostrap <- function(base1, base2, n_states, dist_type){
  
  n1 <- nrow(base1)
  n2 <- nrow(base2)
  n <- n1+n2
  
  trajectoires1 <- list()
  for(i in 1:nrow(base1)) {
    id <- base1$IDENT[i]
    trajectoires1[[as.character(id)]] <- formatage_smp(base1[i,])
  }
  
  trajectoires2 <- list()
  for(i in 1:nrow(base2)) {
    id <- base2$IDENT[i]
    trajectoires2[[as.character(id)]] <- formatage_smp(base2[i,])
  }
  
  likelihood_ratio <- compute_LR(trajectoires1, trajectoires2, n_states, dist_type)
  
  
  return(likelihood_ratio = likelihood_ratio)
}

library(skimr)
skim(data)
#### 2ème test ----
# tentative de test significatif :
# test entre individus similaires, n'étant différents que par le fait d'avoir
# eu une mobilité de commune durant le parcours scolaire ou non (Q31A)
# paramètres
n_states <- 9
dist_type <- "weibull"

base_1 <- data[data$Q1==1&data$Q31==11&data$perefr==1&data$merefr==1&data$Q53==3&data$Q52==3&data$Q31A==1,]
base_2 <- data[data$Q1==1&data$Q31==11&data$perefr==1&data$merefr==1&data$Q53==3&data$Q52==3&data$Q31A==2,]
trajectoires_smp1 <- list()
for(i in 1:nrow(data)) {
  id <- base_1$IDENT[i]
  trajectoires_smp1[[as.character(id)]] <- formatage_smp(data[i,])
}
trajectoires_smp2 <- list()
for(i in 1:nrow(data)) {
  id <- base_2$IDENT[i]
  trajectoires_smp2[[as.character(id)]] <- formatage_smp(data[i,])
}



test_boostrap(base_1, base_2, n_states, dist_type)

p_value_2 <- parametric_bootstrap(trajectoires_smp1, trajectoires_smp2, 1400, 330, n_states, 9, dist_type, 1000)







  
