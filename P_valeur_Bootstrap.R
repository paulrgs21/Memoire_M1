rm(list=ls())
library(MASS)
library(ggplot2)
library(foreach)
library(doParallel)

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
  # Estimation fréquences initiales
  initial_states <- sapply(trajectories, function(x) x$states[1])
  alpha <- table(factor(initial_states, levels=1:n_states))
  alpha <- alpha/sum(alpha)
  
  # Estimation matrice de transitions
  #on extrait sous forme de vecteur chaque changement de trajectoire
  all_from <- unlist(lapply(trajectories, function(traj) traj$states[-length(traj$states)]))
  all_to   <- unlist(lapply(trajectories, function(traj) traj$states[-1]))
  
  M_from <- diag(n_states)[all_from, , drop = FALSE]
  M_to   <- diag(n_states)[all_to, , drop = FALSE]  
  
  #on utilise le produit matriciel
  transition_counts <- t(M_from) %*% M_to             
  
  P <- transition_counts / rowSums(transition_counts)
  
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
    if (length(durations) > 1) {
      if(dist_type == "gamma") {
        fit <- fitdistr(durations, "gamma", start=list(shape=1, rate=1))
        dist_params[state, ] <- c(fit$estimate["shape"], fit$estimate["rate"])
      } else if(dist_type == "weibull") {
        fit <- fitdistr(durations, "weibull", start=list(shape=1, scale=1))
        dist_params[state, ] <- c(fit$estimate["shape"], fit$estimate["scale"])
      } else {
        fit <- fitdistr(durations, "exponential", start=list(rate=1))
        dist_params[state, ] <- fit$estimate["rate"]
      }
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
    
    logL <- logL + log(max(alpha[states[1]], epsilon))
    
    from <- states[-length(states)]
    to   <- states[-1]
    durations <- times[-length(times)]
    
    log_transi <- log(pmax(P[cbind(from, to)], epsilon))
    log_transi[!is.finite(log_transi)] <- log(epsilon)
    logL <- logL + sum(log_transi)
    
    valid <- !is.na(from) & from >= 1 & from <= nrow(dist_params)
    
    if (any(valid)) { #plein de sécurités pour que le code puisse tourner
      if (dist_type == "gamma") {
        dens <- dgamma(durations, shape = dist_params[from, 1], rate = dist_params[from, 2])
        log_dens <- log(pmax(dens, epsilon))
        log_dens[!is.finite(log_dens)] <- log(epsilon)
        logL <- logL + sum(log_dens)
        
      } else if (dist_type == "weibull") {
        dens <- dweibull(durations, shape = dist_params[from, 1], scale = dist_params[from, 2])
        log_dens <- log(pmax(dens, epsilon))
        log_dens[!is.finite(log_dens)] <- log(epsilon)
        logL <- logL + sum(log_dens)
        
      } else if (dist_type == "exponential") {
        dens <- dexp(durations, rate = dist_params[from, 1])
        log_dens <- log(pmax(dens, epsilon))
        log_dens[!is.finite(log_dens)] <- log(epsilon)
        logL <- logL + sum(log_dens)
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

parametric_bootstrap <- function(trajectories1, trajectories2, 
                                 n1, n2, n_states, max_transitions, 
                                 dist_type, R, n_cores = parallel::detectCores() - 1) {
  
  # 📌 Calcul des paramètres du MLE sous H0
  mle_params <- estimate_SMP_params(c(trajectories1, trajectories2), n_states, dist_type)
  
  # 📌 Statistique observée
  T_l <- compute_LR(trajectories1, trajectories2, n_states, dist_type)
  
  # 📌 Lancement du cluster
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  # Export des fonctions + objets nécessaires dans les workers
  parallel::clusterExport(cl, varlist = c(
    "simulate_SMP_bootstrap", 
    "estimate_SMP_params", 
    "compute_LR",
    "log_likelihood",
    "n1", "n2", "n_states", "max_transitions", "dist_type", "mle_params", "T_l"
  ), envir = environment())  # important d'inclure l'environnement
  
  # Sécurité : si une erreur arrive dans le foreach, on stoppe le cluster proprement
  results <- tryCatch({
    
    foreach::foreach(r = 1:R, .combine = '+', .packages = c("stats","MASS")) %dopar% {
      bootstrap_trajectories1 <- simulate_SMP_bootstrap(n1, n_states, mle_params, dist_type, max_transitions)
      bootstrap_trajectories2 <- simulate_SMP_bootstrap(n2, n_states, mle_params, dist_type, max_transitions)
      
      T_star <- compute_LR(bootstrap_trajectories1, bootstrap_trajectories2, n_states, dist_type)
      
      as.integer(T_star[[2]] <= T_l[[2]])
    }
    
  }, finally = {
    # Stoppe le cluster même s’il y a une erreur
    parallel::stopCluster(cl)
  })
  
  # 🔢 Calcul final de la p-valeur
  p_boot <- results / R
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







  
