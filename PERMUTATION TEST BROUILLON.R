data <- read.table("C:/Users/Valentin/OneDrive/Documents/donnees.txt", 
                   header = TRUE,      # La première ligne contient les noms de colonnes
                   sep = "\t",         # Délimiteur tabulation
                   quote = "",         # Pas de guillemets pour délimiter les chaînes
                   na.strings = "",    # Chaînes vides comme NA
                   fill = TRUE,        # Compléter les lignes trop courtes
                   stringsAsFactors = FALSE)

library(MASS)

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
  # Estimation alpha (fréquences initiales)
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


#TEST DE PERMUTATION
permutation_test <- function(base1, base2, n1, n2, R, n_states, dist_type) {
  
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
  
  T_l <- compute_LR(trajectoires1, trajectoires2, n_states, dist_type)
  
  count <- 0  
  
  for (r in 1:R) {
    
    permutation <- sample(1:(n1+n2), n1+n2, replace = FALSE) # Permutation aléatoire
    
    perm_traj1 <- c(trajectoires1, trajectoires2)[permutation[1:n1]]
    perm_traj2 <- c(trajectoires1, trajectoires2)[permutation[(n1+1):(n1+n2)]]
    
    T_star <- compute_LR(perm_traj1, perm_traj2, n_states, dist_type)
    
    if (T_star[[2]] <= T_l[[2]]) {
      count <- count + 1
    }
  }
  
  p_perm <- count / R
  
  return(p_perm)
}


R <- 50
n_states <- 9
dist_type <- "weibull"

hommes <- data[data$Q1==1,]
femmes <- data[data$Q1==2,]
base1 <- hommes
base2 <- femmes
n1 <- nrow(base1)
n2 <- nrow(base2)

permutation_test(base1, base2, n1, n2, R, n_states, dist_type)

