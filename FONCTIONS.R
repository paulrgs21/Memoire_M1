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
  return(list(alpha=alpha, P=round(P,digits=2), dist_params=round(dist_params,digits=2)))
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
  if (!is.finite(logL_H0) || !is.finite(logL_H1)) { #Evite message d'erreur pour trop de permutations
    return(list(LR = NA, log_LR = NA))
  }
  LR <- exp(logL_H0 - logL_H1)
  return(list(LR = LR, log_LR = log(LR), mle_H0 = mle_H0,
              mle_H1_t1 = mle_H1_t1, mle_H1_t2 = mle_H1_t2))
}



### Formatage SMP ----
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



### 1) Chi-2 ----
chi2 <- function(LR_val, n_states, dist_type) {
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




### 2) Permutation ----
permutation <- function(likelihood_ratio, R, n_states, dist_type) {
  
  count <- foreach(r = 1:R,
                   .combine = "+",
                   .packages = c("MASS"),
                   .export = c("formatage_smp", "compute_LR", "estimate_SMP_params", "log_likelihood")) %dopar% {
                     
                     permutation <- sample(1:(n1 + n2), n1 + n2, replace = FALSE)
                     
                     perm_traj1 <- c(trajectoires1, trajectoires2)[permutation[1:n1]]
                     perm_traj2 <- c(trajectoires1, trajectoires2)[permutation[(n1 + 1):(n1 + n2)]]
                     
                     T_star <- compute_LR(perm_traj1, perm_traj2, n_states, dist_type)
                     
                     if (!is.na(T_star$log_LR) && is.finite(T_star$log_LR) &&
                         !is.na(likelihood_ratio$log_LR) && is.finite(likelihood_ratio$log_LR) &&
                         T_star$log_LR <= likelihood_ratio$log_LR) {
                       return(1)
                     } else {
                       return(0)
                     }
                   }
  
  p_perm <- count / R
  
  return(p_perm)
}




### 3) Bootstrap ----






#### Fonction de test global ----
test <- function(base1, base2, n_states, dist_type){
  
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
  
  return(list(likelihood_ratio = likelihood_ratio))
}



