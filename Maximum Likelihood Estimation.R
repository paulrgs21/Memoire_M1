

estimate_SMP_params <- function(trajectories, absorbing_state, dist_type, n_states) {
  # Estimation α (fréquences initiales)
  initial_states <- sapply(trajectories, function(x) x$states[1]) #Extrait le premier état de chaque trajectoire pour identifier les états initiaux
  alpha <- table(factor(initial_states, levels=1:n_states))
  alpha <- as.vector(alpha/sum(alpha)) # proportion de chaque état initial
  
  # Estimation P (matrice de transitions)
  transition_counts <- matrix(0, n_states, n_states) # Création matrice n*n
  for (traj in trajectories) {
    states <- traj$states
    for (i in 1:(length(states)-1)) {
      transition_counts[states[i], states[i+1]] <- transition_counts[states[i], states[i+1]] + 1
    }
  }
  P <- transition_counts/rowSums(transition_counts)
  
  # Estimation ω (MLE par état)
  omega <- list()
  for (state in 1:n_states) {
    durations <- unlist(lapply(trajectories, function(traj) {
      idx <- which(traj$states[-length(traj$states)] == state) #identifie les indices (idx) où l'état courant apparaît avant la dernière transition (on exclut l'état absorbant)
      if(length(idx) > 0) return(traj$times[idx]) 
      else return(NULL)
    }))
    
    # Vérification : si durations est vide, on met des valeurs par défaut
    if (length(durations) == 0) {
      omega[[state]] <- c(NA, NA)  # Met des valeurs NA au lieu d'estimer
      next  # Passe à l'état suivant
    }
    
    if(dist_type == "gamma") {
      fit <- fitdistr(durations, "gamma", start=list(shape=1, rate=1))
      omega[[state]] <- c(fit$estimate["shape"], fit$estimate["rate"])
    } else if(dist_type == "weibull") {
      fit <- fitdistr(durations, "weibull", start=list(shape=1, scale=1))
      omega[[state]] <- c(fit$estimate["shape"], fit$estimate["scale"]) #La fonction fitdistr ajuste les paramètres de la distribution en maximisant la vraisemblance des données de durées observées.
    } else {
      fit <- fitdistr(durations, "exponential", start=list(shape=1))
      omega[[state]] <- c(fit$estimate["rate"])
    }
  }
  return(list(alpha=alpha, P=P, omega=omega))
}