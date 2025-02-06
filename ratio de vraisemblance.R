# trajectories1 : Premier ensemble de trajectoires.
# trajectories2 : Deuxième ensemble de trajectoires.
# absorbing_state : L'état absorbant de la chaîne SMP.
# dist_type : Le type de distribution pour les durées


compute_LR <- function(trajectories1, trajectories2, absorbing_state, dist_type) {
  # Estimation sous H0 (données regroupées)
  mle_H0 <- estimate_SMP_params(c(trajectories1, trajectories2), 
                                absorbing_state, dist_type)
  
  # Estimation sous H1 (données séparées)
  mle_H1_t1 <- estimate_SMP_params(trajectories1, absorbing_state, dist_type)
  mle_H1_t2 <- estimate_SMP_params(trajectories2, absorbing_state, dist_type)
  
  # Calcul des log-vraisemblances
  logL_H0 <- log_likelihood(mle_H0, c(trajectories1, trajectories2),absorbing_state, dist_type)
  
  logL_H1 <- log_likelihood(mle_H1_t1, trajectories1, absorbing_state, dist_type) + log_likelihood(mle_H1_t2, trajectories2, absorbing_state, dist_type) #La log-vraisemblance totale est la somme des log-vraisemblances des deux ensembles
  
  LR <- exp(logL_H0 - logL_H1)
  return(list(LR = LR, log_LR = log(LR)))
}