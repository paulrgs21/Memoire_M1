# 3.2 Approche asymptotique -------------------------------------------------
#On utilise simplement le fait que -2*ln(LR) converge en distribution vers khi2.

#Arguments : LR_val : Le rapport de vraisemblance entre H0 et H1.
#n_states : le nombre d'états dans le modèle SMP

compute_asymptotic_pvalue <- function(LR_val, n_states, absorbing_state, dist_type) {
  # Calcul des degrés de liberté
  k <- ifelse(dist_type %in% c("gamma", "weibull"), 2, 1) #calcul du nombre de paramètres 
  if (is.null(absorbing_state)) { #le degré de liberté depend de la présence ou non d'un état absorbant.
    d <- n_states^2 - n_states - 1 + k*n_states*(n_states-1)
  } else {
    d <- (n_states-1)^2 - 2*(n_states-1) + k*(n_states-1)^2
  }
  
  # Calcul de la p-value
  test_stat <- -2*log(LR_val)
  pval <- 1 - pchisq(test_stat, df = d)
  return(pval)
}
