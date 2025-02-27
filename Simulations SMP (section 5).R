simulate_SMP <- function(n, n_states, P, alpha, dist_type, dist_params, max_transitions) {
  set.seed(123)
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
      
      # Arrêt si état absorbant (ex : dernier état)
      #if (states[t + 1] == n_states) {
      #  states <- states[1:(t + 1)]
      #  times <- times[1:t]
      #  break
      #}
    }
    trajectories[[i]] <- list(states = states, times = times)
  }
  return(trajectories)
}


### Exemple d'utilisation
n <- 100  # Nombre de trajectoires simulées
n_states <- 4  # Nombre d'états
max_transitions <- 10
dist_type <- "gamma"

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
  dist_params <- matrix(c(1.75, 2.95,
                          3.86, 2.19,
                          4.33, 4.26,
                          4.47, 3.61), ncol = 2, byrow = TRUE)
} else if (dist_type == "weibull") {
  dist_params <- matrix(c(1.30, 3.19,
                          1.80, 3.90,
                          2.05, 4.85,
                          4.04, 1.62), ncol = 2, byrow = TRUE)
} else {
  dist_params <- matrix(c(1.11, 0.83, 1.34, 1.86), ncol = 1)
}

smp_trajectories <- simulate_SMP(n, n_states, P, alpha, dist_type, dist_params, max_transitions)


# nombre de ddl sous H0 que suit la stat de test (on considère l'hyp w_lj = w_l)
# pour l'instant, pas d'état absorbant
d <- n_states^2 + n_states - 1