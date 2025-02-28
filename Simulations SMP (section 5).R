# Code effectuant la section 5.1 et 5.2.1 de l'article
# donc on est que sous H_0 pour l'instant
# pas d'état absorbant
# paramètres des lois aléatoires
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



### Estimation theta ----
estimate_SMP_params <- function(trajectories, n_states, dist_type) {
  # Estimation α (fréquences initiales)
  initial_states <- sapply(trajectories, function(x) x$states[1]) #Extrait le premier état de chaque trajectoire pour identifier les états initiaux
  alpha <- table(factor(initial_states, levels=1:n_states))
  alpha <- alpha/sum(alpha) # proportion de chaque état initial
  
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
  dist_params <- list()
  for (state in 1:n_states) {
    durations <- unlist(lapply(trajectories, function(traj) {
      idx <- which(traj$states[-length(traj$states)] == state) #identifie les indices (idx) où l'état courant apparaît avant la dernière transition (on exclut l'état absorbant)
      if(length(idx) > 0) return(traj$times[idx]) 
      else return(NULL)
    }))
    
    if(dist_type == "gamma") {
      fit <- fitdistr(durations, "gamma", start=list(shape=1, rate=1))
      dist_params[[state]] <- c(fit$estimate["shape"], fit$estimate["rate"])
    } else if(dist_type == "weibull") {
      fit <- fitdistr(durations, "weibull", start=list(shape=1, scale=1))
      dist_params[[state]] <- c(fit$estimate["shape"], fit$estimate["scale"])
    } else {
      fit <- fitdistr(durations, "exponential", start=list(shape=1))
      dist_params[[state]] <- c(fit$estimate["rate"])
    }
  }
  return(list(alpha=alpha, P=P, dist_params=dist_params))
}



### Log-vraisemblance ----
log_likelihood <- function(params, trajectories, dist_type) {
  alpha <- params$alpha
  P <- params$P
  dist_params <- params$dist_params
  logL <- 0 # Initialisation de la log-vraisemblance
  epsilon <- 1e-10 # Petit terme pour éviter log(0)
  
  for (traj in trajectories) {
    states <- traj$states # Ordre des états
    times <- traj$times # Durées passées dans chaque état
    
    # Terme initial alpha_j1
    logL <- logL + log(max(alpha[states[1]], epsilon)) # Évite log(0)
    
    m <- length(states) - 1 # Nombre de transitions
    
    for (i in 1:m) { # Boucle sur les transitions entre états
      from <- states[i]
      to <- states[i+1]
      
      # Terme de transition P_ij
      logL <- logL + log(max(P[from, to], epsilon)) # Ajoute log proba de transition
      
      # Terme de durée
      if (dist_type == "gamma") {
        logL <- logL + dgamma(times[i], shape = dist_params[[from]][1], rate = dist_params[[from]][2], log = TRUE)
      } else if (dist_type == "weibull") {
        logL <- logL + dweibull(times[i], shape = dist_params[[from]][1], scale = dist_params[[from]][2], log = TRUE)
      } else if (dist_type == "exponential") {
        logL <- logL + dexp(times[i], rate = dist_params[[from]][1], log = TRUE)
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
        times[t] <- rgamma(1, shape = dist_params[[states[t]]][1],
                           rate = dist_params[[states[t]]][2])
      } else if (dist_type == "weibull") {
        times[t] <- rweibull(1, shape = dist_params[[states[t]]][1],
                             scale = dist_params[[states[t]]][2])
      } else {
        times[t] <- rexp(1, rate = dist_params[[states[t]]][1])
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
R <- 1000

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


# chi-2

n_repetitions <- 500
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

line_df <- data.frame(x = c(0, 1), y = c(0, 1))
ggplot(p_values_df, aes(x = p_value)) +
  stat_ecdf(geom = "step", color = "blue", size = 1) +  # Fonction de répartition empirique
  geom_line(data = line_df, aes(x = x, y = y), color = "red", size = 1) +  # Droite y = x
  labs(title = "Fonction de répartition empirique des p-valeurs",
       x = "P-valeurs",
       y = "F(x)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))



# parametric bootstrap
p_value_2 <- parametric_bootstrap(smp_trajectories1, smp_trajectories2, n1, n2, n_states, max_transitions, dist_type, R)

# permutation
p_value_3 <- permutation_test(smp_trajectories1, smp_trajectories2, n1, n2, R)

# très (trop) long à tourner : à optimiser pour reproduire Fig. 4