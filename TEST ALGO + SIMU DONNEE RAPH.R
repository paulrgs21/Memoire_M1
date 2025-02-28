rm(list = ls())
library(fitdistrplus)
# Création manuelle des trajectoires pour 4 individus
# Trois état : 1=chomage, 2=cadre, 3=ouvrier
# Individu 1
# On va diviser en catégorie pour prévoir déjà le test homme vs femme
# l'individu 1 et 2 repésentent les hommes
# le 3 et 4 femme
# vu que le n est faible il faut que j'exagère le difference entre homme et femme
# pour l'homme il sera essentiellement cadre et la femme au chommage et un peu ouvrier.
individu1 <- list(
  states = c(2, 1, 2),  # États visités
  times = c(8, 1, 40)  # Temps de séjour dans chaque état en année
)

# Individu 2
individu2 <- list(
  states = c(1, 2, 1),  # États visités
  times = c(4, 40, 1)  # Temps de séjour dans chaque état
)

# Individu 3
individu3 <- list(
  states = c(1,3,1),  # États visités
  times = c(22,8,2)  # Temps de séjour dans chaque état
)

# Individu 4
individu4 <- list(
  states = c(3,1,3),  # États visités
  times = c(10, 11, 16)  # Temps de séjour dans chaque état
)

# Création de la liste des trajectoires (avec 4 indivius seulement l'optimisation plante...)
trajectories <- list(individu1, individu2, individu3, individu4,individu1, individu1, individu2, individu3, individu4,individu1, individu2, individu3, individu4,individu1, individu2, individu3, individu4,individu1, individu2, individu3, individu4,individu1, individu2, individu3, individu4,individu2, individu3, individu4,individu1, individu2, individu3, individu4,individu1, individu2, individu3, individu4,individu1, individu2, individu3, individu4,individu1, individu2, individu3, individu4,individu1, individu2, individu3, individu4)

dist_type = "gamma" #je sera à nous de la derteminer à partir des données ou de la fixer en avance comme je fais la ...

n_states=3




# Test de l'algo MLE :
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

params <- estimate_SMP_params(trajectories,0,dist_type,n_states) #sockage dans params pour l'algo de log vraissemblance
# Pour rappel : 
#Params : contient : alpha les proba initiales des états, P la matrice de transition entre les états, 
#et omega la liste des paramètres des lois de durée associées aux états


#Test de la fonction de Log_vraissemblance :
log_likelihood <- function(params, trajectories, absorbing_state, dist_type) {
  alpha <- params$alpha
  P <- params$P
  omega <- params$omega #On extrait les données depuis params
  logL <- 0 # Initialisation de la log vraisemblance
  
  for (traj in trajectories) { # On parcourt toutes les trajectoires
    states <- traj$states #ordre des états
    times <- traj$times #durées passées dans chaque état
    
    # Terme initial alpha_j1
    logL <- logL + log(alpha[states[1]]) #ajoute la proba initiale de commencer dans un état s1 est alpha_j1
    
    for (i in 1:(length(states)-1)) {#Boucle sur les transitions entre états
      from <- states[i] # Etat de départ
      to <- states[i+1] # Etat d'arrivée
      
      # Terme de transition P_ij
      logL <- logL + log(P[from, to]) #ajoute la proba de transition entre les états from et to
      
      # Terme de durée (sauf dernière transition absorbante)
      if (to != absorbing_state) {#Si to est un état absorbant on passe à la transition suivante
        if (dist_type == "gamma") {#Cas ou la loi des durées suit une distribution gamma.
          logL <- logL + dgamma(times[i],shape=omega[[from]][1],rate=omega[[from]][2], log=TRUE)#Shape = paramètre k, rate=parametre teta
        } 
        else if (dist_type == "weibull") {#Cas ou la loi des durées suit une distribution Weibull
          logL <- logL + dweibull(times[i],shape=omega[[from]][1],scale=omega[[from]][2], log=TRUE)
        }
        else if (dist_type == "exponential") {
          logL <- logL + dexp(times[i], rate = omega[[from]][1], log = TRUE)
        }
      }
    }
  }
  return(logL)
}


log <- log_likelihood(params, trajectories, 0, dist_type)



#Maintenant test du ratio de vraissemblance
#il faut donc que je créer deux jeux de trajectoires séparé :
#l'optimisation échoue à chaque fois, je pense que c'est du à l'effectif, je vais l'augmenter.
#Et le diversifier, quand je réplique les meme individu, ca ne marche pas, ca fais des matrice de transition avec des 0 du coup ca fait erreur
individu1BIS <- list(
  states = c(1, 2, 1),  # États visités
  times = c(8, 4, 28)  # Temps de séjour dans chaque état en année
)

# Individu 2
individu2BIS <- list(
  states = c(2, 1, 2),  # États visités
  times = c(10, 3, 50)  # Temps de séjour dans chaque état
)

# Individu 3 
individu3BIS <- list(
  states = c(3,1,3),  # États visités
  times = c(2,4,10)  # Temps de séjour dans chaque état
)

# Individu 4
individu4BIS <- list(
  states = c(1,3,1),  # États visités
  times = c(10, 1, 36)  # Temps de séjour dans chaque état
)


trajectories2 <- list(individu3,individu4,individu3BIS,individu4BIS,individu3,individu4,individu3BIS,individu4BIS,individu3,individu4,individu3BIS,individu4BIS,individu3,individu4,individu3BIS,individu4BIS,individu3,individu4,individu3BIS,individu4BIS,individu3,individu4,individu3BIS,individu4BIS,individu3,individu4,individu3BIS,individu4BIS,individu3,individu4,individu3BIS,individu4BIS,individu3,individu4,individu3BIS,individu4BIS,individu3,individu4,individu3BIS,individu4BIS,individu3,individu4,individu3BIS,individu4BIS,individu3,individu4,individu3BIS,individu4BIS,individu3,individu4,individu3BIS,individu4BIS)
trajectories1 <- list(individu1,individu2,individu1BIS,individu2BIS)

compute_LR <- function(trajectories1, trajectories2, absorbing_state, dist_type, n_states) {
  # Estimation sous H0 (données regroupées)
  mle_H0 <- estimate_SMP_params(c(trajectories1, trajectories2), 
                                absorbing_state, dist_type, n_states)
  
  # Estimation sous H1 (données séparées)
  mle_H1_t1 <- estimate_SMP_params(trajectories1, absorbing_state, dist_type, n_states)
  mle_H1_t2 <- estimate_SMP_params(trajectories2, absorbing_state, dist_type, n_states)
  
  # Calcul des log-vraisemblances
  logL_H0 <- log_likelihood(mle_H0, c(trajectories1, trajectories2),absorbing_state, dist_type)
  
  logL_H1 <- log_likelihood(mle_H1_t1, trajectories1, absorbing_state, dist_type) + log_likelihood(mle_H1_t2, trajectories2, absorbing_state, dist_type) #La log-vraisemblance totale est la somme des log-vraisemblances des deux ensembles
  
  LR <- exp(logL_H0 - logL_H1)
  return(list(LR = LR, log_LR = log(LR)))
}
LR <-compute_LR(trajectories1, trajectories2, 0, dist_type, n_states)
print(LR)
LR_val <- LR$log_LR

print(estimate_SMP_params(trajectories2,0,dist_type,n_states))

print(estimate_SMP_params(trajectories1, 0, dist_type, n_states))

print(estimate_SMP_params(trajectories, 0, dist_type, n_states))



# Bien maintenant on va test les p_value avec les 2 trajectoires très diffèrentes pour vérifié que les résultats sont cohérent

### 1) Chi-2 
#Ca demande un n grand donc c'est normal si le résultat ne fonctionne pas je pense
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

print(compute_asymptotic_pvalue(LR, n_states, dist_type))
#On obtient une p_valeur très très petite donc rejet à priori de H0 good



### 3) Permutation ----
permutation_test <- function(trajectories1, trajectories2,absorbing_state, dist_type, n_states, n1, n2, R) {
  
  T_l <- compute_LR(trajectories1, trajectories2,absorbing_state, dist_type, n_states)
  
  count <- 0  
  
  for (r in 1:R) {
    
    permutation <- sample(1:(n1+n2), n1+n2, replace = FALSE) # Permutation aléatoire
    
    perm_traj1 <- c(trajectories1, trajectories2)[permutation[1:n1]]
    perm_traj2 <- c(trajectories1, trajectories2)[permutation[(n1+1):(n1+n2)]]
    n_states <- n_states
    T_star <- compute_LR(perm_traj1, perm_traj2,absorbing_state, dist_type, n_states)
    
    if (T_star[[2]] <= T_l[[2]]) {
      count <- count + 1
    }
  }
  
  p_perm <- count / R
  
  return(p_perm)
}
print(permutation_test(trajectories1,trajectories2,0,dist_type,n_states,4,54,1000))


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
        times[t] <- rgamma(1, shape=1,rate=1)
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

parametric_bootstrap <- function(trajectories1, trajectories2, n1, n2, n_states, max_transitions, dist_type, R,absorbing_state) {
  
  # On calcule les paramètres du MLE (pour générer les échantillons)
  mle_params <- estimate_SMP_params(c(trajectories1, trajectories2), absorbing_state, dist_type, n_states)
  
  # Calcul de la statistique de test T_l
  T_l <- compute_LR(trajectories1, trajectories2, absorbing_state, dist_type, n_states)
  
  # Compteur pour la p-value
  count <- 0
  for (r in 1:R) {
    
    # Génére n trajectoires indépendantes sous H0 avec les paramètres du MLE
    bootstrap_trajectories1 <- simulate_SMP_bootstrap(n1, n_states, mle_params, dist_type, max_transitions)
    bootstrap_trajectories2 <- simulate_SMP_bootstrap(n2, n_states, mle_params, dist_type, max_transitions)
    
    # Calcule la stat de test T* pour les trajectoires générées
    T_star <- compute_LR(bootstrap_trajectories1, bootstrap_trajectories2, absorbing_state, dist_type, n_states)
    
    # Itérations pour le calcul de la p-value
    if (T_star[[2]] <= T_l[[2]]) {
      count <- count + 1
    }
  }
  
  p_boot <- count / R
  
  return(p_boot)
}

print(parametric_bootstrap(trajectories1,trajectories2,4,54,3,3,dist_type,1000,0))







