
#Fonction qui calcule la log-vraisemblance d'un ensemble de trajectoires d'un SMP.
#Params : contient alpha : les proba initiales des états, P la matrice de transition entre les états, omega la liste des paramètres des lois de durée associées aux états
#trajectories : liste avec les trajectoires observées
#dist_type : type de distribution pour les temps de séjour, ici prend que gamma ou Weibull.
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
