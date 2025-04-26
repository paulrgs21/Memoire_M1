rm(list = ls())
library(tidyverse)
library(survival)
library(survminer)
library(TraMineR)

#### Data management 
data <- read.table("C:/Users/rapha/Desktop/Mémoire/donnees.txt", 
                   header = TRUE,      # La première ligne contient les noms de colonnes
                   sep = "\t",         # Délimiteur tabulation
                   quote = "",         # Pas de guillemets pour délimiter les chaînes
                   na.strings = "",    # Chaînes vides comme NA
                   fill = TRUE,        # Compléter les lignes trop courtes
                   stringsAsFactors = FALSE)


sequence_cols <- paste0("c", 1:94)
seq_data <- seqdef(data[, sequence_cols], states = as.character(1:9),
                   labels = c("Permanent employment", "Fixed-term employment", "work-study programs", "Youth employment", "Agency work", "Unemployed", "Inactive", "Military service", "Studies"))

# Distribution par temps
seqdplot(seq_data, 
         with.legend = FALSE, # Désactive la légende intégrée
         cex.legend = 0.8)

# Avec l'entropie on verifie la predictibilité des changements d'état

base_H <- data[data$Q1==1,]
start <- which(colnames(base_H) == "c1")
end <- which(colnames(base_H) == "c94")
base_H <- base_H[, start:end]

base_F <- data[data$Q1==2,]
start <- which(colnames(base_F) == "c1")
end <- which(colnames(base_F) == "c94")
base_F <- base_F[, start:end]

seq_data <- seqdef(base_H)
entropies <- seqient(seq_data)
hist(entropies, main = "Men Entropy Distribution", xlab = "Entropie")
seq_data <- seqdef(base_F)
entropies <- seqient(seq_data)
hist(entropies, main = "Women Entropy Distribution", xlab = "Entropie")


#Graphique à bulle avec les pourcentages de changement d'ETAT
# Charger les données 
trajectories <- trajectoires_smp

# Créer une matrice de transitions
create_transition_matrix <- function(trajectories) {
  # Initialiser une matrice 9x9
  trans_matrix <- matrix(0, nrow=9, ncol=9, 
                         dimnames=list(paste0("From_",1:9), paste0("To_",1:9)))
  
  # Parcourir toutes les trajectoires
  for(traj in trajectories) {
    states <- traj$states
    if(length(states) < 2) next  # Pas de transitions si seul état
    
    # Compter les transitions successives
    for(i in 1:(length(states)-1)) {
      from <- states[i]
      to <- states[i+1]
      trans_matrix[from, to] <- trans_matrix[from, to] + 1
    }
  }
  
  return(trans_matrix)
}

# Calculer la matrice de transition
transition_counts <- create_transition_matrix(trajectories)

# Calculer les pourcentages
transition_percentages <- prop.table(transition_counts, margin=1) * 100

# Afficher les résultats
cat("Matrice des transitions (comptages):\n")
print(transition_counts)

cat("\nMatrice des transitions (pourcentages):\n")
print(round(transition_percentages, 2))

library(DiagrammeR)

# Créer le graphique des transitions
grViz("
digraph transitions {
  graph [layout = neato, overlap = false,fontname = 'Arial Narrow', size = 1, nodesep = 1.2, ranksep = 10]
  
  node [shape = circle, style = filled, fontsize = 6, width = 0.8, height = 0.5, fixedsize = true]

  1 [label = 'Permanent employment']
  2 [label = 'Fixed-term employment']
  3 [label = 'Alternance']
  4 [label = 'Youth employment']
  5 [label = 'Agency work']
  6 [label = 'Unemployed']
  7 [label = 'Inactive']
  8 [label = 'Military']
  9 [label = 'Studies']
  edge [fontname = 'Arial Narrow', fontsize = 5, arrowsize = 0.6,penwidth = 1.2]
  1 -> 2 [label = '23.9%']
  1 -> 6 [label = '38.7%']
  2 -> 1 [label = '45.5%']
  2 -> 6 [label = '33.4%']
  3 -> 1 [label = '43%']
  4 -> 1 [label = '37%']
  5 -> 2 [label = '28,2%']
  5 -> 6 [label = '31,6%']
  7 -> 2 [label = '25,4%']
  7 -> 6 [label = '25,2%']
  8 -> 6 [label = '25,7%']
  4 -> 6 [label = '31,2%']
  6 -> 2 [label = '40.6%']
  9 -> 6 [label = '20.7%']
  9 -> 7 [label = '35.2%']
  
  edge [fontsize = 14]
}
")

















