rm(list = ls())
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
                   labels = c("CDI", "CDD", "Alt", "Public", "Intérim", "Chômage", "Inactif", "Service", "Études"))

# Distribution par temps
seqdplot(seq_data, with.legend = TRUE, legend.position = "bottom", legend.prop = 0.22)


#Comparer l'entropie avant/après une réforme politique (ex : loi sur les CDD).
# Globalement les sujets ont 23 ans en 1998 et on a eu la réformes 

# Avec l'entropi on verifie la predictibilité des changements d'état
# Exemple on voit bien que en fin de carrière l'etat est beaucoup plus prévisible exemple
start <- which(colnames(data) == "c8")
end <- which(colnames(data) == "c18")
Debut_carriere <- data[, start:end]

start <- which(colnames(data) == "c12")
end <- which(colnames(data) == "c29")
Fin_carriere <- data[, start:end]


seq_data <- seqdef(Debut_carriere)
entropies <- seqient(seq_data)
hist(entropies, main = "Distribution des entropies", xlab = "Entropie")
seq_data <- seqdef(Fin_carriere)
entropies <- seqient(seq_data)
hist(entropies, main = "Distribution des entropies", xlab = "Entropie")
