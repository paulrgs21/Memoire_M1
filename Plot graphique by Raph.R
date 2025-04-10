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


entropies <- seqient(seq_data)
hist(entropies, main = "Distribution des entropies", xlab = "Entropie")
