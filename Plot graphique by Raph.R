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
# Version alternative avec knitr pour un rendu plus propre (si besoin)
# knitr::kable(legend_table, caption = "Légende des états")
seqdplot(seq_data, with.legend = TRUE, legend.position = "bottom", legend.prop = 0.23)


#Comparer l'entropie avant/après une réforme politique (ex : loi sur les CDD).
# Globalement les sujets ont 23 ans en 1998 et on a eu la réformes 

# Avec l'entropi on verifie la predictibilité des changements d'état

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







