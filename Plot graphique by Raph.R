rm(list = ls())

library(MASS)
library(ggplot2)
library(vcd)
library(dplyr)
#### Data management 
data <- read.table("C:/Users/rapha/Desktop/Mémoire/donnees.txt", 
                   header = TRUE,      # La première ligne contient les noms de colonnes
                   sep = "\t",         # Délimiteur tabulation
                   quote = "",         # Pas de guillemets pour délimiter les chaînes
                   na.strings = "",    # Chaînes vides comme NA
                   fill = TRUE,        # Compléter les lignes trop courtes
                   stringsAsFactors = FALSE)


#On va faire quelques graph pour essayer de deviner des corrélations les parametres.
# Sélection des variables numériques
num_vars <- data[, sapply(data, is.numeric)]

# Calcul de la matrice de corrélation
cor_matrix <- cor(num_vars, use = "pairwise.complete.obs") #Trop de variables essayer de séparer ?



# je vais essayer de deviner des corrélations possible, d'abord entre les variables qualitative.
# corrélation entre sexe et petits boulots pendant les études:
ggplot(data, aes(x = as.factor(Q1), fill = as.factor(Q52))) +
  geom_bar(position = "fill") +  # Position "fill" pour voir la proportion
  labs(title = "Lien entre le sexe et les petits boulots", 
       x = "Sexe (1 = Homme, 2 = Femme)", 
       y = "Proportion", 
       fill = "Petits boulots") +
  scale_fill_manual(values = c("red", "yellow", "blue")) +
  theme_minimal()





# Maintenant on va plot quelques relations entre les variables qualitative et quantitative 
ggplot(data, aes(x = as.factor(Q1), fill = as.factor(c1))) +
  geom_bar(position = "fill") +  # Position "fill" pour voir la proportion
  labs(title = "Lien entre le sexe et les petits boulots", 
       x = "Sexe (1 = Homme, 2 = Femme)", 
       y = "Proportion", 
       fill = "Petits boulots") +
  scale_fill_manual(values = c("red", "yellow", "blue")) +
  theme_minimal()
ggplot(data(aes(x = as.factor)))


library(ggplot2)
library(tidyr)
library(dplyr)

  
# Transformer les données en format long
donnees_long <- data %>% 
  pivot_longer(
      cols = c1:c30,
      names_to = "variable",
      values_to = "valeur"
  ) %>% 
  mutate(variable = factor(variable, levels = paste0("c", 1:30)))
  
# Créer le graphique
ggplot(donnees_long, aes(x = variable, y = valeur, group = id)) +
    geom_line(alpha = 0.2, color = "steelblue") + 
    scale_y_continuous(breaks = 1:9, limits = c(1, 9)) +
    labs(x = "Variables c1 à c30", 
         y = "Valeurs",
         title = "Profil des individus à travers les variables") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))




