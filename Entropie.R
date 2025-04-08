data <- read.table("C:/Users/Valentin/OneDrive/Documents/donnees.txt", 
                   header = TRUE,      # La première ligne contient les noms de colonnes
                   sep = "\t",         # Délimiteur tabulation
                   quote = "",         # Pas de guillemets pour délimiter les chaînes
                   na.strings = "",    # Chaînes vides comme NA
                   fill = TRUE,        # Compléter les lignes trop courtes
                   stringsAsFactors = FALSE)

library(TraMineR)
library(ggplot2)
library(dplyr)
library(tidyr)

seq_data <- data[, paste0("c", 1:94)]

# Définition des labels et états
state_labels <- c(
  "CDI/Indépendant", "CDD", "Alternance", "Mesures publiques", "Intérim", 
  "Chômage", "Inactif", "Service militaire", "Études"
)
state_codes <- 1:9

# Création de l'objet de séquence
seq_object <- seqdef(seq_data, states = state_codes, labels = state_labels)

# Calcul de l'entropie
entropies <- seqient(seq_object)

# Ajout au jeu de données
data$entropy <- entropies

# Tentative de découpage en quartiles
breaks <- unique(quantile(entropies, probs = seq(0, 1, 0.25), na.rm = TRUE))

if (length(breaks) > 2) {
  entropy_groups <- cut(entropies, breaks = breaks, include.lowest = TRUE)
} else {
  entropy_groups <- factor(rep("Entropie constante", length(entropies)))
}

# Conversion en format long pour ggplot
seq_long <- as.data.frame(seq_data)
seq_long$ID <- 1:nrow(seq_long)

seq_long_long <- pivot_longer(seq_long, cols = starts_with("c"), names_to = "Month", values_to = "State")
seq_long_long$Month <- as.numeric(gsub("c", "", seq_long_long$Month))

# Ajout des groupes d'entropie
seq_long_long$entropy_group <- entropy_groups[seq_long_long$ID]

# Affichage avec ggplot2
ggplot(seq_long_long, aes(x = Month, y = reorder(ID, -as.numeric(entropy_group)), fill = factor(State))) +
  geom_tile() +
  facet_wrap(~ entropy_group, scales = "free_y", ncol = 1) +
  labs(
    title = "Trajectoires groupées par quartile d'entropie",
    x = "Mois", y = "Individu", fill = "État"
  ) +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())