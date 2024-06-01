# bombus preference analysis for using proportion and invasive interactions
# using lme4 

library(tidyverse)
library(lme4)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

setwd("C:/Users/mplat/Documents/thesis")



#------------------------------------
#  Logistic regression 
#------------------------------------

log_reg <- read.csv("log_reg_analyses_file.csv")
pref_data <- log_reg %>% filter(species %in% c('bom_fla', 'bom_mel', 'bom_mix', 'bom_vos', 'bom_imp'))

# Assuming 'species', 'site', and 'date' are categorical variables, you need to convert them to factors
pref_data$species <- as.factor(pref_data$species)
pref_data$site <- as.factor(pref_data$site)
pref_data$date <- as.Date(pref_data$date)
pref_data$proportion_invasive <- as.numeric(pref_data$proportion_invasive)
pref_data$invasive_interaction <- as.numeric(pref_data$invasive_interaction)

pref_data$species <- relevel(pref_data$species, ref = "bom_imp")

# Logistic regression models
#fit0 <- glmer(invasive_interaction ~ proportion_invasive + (1|site), data = pref_data, family = binomial)

fit2 <- glmer(invasive_interaction ~ proportion_invasive + species + proportion_invasive*species + (1|site), data = pref_data, family = binomial)

# Print the summary of the model
#summary(fit0)
summary(fit2)

fit2_anova <- Anova(fit2, type = 3)
fit2_anova

# Goodness of fit?

# plot of the 5 species combined
species_combined_plot <- ggplot(pref_data, aes(x = proportion_invasive, y = invasive_interaction)) + geom_point() + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)
species_combined_plot


#-------------------------------------------------------------------------------
# Plot fit 2 with the interaction term
#-------------------------------------------------------------------------------

n_species = 5
species_names <- c("bom_fla", "bom_mel", "bom_mix", "bom_vos", "bom_imp")
my_palette_reduced <- brewer.pal(n_species, "Paired")

coefficients <- fixef(fit2)

# Generate data for plotting
site_levels <- unique(pref_data$site)
plot_data <- expand.grid(
  proportion_invasive = seq(0, 1, by = 0.01), # adjust the sequence as needed
  species = c('bom_fla', 'bom_mel', 'bom_mix', 'bom_vos', 'bom_imp'),
  site = sample(site_levels, size = length(seq(0, 1, by = 0.01)), replace = TRUE) # Placeholder for site
)

# # Calculate predicted probabilities for each species
# plot_data$predicted_prob <- predict(
#   fit2, 
#   newdata = plot_data, 
#   type = "response"
# )

# Generate data for plotting
plot_data <- expand.grid(
  proportion_invasive = seq(0, 1, by = 0.1), # adjust the step size as needed
  species = c('bom_fla', 'bom_mel', 'bom_mix', 'bom_vos', 'bom_imp'),
  site = rep(site_levels, each = length(seq(0, 1, by = 0.1))) # Repeat each level of 'site'
)

# Calculate predicted probabilities for each species
plot_data$predicted_prob <- predict(
  fit2, 
  newdata = plot_data, 
  type = "response"
)


k <- ggplot(plot_data, aes(x = proportion_invasive, y = predicted_prob, color = species)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.4) + 
  geom_point(data = pref_data, aes(x = proportion_invasive, y = invasive_interaction, color = species), size = 2) +
  scale_color_manual(name = "Species",   # Use scale_color_manual to customize point and line colors
                     labels = c("B. flavifrons", "B. melanopygus", "B. mixtus", "B. vosnesenskii", "B. impatiens"),
                     values = my_palette_reduced) +  # Apply custom palette
  labs(title = "Predicted Probability of Invasive Interaction by Species",
       x = "Proportion Invasive",
       y = "Predicted Probability") +
  theme_minimal()
k


#-------------------------------------------------------------------------------
# table creation

library(rempsyc)

phen_table_data <- data.frame(
  Term = c("(Intercept)", "Invasive flower proportion", "B. flavifrons", "B. melanopygus", "B. mixtus", "B. vosnesenskii"),
  'b' = c(-1.6868, 4.8963, -0.5563, -0.4512, -0.7825, 0.1848),
  'SE' = c(0.3791, 0.2980, 0.2081, 0.3183, 0.1964, 0.2431),
  'z' = c(-4.449, 16.432, -2.674, -1.418, -3.985, 0.760),
  'p' = c(0.001, 0.001, 0.0075, 0.1563, 0.001, 0.4473)
)
print(phen_table_data)
phen_table_data$p <- as.numeric(gsub("[^0-9.]+", "", phen_table_data$p))
nice_table(phen_table_data)
