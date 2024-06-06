library(rstan)
library(tidyverse)

## predictor center scaling function
center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

## predictor link function
inv_logit <- function(x) exp(x)/(1+exp(x))

#------------------------------------------------------------------------------
# load and clean the data

log_reg <- read.csv("./data/melissa_2022_2023_logistic_analyses_file.csv")
pref_data <- log_reg %>% filter(species %in% c('bom_fla', 'bom_mel', 'bom_mix', 'bom_vos', 'bom_imp')) %>%
  
  # currently there is NA data? Why?
  filter(!is.na(proportion_invasive))

# Assuming 'species', 'site', and 'date' are categorical variables, you need to convert them to factors
pref_data$species <- as.factor(pref_data$species)
pref_data$site <- as.factor(pref_data$site)
pref_data$date <- as.Date(pref_data$date)
pref_data$proportion_invasive <- as.numeric(pref_data$proportion_invasive)
pref_data$invasive_interaction <- as.numeric(pref_data$invasive_interaction)

pref_data$species <- relevel(pref_data$species, ref = "bom_imp")
n_species <- length(unique(pref_data$species))

#------------------------------------------------------------------------------
# Scale the proportion of invasive plants
# To get scaled dates we need to have just one event per site visit

pref_data <- pref_data %>%
  group_by(site, month_day) %>%
  mutate(survey = cur_group_id()) %>%
  ungroup()

pref_data <- pref_data %>%
  group_by(site, month_day) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(scaled_prop_nvsv = center_scale(proportion_invasive)) %>%
  select(survey, scaled_prop_nvsv) %>%
  left_join(pref_data, ., by="survey")

#-----------------------------------------------------------------------
# Use the estimates from the fitted poisson model with overdispersion in counts
# read in the saved RDS file containing the model outputs
fit_summary <- rstan::summary(readRDS("./diet_choice/model_outputs/stan_out_diet_choice.RDS"))
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

# get original proportion of invasive plants back from the data
pred_length <- 100 # divide the spectrum by some interval
lower_prop_nvsv <- 0
upper_prop_nvsv <- 1

temp1 <- pref_data %>%
  group_by(site, month_day) %>%
  slice(1) %>%
  ungroup() %>%
  select(proportion_invasive)

mean <- mean(temp1$proportion_invasive)
sd <- sd(temp1$proportion_invasive)

original_nvsv <- seq(lower_prop_nvsv, upper_prop_nvsv, length.out = pred_length) # reasonable prediction range
nvsv_pred <- (original_nvsv - mean) / sd # unscale the dates

choice_mean <- array(dim = c(n_species, pred_length))
choice_lower95 <- array(dim = c(n_species, pred_length))
choice_upper95 <- array(dim = c(n_species, pred_length))
choice_lower50 <- array(dim = c(n_species, pred_length))
choice_upper50 <- array(dim = c(n_species, pred_length))

# expected choices
for(i in 1:n_species){
  for(j in 1:pred_length){
    
    if(i != 1){ # i != b. impatiens, ie. the reference level
      # compute for a comparison species
      choice_mean[i,j] =
        inv_logit(
            # intercept +
            fit_summary$summary[1,1] +
            # species differences in choice at mean level of invasive plants
            fit_summary$summary[1+(i-1),1] + 
            # species differences in choice as invasive plants increase
            (fit_summary$summary[6,1] + fit_summary$summary[6+(i-1),1]) * nvsv_pred[j]
        )
      
      choice_lower95[i,j] =
        inv_logit(
          # intercept +
          fit_summary$summary[1,4] +
            # species differences in choice at mean level of invasive plants
            fit_summary$summary[1+(i-1),4] + 
            # species differences in choice as invasive plants increase
            (fit_summary$summary[6,4] + fit_summary$summary[6+(i-1),4]) * nvsv_pred[j]
        )
      
      choice_upper95[i,j] =
        inv_logit(
          # intercept +
          fit_summary$summary[1,8] +
            # species differences in choice at mean level of invasive plants
            fit_summary$summary[1+(i-1),8] + 
            # species differences in choice as invasive plants increase
            (fit_summary$summary[6,8] + fit_summary$summary[6+(i-1),8]) * nvsv_pred[j]
        )
      
      choice_lower50[i,j] =
        inv_logit(
          # intercept +
          fit_summary$summary[1,5] +
            # species differences in choice at mean level of invasive plants
            fit_summary$summary[1+(i-1),5] + 
            # species differences in choice as invasive plants increase
            (fit_summary$summary[6,5] + fit_summary$summary[6+(i-1),5]) * nvsv_pred[j]
        )
      
      choice_upper50[i,j] =
        inv_logit(
          # intercept +
          fit_summary$summary[1,7] +
            # species differences in choice at mean level of invasive plants
            fit_summary$summary[1+(i-1),7] + 
            # species differences in choice as invasive plants increase
            (fit_summary$summary[6,7] + fit_summary$summary[6+(i-1),7]) * nvsv_pred[j]
        )
      
    } else { 
      
      # else compute for the reference level
      choice_mean[i,j] =
        inv_logit(
            # intercept +
            fit_summary$summary[1,1] +
            # species differences in choice as invasive plants increase
            (fit_summary$summary[6,1]) * nvsv_pred[j]
        )
      
      choice_lower95[i,j] =
        inv_logit(
          # intercept +
          fit_summary$summary[1,4] +
            # species differences in choice as invasive plants increase
            (fit_summary$summary[6,4]) * nvsv_pred[j]
        )
      
      choice_upper95[i,j] =
        inv_logit(
          # intercept +
          fit_summary$summary[1,8] +
            # species differences in choice as invasive plants increase
            (fit_summary$summary[6,8]) * nvsv_pred[j]
        )
      
      choice_lower50[i,j] =
        inv_logit(
          # intercept +
          fit_summary$summary[1,5] +
            # species differences in choice as invasive plants increase
            (fit_summary$summary[6,5]) * nvsv_pred[j]
        )
      
      choice_upper50[i,j] =
        inv_logit(
          # intercept +
          fit_summary$summary[1,7] +
            # species differences in choice as invasive plants increase
            (fit_summary$summary[6,7]) * nvsv_pred[j]
        )
    }
    
  }
}

# unite the means and 95% CI's with species number and julian date
species <- rep(1:n_species, each=pred_length) # species number
# empty vecs to fill with predictions
mean <- vector() 
lower_95 <- vector()
upper_95 <- vector()
lower_50 <- vector()
upper_50 <- vector()

# fill predictions by species and add sequentially to the vector
for(i in 1:n_species){
  mean <- c(mean, choice_mean[i,]) 
  lower_95 <- c(lower_95, choice_lower95[i,])
  upper_95 <- c(upper_95, choice_upper95[i,])
  lower_50 <- c(lower_50, choice_lower50[i,])
  upper_50 <- c(upper_50, choice_upper50[i,])
}
# unite as a dataframe and rename the repped julian dates
new_df <- as.data.frame(cbind(species, nvsv_pred, rep(original_nvsv, n_species),
                              mean, lower_95, upper_95, lower_50, upper_50)) %>%
  rename_with(.cols = 3, ~"prop_nvsv") %>%
  mutate(species_new = as.factor(species)) %>%
  mutate(species_new_ordered = factor(species_new, 
                                      levels = c("1","4", "2", "5", "3")))

temp <- pref_data %>%
  mutate(species_new = as.factor(as.numeric(species))) %>%
  mutate(species_new_ordered = factor(species_new, 
                                      levels = c("1","4", "2", "5", "3")))

library(RColorBrewer)

#my_palette <- brewer.pal(5, "Blues")
#my_palette[1] <- "#7C0000"

my_palette <- hcl.colors(5, "Inferno")

# plot just the expected means (plot p)
p <- ggplot(data = new_df, aes(nvsv_pred, mean, fill=as.factor(species_new_ordered))) +
  
  geom_ribbon(aes(
  ymin=lower_50, ymax=upper_50), alpha=0.5) +
  #geom_ribbon(aes(
  #  ymin=lower_95, ymax=upper_95), alpha=0.7) +
 #geom_line(size=2) +
  #geom_ribbon(aes(
  #  ymin=mean-0.01, ymax=mean), alpha=1) +
  #xlim(c(lower_prop_nvsv, upper_prop_nvsv)) +
  #xlim(-1, 1) +
  geom_point(data = temp, aes(
    x = scaled_prop_nvsv, y = invasive_interaction, 
    color = species_new_ordered), size = 2) +
  scale_x_continuous(limits = c(-0.81, 2.05), 
                     breaks = c(-0.5, 0, 0.5, 1, 1.5, 2),
                     labels = signif(x = c(-0.5, 0, 0.5, 1, 1.5, 2), digits=2)) +
  theme_bw() +
  ylab("Prob(interacting with invasive plant)\n(50% BCI)") +
  xlab("Proportion of invasive species in the plant community\n(z-score scaled)") +
  scale_fill_manual(name = "Species",
                    labels=c("B. impatiens",
                             "B. mixtus",
                             "B. flavifrons",
                             "B. vosnesenskii", 
                             "B. melanopygus"),
                    values = my_palette) +
  scale_colour_manual(name = "Species",
                    labels=c("B. impatiens",
                             "B. mixtus",
                             "B. flavifrons",
                             "B. vosnesenskii", 
                             "B. melanopygus"),
                    values = my_palette) +
  scale_y_continuous(limits = c(0, 1), 
                     breaks = c(0, 0.5, 1),
                     labels = scales::percent) +
  theme(legend.position = c(0.15, 0.7),
        legend.title=element_text(size=18),
        legend.text=element_text(size=16),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

#-------------------------------------------------------------------------------
# make caterpillar plots
## --------------------------------------------------
## Plot ecological paramter means and variation

# parameter means
# let's plot species effects of julian date (for 5 species)
X_eco <- c(1, 2, 3, 4, 5, 6, 7, 8) # params of interest

# mean of eco params
Y_eco <- c(
  fit_summary$summary[2,1],
  fit_summary$summary[3,1],
  fit_summary$summary[4,1], 
  fit_summary$summary[5,1],
  fit_summary$summary[7,1],
  fit_summary$summary[8,1],
  fit_summary$summary[9,1], 
  fit_summary$summary[10,1]
)

# confidence intervals
lower_95_eco <- c( 
  fit_summary$summary[2,4],
  fit_summary$summary[3,4],
  fit_summary$summary[4,4], 
  fit_summary$summary[5,4],
  fit_summary$summary[7,4],
  fit_summary$summary[8,4],
  fit_summary$summary[9,4], 
  fit_summary$summary[10,4]
)

upper_95_eco <- c(
  fit_summary$summary[2,8],
  fit_summary$summary[3,8],
  fit_summary$summary[4,8], 
  fit_summary$summary[5,8],
  fit_summary$summary[7,8],
  fit_summary$summary[8,8],
  fit_summary$summary[9,8], 
  fit_summary$summary[10,8]
)

# confidence intervals
lower_50_eco <- c(
  fit_summary$summary[2,5],
  fit_summary$summary[3,5],
  fit_summary$summary[4,5], 
  fit_summary$summary[5,5],
  fit_summary$summary[7,5],
  fit_summary$summary[8,5],
  fit_summary$summary[9,5], 
  fit_summary$summary[10,5]
)

upper_50_eco <- c(
  fit_summary$summary[2,7],
  fit_summary$summary[3,7],
  fit_summary$summary[4,7], 
  fit_summary$summary[5,7],
  fit_summary$summary[7,7],
  fit_summary$summary[8,7],
  fit_summary$summary[9,7], 
  fit_summary$summary[10,7]
)


df_estimates_eco <- as.data.frame(cbind(X_eco, Y_eco, 
                                        lower_95_eco, upper_95_eco,
                                        lower_50_eco, upper_50_eco))

df_estimates_eco$X_eco <- as.factor(df_estimates_eco$X_eco)

df_estimates_eco1 <- df_estimates_eco[1:4,]
df_estimates_eco2 <- df_estimates_eco[5:8,]

## --------------------------------------------------
## Draw ecological parameter plot

(s <- ggplot(df_estimates_eco1) +
   theme_bw() +
   # scale_color_viridis(discrete=TRUE) +
   scale_x_discrete(name="", breaks = c(1, 2, 3, 4),
                    labels=c(
                      "B. flavifrons",
                      "B. melanopygus",
                      "B. mixtus", 
                      "B. vosnesenskii")
   ) +
   scale_y_continuous("Effect of species identity on intercept\n(relative to B. impatiens) (logit-scaled)",
                      limits = c(-2, 1), breaks = c(-2, -1, 0, 1)) +
   guides(color = guide_legend(title = "")) +
   geom_hline(yintercept = 0, lty = "dashed") +
   theme(legend.text=element_text(size=10),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18, angle=0, vjust=0.5),
         axis.title.x = element_text(size = 18),
         axis.title.y = element_text(size = 18),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   coord_flip()
)

s <- s +
  geom_errorbar(aes(x=X_eco, ymin=lower_95_eco, ymax=upper_95_eco),
                color="black",width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X_eco, ymin=lower_50_eco, ymax=upper_50_eco),
                color="black",width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X_eco, y=Y_eco),
             size = 5, alpha = 0.8) 


s

(t <- ggplot(df_estimates_eco2) +
    theme_bw() +
    # scale_color_viridis(discrete=TRUE) +
    scale_x_discrete(name="", breaks = c(5, 6, 7, 8),
                     labels=c(
                       "B. flavifrons",
                       "B. melanopygus",
                       "B. mixtus", 
                       "B. vosnesenskii")
    ) +
    scale_y_continuous("Effect of species identity on response to increasing invasive plants\n(relative to B. impatiens) (logit-scaled)",
                       limits = c(-2, 1), breaks = c(-2, -1, 0, 1)) +
    guides(color = guide_legend(title = "")) +
    geom_hline(yintercept = 0, lty = "dashed") +
    theme(legend.text=element_text(size=10),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18, angle=0, vjust=0.5),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    coord_flip()
)

t <- t +
  geom_errorbar(aes(x=X_eco, ymin=lower_95_eco, ymax=upper_95_eco),
                color="black",width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X_eco, ymin=lower_50_eco, ymax=upper_50_eco),
                color="black",width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X_eco, y=Y_eco),
             size = 5, alpha = 0.8) 


t

#------------------------------------------------------------------------------
# Create tables with model data

library(rempsyc)

diet_table_data <- data.frame(
  Term = c("Intercept", "B. flavifrons", "B. mixtus", "B. melanopygus", "B. vosnesenskii"),
  'mean' = c(fit_summary$summary[1,1], fit_summary$summary[2,1], 
             fit_summary$summary[3,1], fit_summary$summary[4,1], fit_summary$summary[5,1]),
  'lower 95 BCI' = c(fit_summary$summary[1,4], fit_summary$summary[2,4], 
                     fit_summary$summary[3,4], fit_summary$summary[4,4], fit_summary$summary[5,4]),
  'upper 95 BCI' = c(fit_summary$summary[1,8], fit_summary$summary[2,8], 
                     fit_summary$summary[3,8], fit_summary$summary[4,8], fit_summary$summary[5,8])
)
print(diet_table_data)
nice_table(diet_table_data)

diet_table_data2 <- data.frame(
  Term = c("Prop.\ninvasive plants", "B. flavifrons", "B. mixtus", "B. melanopygus", "B. vosnesenskii"),
  'mean' = c(fit_summary$summary[6,1], fit_summary$summary[7,1], 
             fit_summary$summary[8,1], fit_summary$summary[9,1], fit_summary$summary[10,1]),
  'lower 95 BCI' = c(fit_summary$summary[6,4], fit_summary$summary[7,4], 
                     fit_summary$summary[8,4], fit_summary$summary[9,4], fit_summary$summary[10,4]),
  'upper 95 BCI' = c(fit_summary$summary[6,8], fit_summary$summary[7,8], 
                     fit_summary$summary[8,8], fit_summary$summary[9,8], fit_summary$summary[10,8])
)
print(diet_table_data2)
nice_table(diet_table_data2)

diet_table_data3 <- data.frame(
  Term = c("(Intercept)", "B. flavifrons", "B. mixtus", "B. melanopygus", "B. vosnesenskii",
           "(Prop.\ninvasive plants)", "B. flavifrons", "B. mixtus", "B. melanopygus", "B. vosnesenskii"),
  'mean' = c(fit_summary$summary[1,1], fit_summary$summary[2,1], 
             fit_summary$summary[3,1], fit_summary$summary[4,1], fit_summary$summary[5,1],
             fit_summary$summary[6,1], fit_summary$summary[7,1], 
             fit_summary$summary[8,1], fit_summary$summary[9,1], fit_summary$summary[10,1]),
  'lower 95 BCI' = c(fit_summary$summary[1,4], fit_summary$summary[2,4], 
                     fit_summary$summary[3,4], fit_summary$summary[4,4], fit_summary$summary[5,4],
                     fit_summary$summary[6,4], fit_summary$summary[7,4], 
                     fit_summary$summary[8,4], fit_summary$summary[9,4], fit_summary$summary[10,4]),
  'upper 95 BCI' = c(fit_summary$summary[1,8], fit_summary$summary[2,8], 
                     fit_summary$summary[3,8], fit_summary$summary[4,8], fit_summary$summary[5,8],
                     fit_summary$summary[6,8], fit_summary$summary[7,8], 
                     fit_summary$summary[8,8], fit_summary$summary[9,8], fit_summary$summary[10,8])
)
print(diet_table_data3)
nice_table(diet_table_data3)
