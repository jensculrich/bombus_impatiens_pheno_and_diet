library(rstan)
library(tidyverse)

## predictor center scaling function
center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

#-----------------------------------------------------------------------
# To visualize species-specific phenology differences:
# Plot bumble bee counts (by species) across date with the upper and lower 95 BCI. 
# We could also make a caterpillar plot with all parameter estimates.

# First we need our original data to get the prediction interval
df <- as.data.frame(read.csv("./data/melissa_2022_2023_phenology_data.csv"))

# see how many detections per species
unique_detections_by_species <- df %>%
  group_by(date, site, species) %>%
  slice(1) %>%
  ungroup() %>%
  group_by(species) %>%
  add_tally() %>%
  slice(1) %>%
  select(species, n)

# IMPORTANT: enter names in alphabetical order
species_names <- c("bom_fla", "bom_imp", "bom_mel", "bom_mix", "bom_vos")
# how many species are we going to consider? Just the main 5 for now
n_species = length(species_names)

df_subset_with_counts <- subset(df, species %in% species_names) %>%
  # let's try looking at just the workers
  filter(caste == "w") %>%
  group_by(date, site, species) %>%
  add_tally() %>%
  slice(1) %>%
  ungroup() %>%
  select(date, year, site, species, n) 

# format data types as factors / scale the dates
df_site_visits <- df %>%
  group_by(date, site) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    site_factor = as.factor(site),
    year_factor = as.factor(year),
    julian = as.numeric(as.POSIXlt(date)$yday), 
    julian_scaled = center_scale(julian),
    julian_scaled_sq = julian_scaled^2)

unique_site_dates <- df_site_visits %>%
  select(-species) # get rid of species because we just have first aplhabetical per site/visit here 

# expand visits by number of species
unique_site_dates <- unique_site_dates[rep(seq_len(nrow(unique_site_dates)), each = n_species), ]

# make a vector of species names possible at each site/visit
species_names_full <- rep(species_names, (nrow(unique_site_dates) / n_species)) 

# add species names to possible visits
all_possible_counts <- as.data.frame(cbind(unique_site_dates, species_names_full)) %>%
  rename("species" = "species_names_full")

# join actual counts to the possible species/site/visit combos
df <- left_join(all_possible_counts, df_subset_with_counts) %>%
  # replace unobserved NA with 0
  mutate(n = replace_na(n, 0),
         #species_factor = relevel(factor(species), ref = 'bom_imp'))
         species_factor = as.factor(species))



#-----------------------------------------------------------------------
# Use the estimates from the fitted poisson model with overdispersion in counts
# read in the saved RDS file containing the model outputs
fit_summary <- rstan::summary(readRDS("./phenology/model_outputs/stan_out_phenology.RDS"))
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

# get original julian dates back
pred_length <- 100 # divide the calendar year by some number of days
lower_date <- 85
upper_date <- 300
original_julian <- seq(lower_date, upper_date, length.out = pred_length) # reasonable prediction range
julian_pred <- (original_julian - mean(df$julian)) / sd(df$julian) # unscale the dates

count_mean <- array(dim = c(n_species, pred_length))
count_lower95 <- array(dim = c(n_species, pred_length))
count_upper95 <- array(dim = c(n_species, pred_length))

# expected counts (for year 2022 at average site)
for(i in 1:n_species){
  for(j in 1:pred_length){
    
      # compute for a comparison species
      count_mean[i,j] =
        exp(
            # intercept +
            fit_summary$summary[1,1] +
            # species-specific intercept
            fit_summary$summary[22+(i-1),1] + 
            # species-specific julian date
            fit_summary$summary[2+(i-1),1] * julian_pred[j] +
            # species-specific julian date^2
            #fit_summary$summary[7+(i-1),1] * julian_pred[j]^2
            (fit_summary$summary[7+(i-1),1] + fit_summary$summary[28,1]) * julian_pred[j]^2 
        )
      
      count_lower95[i,j] =
        exp(
            # intercept +
            fit_summary$summary[1,4] +
            # species-specific intercept
            fit_summary$summary[22+(i-1),4] + 
            # species-specific julian date
            fit_summary$summary[2+(i-1),4] * julian_pred[j] +
            # species-specific julian date^2
            #fit_summary$summary[7+(i-1),4] * julian_pred[j]^2 
            (fit_summary$summary[7+(i-1),4] + fit_summary$summary[28,4]) * julian_pred[j]^2 
        )
      
      count_upper95[i,j] =
        exp(
            # intercept +
            fit_summary$summary[1,8] +
            # species-specific intercept
            fit_summary$summary[22+(i-1),8] + 
            # species-specific julian date
            fit_summary$summary[2+(i-1),8] * julian_pred[j] +
            # species-specific julian date^2
            #fit_summary$summary[7+(i-1),8] * julian_pred[j]^2 
            (fit_summary$summary[7+(i-1),8] + fit_summary$summary[28,8]) * julian_pred[j]^2 
        )
      
  }
}


test <- # intercept +
  fit_summary$summary[1,1] +
  # species-specific intercept
  fit_summary$summary[22+(3-1),1] + 
  # species-specific julian date
  fit_summary$summary[2+(3-1),1] * julian_pred[1] +
  # species-specific julian date^2
  fit_summary$summary[7+(3-1),1] * julian_pred[1]^2 

# unite the means and 95% CI's with species number and julian date
species <- rep(1:n_species, each=pred_length) # species number
# empty vecs to fill with predictions
mean <- vector() 
lower_95 <- vector()
upper_95 <- vector()
# fill predictions by species and add sequentially to the vector
for(i in 1:n_species){
  mean <- c(mean, count_mean[i,]) 
  lower_95 <- c(lower_95, count_lower95[i,])
  upper_95 <- c(upper_95, count_upper95[i,])
}
# unite as a dataframe and rename the repped julian dates
new_df <- as.data.frame(cbind(species, julian_pred, rep(original_julian, n_species),
                              mean, lower_95, upper_95)) %>%
  rename_with(.cols = 3, ~"julian")

# plotting palette
library(RColorBrewer)
my_palette_reduced <- brewer.pal(n_species, "Set2")

# plot just the expected means (plot p)
p <- ggplot(data = new_df, aes(julian, mean, fill=as.factor(species))) +
  geom_ribbon(aes(
    ymin=lower_95, ymax=upper_95), alpha=0.8) +
  #geom_line(size=2) +
  xlim(c(lower_date, upper_date)) +
  ylim(c(0, 225)) +
  theme_bw() +
  ylab("abundance per survey \n(with 95% CI)") +
  xlab("julian date") +
  scale_fill_manual(name = "Species",
                    labels=c("B. flavifrons",
                             "B. impatiens",
                             "B. melanopygus",
                             "B. mixtus",
                             "B. vosnesenskii"),
                    values = my_palette_reduced) +
  theme(legend.position = c(0.15, 0.8),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p










# old
df1 <- as.data.frame(cbind((rep(1, pred_length)), pred_seq,
                           count_mean[1,],
                           #lower_50[1,], upper_50[1,],
                           count_lower95[1,], count_upper95[1,]))  %>%
  rename("species" = "V1",
         "mean" = "V3",
         #"lower_50" = "V4",
         #"upper_50" = "V5",
         "lower_95" = "V4",
         "upper_95" = "V5") 

df2 <- as.data.frame(cbind((rep(2, pred_length)), pred_seq,
                           count_mean[2,],
                           #lower_50[2,], upper_50[2,],
                           count_lower95[2,], count_upper95[2,]))  %>%
  rename("species" = "V1",
         "mean" = "V3",
         #"lower_50" = "V4",
         #"upper_50" = "V5",
         "lower_95" = "V4",
         "upper_95" = "V5") 

df3 <- as.data.frame(cbind((rep(3, pred_length)), pred_seq,
                           count_mean[3,],
                           #lower_50[3,], upper_50[3,],
                           count_lower95[3,], count_upper95[3,]))  %>%
  rename("species" = "V1",
         "mean" = "V3",
         #"lower_50" = "V4",
         #"upper_50" = "V5",
         "lower_95" = "V4",
         "upper_95" = "V5") 

df4 <- as.data.frame(cbind((rep(4, pred_length)), pred_seq,
                           count_mean[4,],
                           #lower_50[4,], upper_50[4,],
                           count_lower95[4,], count_upper95[4,]))  %>%
  rename("species" = "V1",
         "mean" = "V3",
         #"lower_50" = "V4",
         #"upper_50" = "V5",
         "lower_95" = "V4",
         "upper_95" = "V5")

df5 <- as.data.frame(cbind((rep(5, pred_length)), pred_seq,
                           count_mean[5,],
                           #lower_50[4,], upper_50[4,],
                           count_lower95[5,], count_upper95[5,]))  %>%
  rename("species" = "V1",
         "mean" = "V3",
         #"lower_50" = "V4",
         #"upper_50" = "V5",
         "lower_95" = "V4",
         "upper_95" = "V5")

new_df <- rbind(df1, df2, df3, df4, df5) 

df <- df %>%
  mutate(pred_seq = julian_scaled,
         mean = as.numeric(n)) 

library(RColorBrewer)
my_palette_reduced <- brewer.pal(n_species, "Blues")

p <- ggplot(data = new_df, aes(pred_seq, mean)) +
  geom_line(data=df1, size=2) +
  geom_ribbon(data=df1, aes(ymin=lower_95, ymax=upper_95, fill = as.factor(species)), alpha=0.8) +
  geom_line(data=df2, size=2) +
  geom_ribbon(data=df2, aes(ymin=lower_95, ymax=upper_95, fill = as.factor(species)), alpha=0.8) + 
  geom_line(data=df3, size=2) +
  geom_ribbon(data=df3, aes(ymin=lower_95, ymax=upper_95, fill = as.factor(species)), alpha=0.8) +  
  geom_line(data=df4, size=2) +
  geom_ribbon(data=df4, aes(ymin=lower_95, ymax=upper_95, fill = as.factor(species)), alpha=0.8) +  
  xlim(c(-3, 3.5)) +
  ylim(c(0, 600)) +
  theme_bw() +
  ylab("expected mean abundance per survey \n(with 10 and 90% quantiles)") +
  xlab("julian date scaled") +
  scale_fill_manual(name = "Species",
                    labels=c("B. flavifrons",
                             "B. impatiens",
                             "B. melanopygus",
                             "B. mixtus",
                             "B. vosnesenskii"
                             ),
                    values=my_palette_reduced) +
  theme(legend.position = c(0.15, 0.8),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size = 20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p
