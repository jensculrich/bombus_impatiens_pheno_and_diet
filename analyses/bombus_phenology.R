# bombus phenology analysis and plot for melissa's data
# using rstanarm to estimate parameter distributions with Bayesian implementation

library(tidyverse)
library(rstanarm)
library(bayesplot)

## predictor center scaling function
center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

#------------------------------------------------------------------------------
# initial data cleaning

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

# how many species are we going to consider? Just the main 4 for now
n_species = 8
species_names <- c("bom_fla", "bom_mel", "bom_imp", "bom_mix", 
                   "bom_sit", "bom_cal", "bom_ruf", "bom_vos")


df_subset_with_counts <- subset(df, species %in% species_names) %>%
  # let's try looking at just the workers
  #filter(caste == "w") %>%
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

#------------------------------------------------------------------------------
# Add zeroes here for zero counts on surveys (because you looked for and found zero when there is no row)

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
         species_factor = relevel(factor(species), ref = 'bom_imp'))

#------------------------------------------------------------------------------
# fit model

# set plotting colour scheme
color_scheme_set("pink")

# poisson model
stan_out <- rstanarm::stan_glmer(data = df, 
                                 n ~ (1|site_factor) + year_factor +
                                   species_factor +
                                   julian_scaled*species_factor +
                                   julian_scaled_sq*species_factor,
                                 family = "poisson")

pp_check(stan_out, plotfun = "stat_grouped", group = "species_factor")
pp_check(stan_out, plotfun = "stat_grouped",  stat = "max", group = "species_factor")
pp_check(stan_out, plotfun = "stat_grouped",  stat = "min", group = "species_factor")


mcmc_areas(stan_out, regex_pars = c("species_factor", 
                                    "julian_scaled",
                                    "julian_scaled_sq",
                                    "julian_scaled*species_factor",
                                    "julian_scaled_sq*species_factor"), 
           prob = 0.89)


# negative binomial model
stan_out_nb <- rstanarm::stan_glmer.nb(data = df, 
                                       n ~ (1|site_factor) + year_factor +
                                         species_factor +
                                         julian_scaled*species_factor +
                                         julian_scaled_sq*species_factor)

pp_check(stan_out_nb, plotfun = "stat_grouped", group = "species_factor")
pp_check(stan_out_nb, plotfun = "stat_grouped",  stat = "max", group = "species_factor")
pp_check(stan_out_nb, plotfun = "stat_grouped",  stat = "min", group = "species_factor")


mcmc_areas(stan_out_nb, regex_pars = c("species_factor", 
                                       "julian_scaled",
                                       "julian_scaled_sq",
                                       "julian_scaled*species_factor",
                                       "julian_scaled_sq*species_factor"), 
           prob = 0.89)


#-----------------------------------------------------------------------
# Plot means with the upper and lowers 90 [,4] [,6] 
# uppers and lowers for all effects not just the species effects

# Just use the poisson fit for now
fit_summary <- rstan::summary(stan_out)
View(cbind(1:nrow(fit_summary), fit_summary)) # View to see which row corresponds to the parameter of interest

pred_seq <- seq(-3, 3, 0.1)
pred_length <- length(pred_seq)

count_mean <- array(dim = c(n_species, pred_length))
count_lower95 <- array(dim = c(n_species, pred_length))
count_upper95 <- array(dim = c(n_species, pred_length))


# expected counts
for(i in 1:n_species){
  for(j in 1:pred_length){
    
    if(i != n_species){ # i == n_species is the last species (b. impatiens, ie. the reference level)
      # compute for a comparison species
      count_mean[i,j] =
        exp(
          # intercept +
          fit_summary[1,1] +
            # year * 0 for 2022 +
            fit_summary[2,1] * 0 + 
            # start at first row of species effects
            # then each next species will be + 4 + i
            fit_summary[3+(i-1),1] + 
            # julian scaled
            fit_summary[6,1] * pred_seq[j] +
            # julian scaled squared
            fit_summary[7,1] * pred_seq[j]^2 +
            # julian scaled * species
            fit_summary[8+(i-1),1] * pred_seq[j] +
            # julian scaled squared * species
            fit_summary[11+(i-1),1] * pred_seq[j]^2
        )
      
      count_lower95[i,j] =
        exp(
          # intercept +
          fit_summary[1,4] +
            # year * 0 for 2022 +
            fit_summary[2,4] * 0 + 
            # start at first row of species effects
            # then each next species will be + 4 + i
            fit_summary[3+(i-1),4] + 
            # julian scaled
            fit_summary[6,4] * pred_seq[j] +
            # julian scaled squared
            fit_summary[7,4] * pred_seq[j]^2 +
            # julian scaled * species
            fit_summary[8+(i-1),4] * pred_seq[j] +
            # julian scaled squared * species
            fit_summary[11+(i-1),4] * pred_seq[j]^2
        )
      
      count_upper95[i,j] =
        exp(
          # intercept +
          fit_summary[1,6] +
            # year * 0 for 2022 +
            fit_summary[2,6] * 0 + 
            # start at first row of species effects
            # then each next species will be + 4 + i
            fit_summary[3+(i-1),6] + 
            # julian scaled
            fit_summary[6,6] * pred_seq[j] +
            # julian scaled squared
            fit_summary[7,6] * pred_seq[j]^2 +
            # julian scaled * species
            fit_summary[8+(i-1),6] * pred_seq[j] +
            # julian scaled squared * species
            fit_summary[11+(i-1),6] * pred_seq[j]^2
        )
      
    } else { 
      # else compute for the reference level
      count_mean[i,j] =
        exp(
          # intercept +
          fit_summary[1,1] +
            # year * 0 for 2022 +
            fit_summary[2,1] * 0 + 
            # start at first row of species effects
            # then each next species will be + 4 + i
            #fit_summary[3+(i-1),1] + 
            # julian scaled
            fit_summary[6,1] * pred_seq[j] +
            # julian scaled squared
            fit_summary[7,1] * pred_seq[j]^2 #+
          # julian scaled * species
          #fit_summary[8+(i-1),1] * pred_seq[j] +
          # julian scaled squared * species
          #fit_summary[11+(i-1),1] * pred_seq[j]^2
        )
      
      count_lower95[i,j] =
        exp(
          # intercept +
          fit_summary[1,4] +
            # year * 0 for 2022 +
            fit_summary[2,4] * 0 + 
            # start at first row of species effects
            # then each next species will be + 4 + i
            #fit_summary[3+(i-1),1] + 
            # julian scaled
            fit_summary[6,4] * pred_seq[j] +
            # julian scaled squared
            fit_summary[7,4] * pred_seq[j]^2 #+
          # julian scaled * species
          #fit_summary[8+(i-1),1] * pred_seq[j] +
          # julian scaled squared * species
          #fit_summary[11+(i-1),1] * pred_seq[j]^2
        )
      
      count_upper95[i,j] =
        exp(
          # intercept +
          fit_summary[1,6] +
            # year * 0 for 2022 +
            fit_summary[2,6] * 0 + 
            # start at first row of species effects
            # then each next species will be + 4 + i
            #fit_summary[3+(i-1),1] + 
            # julian scaled
            fit_summary[6,6] * pred_seq[j] +
            # julian scaled squared
            fit_summary[7,6] * pred_seq[j]^2 #+
          # julian scaled * species
          #fit_summary[8+(i-1),1] * pred_seq[j] +
          # julian scaled squared * species
          #fit_summary[11+(i-1),1] * pred_seq[j]^2
        )
    }
    
  }
}

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

new_df <- rbind(df1, df2, df3, df4) 

df <- df %>%
  mutate(pred_seq = julian_scaled,
         mean = as.numeric(n)) 

library(RColorBrewer)
my_palette_reduced <- brewer.pal(4, "Blues")

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
  ylim(c(0, 40)) +
  theme_bw() +
  ylab("expected mean abundance per survey \n(with 10 and 90% quantiles)") +
  xlab("julian date scaled") +
  scale_fill_manual(name = "Species",
                    labels=c("B. flavifrons",
                             "B. mixtus",
                             "B. vosnesenskii",
                             "B. impatiens"),
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

# would be good to add counts from actual surveys to plot
# but fillled by species colour
# and also to backtransform the dates to real julian days
