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
n_species = 5
species_names <- c("bom_fla",  "bom_imp", "bom_mel", "bom_mix", "bom_vos")


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

# get original julian dates back
pred_length <- 100 # divide the calendar year by some number of days
lower_date <- 80
upper_date <- 325
original_julian <- seq(80, 325, length.out = pred_length) # reasonable prediction range
julian_pred <- (original_julian - mean(df$julian)) / sd(df$julian) # unscale the dates

# prep arrays to fill with predictions
count_mean <- array(dim = c(n_species, pred_length))
count_lower95 <- array(dim = c(n_species, pred_length))
count_upper95 <- array(dim = c(n_species, pred_length))


# expected counts
for(i in 1:n_species){
  for(j in 1:pred_length){
    
    if(i != 5){ # i == n_species is the last species (b. impatiens, ie. the reference level)
      # compute for a comparison species
      count_mean[i,j] =
      exp(
        # intercept +
        fit_summary[1,1] +
          # year * 0 for 2022 +
          fit_summary[2,1] * 0 + 
          # species differences in overall abundance
          fit_summary[3+(i-1),1] + 
          # julian scaled with species effect
          (fit_summary[7,1] + fit_summary[9+(i-1),1]) * julian_pred[j] +
          # julian scaled squared with species effect
          (fit_summary[8,1] + fit_summary[13+(i-1),1]) * julian_pred[j]^2 #+
      )
      
      count_lower95[i,j] =
        exp(
          # intercept +
          fit_summary[1,4] +
            # year * 0 for 2022 +
            fit_summary[2,4] * 0 + 
            # species differences in overall abundance
            fit_summary[3+(i-1),4] + 
            # julian scaled with species effect
            (fit_summary[7,4] + fit_summary[9+(i-1),4]) * julian_pred[j] +
            # julian scaled squared with species effect
            (fit_summary[8,4] + fit_summary[13+(i-1),4]) * julian_pred[j]^2 #+
        )
      
      count_upper95[i,j] =
        exp(
          # intercept +
          fit_summary[1,6] +
            # year * 0 for 2022 +
            fit_summary[2,6] * 0 + 
            # species differences in overall abundance
            fit_summary[3+(i-1),6] + 
            # julian scaled with species effect
            (fit_summary[7,6] + fit_summary[9+(i-1),6]) * julian_pred[j] +
            # julian scaled squared with species effect
            (fit_summary[8,6] + fit_summary[13+(i-1),6]) * julian_pred[j]^2 #+
        )
      
    } else { 
      # else compute for the reference level
      count_mean[i,j] =
        exp(
          # intercept +
          fit_summary[1,1] +
            # year * 0 for 2022 +
            fit_summary[2,1] * 0 + 
            # julian scaled
            fit_summary[7,1] * julian_pred[j] +
            # julian scaled squared
            fit_summary[8,1] * julian_pred[j]^2 
        )
      
      count_lower95[i,j] =
        exp(
          # intercept +
          fit_summary[1,4] +
            # year * 0 for 2022 +
            fit_summary[2,4] * 0 + 
            # julian scaled
            fit_summary[7,4] * julian_pred[j] +
            # julian scaled squared
            fit_summary[8,4] * julian_pred[j]^2 
        )
      
      count_upper95[i,j] =
        exp(
          # intercept +
          fit_summary[1,6] +
            # year * 0 for 2022 +
            fit_summary[2,6] * 0 + 
            # julian scaled
            fit_summary[7,6] * julian_pred[j] +
            # julian scaled squared
            fit_summary[8,6] * julian_pred[j]^2 
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

library(RColorBrewer)
my_palette_reduced <- brewer.pal(5, "Blues")

# plot just the expected means (plot p)
p <- ggplot(data = new_df, aes(julian, mean, fill=as.factor(species))) +
  geom_line(size=2) +
  geom_ribbon(aes(
    ymin=lower_95, ymax=upper_95), alpha=0.8) +
  xlim(c(lower_date, upper_date)) +
  ylim(c(0, 20)) +
  theme_bw() +
  ylab("abundance per survey \n(with 95% CI)") +
  xlab("julian date") +
  scale_fill_manual(name = "Species",
                    labels=c("B. flavifrons",
                             "B. melanopygus",
                             "B. mixtus",
                             "B. vosnesenskii", 
                             "B. impatiens"),
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

# would be good to add counts from actual surveys to plot
# but fillled by species colour
# and also to backtransform the dates to real julian days
