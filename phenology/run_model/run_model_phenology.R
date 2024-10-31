library(rstan)
library(tidyverse)

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

# how many species are we going to consider? Just the main 5 for now
species_names <- c("bom_fla",  "bom_imp", "bom_mel", "bom_mix", "bom_vos")
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
# get unique species and sites

n_sites <- length(unique(sites <- df %>% 
  group_by(site) %>%
  slice(1) %>%
  pull(site)))

n_species <- length(unique(species <- df %>% 
                           group_by(species) %>%
                           slice(1) %>%
                           pull(species)))

#------------------------------------------------------------------------------
# Prepare data for stan model

# data to feed to the model
# dummy_variables <- model.matrix(~ species_factor, data = df)
X <- model.matrix(n ~ species_factor, data = df)

N <- nrow(df) # number of pairs
y <- df$n # outcomes (counts)
n_species <- n_species # number of species
n_survey_events <- N/n_species # for PPC
julian_scaled <- df$julian_scaled
julian_scaled_sq <- df$julian_scaled_sq
n_sites <- n_sites
site_names <- sites
sites = as.numeric(as.factor(df$site)) # numeric site names
years <- as.numeric(as.factor(df$year)) - 1
n_years <- length(unique(years))

species_vector <- as.numeric((df$species_factor))

stan_data <- c("N", "y", "n_survey_events",
               "n_species", "X", "species_vector",
               "n_sites", "sites",
               "n_years", "years",
               "julian_scaled",
               "julian_scaled_sq"
)

# Parameters monitored
params <- c("beta",
            "beta_julian",
            "beta_julian_sq",
            "beta_site",
            "sigma_site",
            "beta_year",
            "sigma",
            "mean_y_rep_species",
            "max_y_rep_species"
)

# MCMC settings
n_iterations <- 4000
n_thin <- 1
n_burnin <- 0.5*n_iterations
n_chains <- 4
n_cores <- n_chains

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(beta0 = runif(1, -1, 1)
  )
)

## --------------------------------------------------
### Run model

stan_model <- "./phenology/models/final_model_phenology.stan"

## Call Stan from R
stan_out <- stan(stan_model,
                 data = stan_data, 
                 #init = inits, 
                 pars = params,
                 chains = n_chains, iter = n_iterations, 
                 warmup = n_burnin, thin = n_thin,
                 seed = 1,
                 open_progress = FALSE,
                 cores = n_cores,
                 control = list(adapt_delta = 0.9))

print(stan_out, digits = 3)

saveRDS(stan_out, "./phenology/model_outputs/stan_out_phenology.RDS")

stan_out <- readRDS("./phenology/model_outputs/stan_out_phenology.RDS")
fit_summary <- rstan::summary(stan_out)

## --------------------------------------------------
### Simple diagnostic plots

# traceplot
traceplot(stan_out, pars = c(
  "beta",
  "beta_julian",
  "beta_julian_sq"
))

# pairs plot
pairs(stan_out, pars = c(
  "beta"
))

pairs(stan_out, pars = c(
  "sigma"
))

## --------------------------------------------------
# posterior predictive check

fit_summary <- rstan::summary(stan_out)

View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")


# Summarize counts by species -> to get W_mean and W_max per species
W_species <- df %>%
  group_by(species) %>%
  mutate(W_mean_species = mean(n),
         W_max_species = max(n)) %>%
  slice(1) %>%
  ungroup() %>%
  select(species_factor, W_mean_species, W_max_species)

stan_fit_first_W_mean <- 32 # which row in the tracked parameters is the first W mean

df_estimates <- data.frame(X = numeric(), 
                           Y = numeric(), 
                           lower_95 = numeric(),
                           upper_95 = numeric(),
                           lower_50 = numeric(),
                           upper_50 = numeric()
) 

for(i in 1:n_species){
  
  row <- c((i - 1), 
           fit_summary$summary[(stan_fit_first_W_mean+(i-1)),1],
           fit_summary$summary[(stan_fit_first_W_mean+(i-1)),4],
           fit_summary$summary[(stan_fit_first_W_mean+(i-1)),8],
           fit_summary$summary[(stan_fit_first_W_mean+(i-1)),5],
           fit_summary$summary[(stan_fit_first_W_mean+(i-1)),7])
  
  df_estimates[i,] <- row
  
}

labels=species_names
ylims = c(0,(max(df_estimates$upper_95)+5))
end_point  = 0.5 + nrow(df_estimates) + nrow(df_estimates) - 1 #

par(mar = c(9,4,1,2))
plot(1, type="n",
     xlim=c(-0.5, n_species - 1 + 0.5), 
     xlab="",
     xaxt = "n",
     ylim=ylims, 
     ylab="50% and 95% Marginal Posterior Quantiles",
     main="Mean Detections vs. Model Expectations of Mean Detections")

#axis(1, at=start:(start+n), labels=labels, las = 2, cex.axis=.75)
text(seq(0, n_species - 1, by = 1), par("usr")[3]-0.25, 
     srt = 60, adj = 1, xpd = TRUE,
     labels = labels, cex = 1)

for(i in 1:n_species){
  sliced <- df_estimates[i,]
  W_sliced <- W_species[i, 2]
  
  rect(xleft = (sliced$X-0.35), xright=(sliced$X+0.35), 
       ytop = sliced$lower_95, ybottom = sliced$upper_95,
       col = c_mid, border = NA
  )
  
  rect(xleft =(sliced$X-0.35), xright=(sliced$X+0.35), 
       ytop = sliced$lower_50, ybottom = sliced$upper_50,
       col = c_mid_highlight, border = NA
  )
  
  points(x=sliced$X, y=W_sliced, pch=1)
}

#------------------------------------------------------------------------------
## PPC for max detections

stan_fit_first_W_max <- 37 # which row in the tracked parameters is the first W mean

df_estimates <- data.frame(X = numeric(), 
                           Y = numeric(), 
                           lower_95 = numeric(),
                           upper_95 = numeric(),
                           lower_50 = numeric(),
                           upper_50 = numeric()
) 

for(i in 1:n_species){
  
  row <- c((i - 1), 
           fit_summary$summary[(stan_fit_first_W_max+(i-1)),1],
           fit_summary$summary[(stan_fit_first_W_max+(i-1)),4],
           fit_summary$summary[(stan_fit_first_W_max+(i-1)),8],
           fit_summary$summary[(stan_fit_first_W_max+(i-1)),5],
           fit_summary$summary[(stan_fit_first_W_max+(i-1)),7])
  
  df_estimates[i,] <- row
  
}

labels=species_names
ylims = c(0, #max(W_species$W_max_species) + 5)
          (max(df_estimates$upper_95)+5))
end_point  = 0.5 + nrow(df_estimates) + nrow(df_estimates) - 1 #

par(mar = c(9,4,1,2))
plot(1, type="n",
     xlim=c(-0.5, n_species - 1 + 0.5), 
     xlab="",
     xaxt = "n",
     ylim=ylims, 
     ylab="50% and 95% Marginal Posterior Quantiles",
     main="Max Detections vs. Model Expectations of Max Detections")

#axis(1, at=start:(start+n), labels=labels, las = 2, cex.axis=.75)
text(seq(0, n_species - 1, by = 1), par("usr")[3]-0.25, 
     srt = 60, adj = 1, xpd = TRUE,
     labels = labels, cex = 1)

for(i in 1:n_species){
  sliced <- df_estimates[i,]
  W_sliced <- W_species[i, 3]
  
  rect(xleft = (sliced$X-0.35), xright=(sliced$X+0.35), 
       ytop = sliced$lower_95, ybottom = sliced$upper_95,
       col = c_mid, border = NA
  )
  
  rect(xleft =(sliced$X-0.35), xright=(sliced$X+0.35), 
       ytop = sliced$lower_50, ybottom = sliced$upper_50,
       col = c_mid_highlight, border = NA
  )
  
  points(x=sliced$X, y=W_sliced, pch=1)
}

