library(tidyverse)
library(rstan)

## predictor center scaling function
center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

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

#------------------------------------------------------------------------------
# GLM fit by maximum likelihood or with canned bayesian regression
# this formula is essentially the same as I specify in the custom model

fit2 <- glmer(invasive_interaction ~ 
                scaled_prop_nvsv + species + scaled_prop_nvsv*species + (1|site), 
              data = pref_data, family = binomial)

summary(fit2)

fit1 <- rstanarm::stan_glmer(data = pref_data, 
                             invasive_interaction ~ 
                               scaled_prop_nvsv + species + 
                               scaled_prop_nvsv*species + (1|site),
                                 family = "binomial")

summary(fit1)

rstanarm::pp_check(fit1)
rstanarm::pp_check(fit1, plotfun = "stat_grouped", group = "species")
rstanarm::pp_check(fit1, plotfun = "stat_grouped",  stat = "max", group = "species")

#------------------------------------------------------------------------------
# Prepare data for stan model

# data to feed to the model
# dummy_variables <- model.matrix(~ species_factor, data = df)
X <- model.matrix(invasive_interaction ~ species, data = pref_data)
n_species <- length(unique(pref_data$species))
sites <- as.numeric(as.factor(pref_data$site))
n_sites <- length(unique(pref_data$site))
N <- nrow(pref_data) # number of pairs
y <- pref_data$invasive_interaction # outcomes (counts)
prop_nvsv <- pref_data$scaled_prop_nvsv

stan_data <- c("N", "y", 
               "X", "n_species",
               "sites", "n_sites",
               "prop_nvsv"
)

# Parameters monitored
params <- c("beta",
            "beta_prop_nvsv",
            "beta_site",
            "sigma_site",
            "sum_y_rep"
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
  
  list(beta = runif(1, -1, 1)
  )
)

## --------------------------------------------------
### Run model

stan_model <- "./diet_choice/models/model1_diet_choice.stan"

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

saveRDS(stan_out, "./diet_choice/model_outputs/stan_out_diet_choice.RDS")
stan_out <- readRDS("./diet_choice/model_outputs/stan_out_diet_choice.RDS")
fit_summary <- rstan::summary(stan_out)

## --------------------------------------------------
### Simple diagnostic plots

# traceplot
traceplot(stan_out, pars = c(
  "beta",
  "beta_prop_nvsv",
  "sigma_site"
))

# pairs plot
pairs(stan_out, pars = c(
  "beta"
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

(W_overall <- sum(pref_data$invasive_interaction))

stan_fit_first_W_mean <- 21 # which row in the tracked parameters is the first W mean

df_estimates <- data.frame(X = numeric(), 
                           Y = numeric(), 
                           lower_95 = numeric(),
                           upper_95 = numeric(),
                           lower_50 = numeric(),
                           upper_50 = numeric()
) 

for(i in 1:length(W_overall)){
  
  row <- c((i - 1), 
           fit_summary$summary[(stan_fit_first_W_mean+(i-1)),1],
           fit_summary$summary[(stan_fit_first_W_mean+(i-1)),4],
           fit_summary$summary[(stan_fit_first_W_mean+(i-1)),8],
           fit_summary$summary[(stan_fit_first_W_mean+(i-1)),5],
           fit_summary$summary[(stan_fit_first_W_mean+(i-1)),7])
  
  df_estimates[i,] <- row
  
}

labels=""
ylims = c((min(df_estimates$lower_95)-100),(max(df_estimates$upper_95)+100))
end_point  = 0.5 + nrow(df_estimates) + nrow(df_estimates) - 1 #

par(mar = c(6,4,3,2))
plot(1, type="n",
     xlim=c(-0.5, 0.5), 
     xlab="",
     xaxt = "n",
     ylim=ylims, 
     ylab="50% and 95% Marginal Posterior Quantiles",
     main="Interactions with invasive plants vs.\nmodel expectations of interactions with inv. plants")

#axis(1, at=start:(start+n), labels=labels, las = 2, cex.axis=.75)
text(seq(0, 1 - 1, by = 1), par("usr")[3]-0.25, 
     srt = 60, adj = 1, xpd = TRUE,
     labels = labels, cex = 1)

for(i in 1:1){
  sliced <- df_estimates[i,]
  W_sliced <- W_overall
  
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

