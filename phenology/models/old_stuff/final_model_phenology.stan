// STAN MODEL FOR BUMBLE BEE PHENOLOGY
// started June, 2024, J Ulrich

// Model is intended to estimate abundance by day of year for the different species

data {
  
  int<lower=0> N; // number of species*site*surveys
  
  int<lower=0> y[N]; // counts
  
  int<lower=1> n_sites;  // number of sites 
  int<lower=1, upper=n_sites> sites[N];  // vector of site names 
  
  int<lower=1> n_species;  // number of sites // number of level-3 clusters
  int<lower=1, upper=n_species> species[N];  // vector of site names // level-3 clusters
    
  int<lower=1> n_years; // number of years  
  int<lower=1, upper=n_years> years[N]; // vector of year names
  
  real julian_scaled[N]; // vector of julian dates
  real julian_scaled_sq[N]; // vector of julian dates squared
  
  int n_survey_events;
}

parameters {
  
  real beta0; // global intercept 
  
  // site random effect
  // site specific intercept allows some sites to have lower abundance than others, 
  // but with overall estimates for abundance partially informed by the data pooled across all sites.
  vector[n_sites] beta_site; // site specific intercept for count outcomes
  real<lower=0> sigma_site; // variance in site intercepts
  
  // species random effect
  // site specific intercept allows some sites to have lower abundance than others, 
  // but with overall estimates for abundance partially informed by the data pooled across all sites.
  vector[n_species] beta_species; // species specific intercept for count outcomes
  real<lower=0> sigma_species; // variance in site intercepts
  
  //vector[n_years] beta_year; // effect of year (2022 versus 2023) on overall abundance
  
  vector[n_species] beta_julian; // effect of julian date on abundance
  //real mu_julian;
  //real<lower=0> sigma_julian;
  vector<upper=0>[n_species] beta_julian_sq; // effect of julian date squared on abundance
  real mu_julian_sq;
  real<lower=0> sigma_julian_sq;
  
  vector[N] epsilon; // overdispersion parameter
  vector<lower=0>[n_species] sigma; // species-specific variation in overdispersion
}

transformed parameters{

  // the linear predictor for the individual observations
  // poisson process "intensity" that underlies the counts
  real lambda[N];

  // Individual flower mean
  for(i in 1:N){
    
      lambda[i] = exp( // exponential link function
             beta0 + // a global intercept
             beta_site[sites[i]] + // a site specific intercept
             beta_species[species[i]] + // a species specific intercept
             //beta_year[years[i]] + // a year specific intercept
             beta_julian[species[i]] * julian_scaled[i] + // an effect of julian date
             beta_julian_sq[species[i]] * julian_scaled_sq[i] + // an effect of julian date squared 
             epsilon[i] // plus a residual dispersion effect
             )
            ; // end lambda
  }
  
}

model {
  
  // PRIORS
  
  beta0 ~ normal(0, 2); // weakly informative prior for global intercept
  
  // site random effect
  beta_site ~ normal(0, sigma_site); // hierarchical prior, i.e., random effect
  sigma_site ~ normal(0, 2); // weakly informative prior
  
  // species effect
  beta_species ~ normal(0, sigma_species); // hierarchical prior, i.e., random effect
  sigma_species ~ normal(0, 2); // weakly informative prior

  // year effect
  //beta_year ~ normal(0, 2); // no hierarchical prior (only two groups; hard to estimate among group variance)
  
  // date effect
  beta_julian ~ normal(0, 1); // weakly informative prior for effect of date
  //mu_julian ~ normal(0, 2);
  //sigma_julian ~ normal(0, 0.5);
  beta_julian_sq ~ normal(mu_julian_sq, sigma_julian_sq); // weakly informative prior for effect ofdate squared
  mu_julian_sq ~ normal(0, 2);
  sigma_julian_sq ~ normal(0, 0.5);
  
  // residual effect
  for(i in 1:N){ // a residual dispersion for each count
    epsilon[i] ~ normal(0, sigma[species[i]]); // with the variance of the dispersion allowed to change among species
  }
  sigma ~ normal(0, 1); // with the variance for each species drawn from a half-normal distribution
  
  // LIKELIHOOD
  
  y ~ poisson(lambda); // y (count) is the outcome of a poisson process with intensity lambda
  
}

generated quantities {
  
  // posterior predictive check
  // group by species and predict the abundance
  // compare the mean abundance and max abundance in the predicted data
  // to the mean and max in our real data set
  // comparison is made post hoc, but we will generate the predictive distributions as such:s
  
  // define abundances to predict
  int<lower=0> y_rep[N]; // simulated abundance outcomes
  int<lower=0> y_rep_species[n_species, n_survey_events]; // abundance 
  // and abundance summary statistics to make
  real mean_y_rep_species[n_species];
  real max_y_rep_species[n_species];
  
  // Now generate posterior predictive distribution:
  // predict abundance for each speciesXsurvey event
  for(i in 1:N) { // loop across all possible events
    y_rep[i] = poisson_rng(lambda[i]);
  } 
  
  //fill in a matrix with rows for species and columns for possible survey events
  for(i in 1:n_species){
    for(j in 1:n_survey_events){
      y_rep_species[i,j] = y_rep[(j-1)*n_species + i];
    }
  }
  
  // calculate the MEAN abundance by species
  for(i in 1:n_species){
    mean_y_rep_species[i] = mean(y_rep_species[i,]);
  }
  
  // calculate the MAX abundance by species
  for(i in 1:n_species){
    max_y_rep_species[i] = max(y_rep_species[i,]);
  }
  
}
