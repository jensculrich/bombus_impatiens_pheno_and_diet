// STAN MODEL FOR BUMBLE BEE PHENOLOGY
// started Dec. 14, 2022, J Ulrich

// Model is intended to estimate abundance by day of year for the different species

data {
  
  int<lower=0> N; // number of species*site*surveys
  
  int<lower=0> y[N]; // counts
  
  int<lower=1> n_sites;  // number of sites 
  int<lower=1, upper=n_sites> sites[N];  // vector of site names 
  
  int<lower=1> n_species;  // number of sites // number of level-3 clusters
  int<lower=1, upper=n_species> species[N];  // vector of site names // level-3 clusters
  
  real julian_scaled[N]; // vector of julian dates
  real julian_scaled_sq[N]; // vector of julian dates squared
  
  int n_survey_events;
}

parameters {
  
  real beta0; // global intercept 
  
  // site random effect
  // site specific intercept allows some sites to have lower success than others, 
  // but with overall estimates for success partially informed by the data pooled across all sites.
  real beta_site[n_sites]; // site specific intercept for count outcomes
  real<lower=0> sigma_site; // variance in site intercepts
  
  vector[n_species] beta_species;
  
  //real beta_julian; // effect of julian date on abundance
  //real beta_julian_sq; // effect of julian date squeared on abundance
  
  vector[n_species] beta_julian; // effect of julian date on abundance
  vector[n_species] beta_julian_sq; // effect of julian date squeared on abundance
  //real mu_beta_julian;
  //real<lower=0> sigma_beta_julian;
  //real mu_beta_julian_sq;
  //real<lower=0> sigma_beta_julian_sq;
  
  real<lower=0, upper=1> dispersion;
  
}

transformed parameters{

  // the linear predictor for the individual observations
  // prob. of being pollen limited
  real lambda[N];

  // Individual flower mean
  for(i in 1:N){
    
      lambda[i] = exp( // exponential link function
             beta0 + // a global intercept
             beta_site[sites[i]] + // a site specific intercept
             beta_species[species[i]] + // a species specific intercept
             //beta_julian * julian_scaled[i] + // an effect of julian date
             //beta_julian_sq * julian_scaled_sq[i] // an effect of julian date squared 
             beta_julian[species[i]] * julian_scaled[i] + // an effect of julian date
             beta_julian_sq[species[i]] * julian_scaled_sq[i] // an effect of julian date squared 
             )
            ; // end lambda
  }
  
}

model {
  
  // PRIORS
  
  beta0 ~ normal(0, 2); // weakly informative prior for global intercept
  
  // site random effect
  beta_site ~ normal(0, sigma_site); 
  // prob of success intercept for each site drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_site ~ normal(0, 1); // weakly informative prior
  
  // species effect
  beta_species ~ normal(0, 2);
  
  // date effect
  //beta_julian ~ normal(0, 2); // weakly informative prior for effect of date
  //beta_julian_sq ~ normal(0, 2); // weakly informative prior for effect ofdate squared
  
  beta_julian ~ normal(0, 2); // weakly informative prior for effect of date
  //mu_beta_julian ~ normal(0, 2);
  //sigma_beta_julian ~ normal(0, 2);
  beta_julian_sq ~ normal(0, 2); // weakly informative prior for effect ofdate squared
  //mu_beta_julian ~ normal(0, 0.5);
  //sigma_beta_julian ~ normal(0, 0.5);
  
  dispersion ~ beta(2, 2);
  
  // LIKELIHOOD
  
  //y ~ poisson(lambda);
  y ~ neg_binomial(lambda, dispersion);
}

generated quantities {
  
  // posterior predictive check
  // group by site and predict the number of pollen limited plants
  // compare the number of pollen limited plants at each site to the
  // number of pollen limited plants in the real data
  
  int<lower=0> y_rep[N]; // simulated abundance outcomes
  int<lower=0> y_rep_species[n_species, n_survey_events]; // abundance 
  real mean_y_rep_species[n_species];
  real max_y_rep_species[n_species];
  
  // generating posterior predictive distribution
  // Predict abundance for each species survey event
  for(i in 1:N) { // loop across all possible events
    //y_rep[i] = poisson_rng(lambda[i]);
    y_rep[i] = neg_binomial_rng(lambda[i], dispersion);
  } 
  
  // now fill in a matrix with rows for species and columns for possible survey events
  for(i in 1:n_species){
    for(j in 1:n_survey_events){
      y_rep_species[i,j] = y_rep[(j-1)*n_species + i];
    }
  }
  
  for(i in 1:n_species){
    mean_y_rep_species[i] = mean(y_rep_species[i,]);
  }
  
  for(i in 1:n_species){
    max_y_rep_species[i] = max(y_rep_species[i,]);
  }
  
}
