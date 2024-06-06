// STAN MODEL FOR BUMBLE BEE PHENOLOGY
// started June, 2024, J Ulrich

// Model is intended to estimate abundance by day of year for the different species

data {
  
  int<lower=0> N; // number of species*site*surveys
  
  int<lower=0> y[N]; // counts
  
  int<lower=1> n_species;  // number of species
  matrix[N,n_species]  X;  // species matrix
  int species_vector[N]; // species vector
  
  int<lower=1> n_sites;  // number of sites 
  int<lower=1, upper=n_sites> sites[N];  // vector of site names 
    
  int<lower=1> n_years; // number of years  
  int<lower=0, upper=1> years[N]; // vector of year names
  
  real julian_scaled[N]; // vector of julian dates
  real julian_scaled_sq[N]; // vector of julian dates squared
  
  int n_survey_events;
}

parameters {
  
  vector[n_species] beta;
  
  // site random effect
  vector[n_sites] beta_site; // site specific intercept for count outcomes
  real<lower=0> sigma_site; // variance in site intercepts
  
  real beta_year; // effect of year (2022 versus 2023) on overall abundance
  
  vector[n_species] beta_julian;
  vector[n_species] beta_julian_sq;
  
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
             X[i] * beta + // species-specific intercepts
             beta_site[sites[i]] + // a site-effect on intercepts
             beta_year * years[i] + // an effect of 2023 versus 2022
             (X[i] * beta_julian) * julian_scaled[i] + // an effect of julian date
             (X[i] * beta_julian_sq) * julian_scaled_sq[i] + // an effect of julian date squared 
             epsilon[i] // plus a residual dispersion effect
             )
            ; // end lambda
  }
  
}

model {
  
  // PRIORS
  
  // species-specific intercept
  beta ~ normal(0, 2); // weakly informative prior for global intercept
  
  // site random effect
  beta_site ~ normal(0, sigma_site); // hierarchical prior, i.e., random effect
  sigma_site ~ normal(0, 2); // weakly informative prior
  
  // year effect
  beta_year ~ normal(0, 1);
  
  // date effect
  beta_julian ~ normal(0, 2); // weakly informative prior for effect of date
  
  beta_julian_sq[1] ~ normal(0, 2); // weakly informative prior for effect of date squared
  beta_julian_sq[2] ~ normal(0, 1); // weakly informative prior for effect of date squared
  beta_julian_sq[3] ~ normal(0, 0.25); // weakly informative prior for effect of date squared
  beta_julian_sq[4] ~ normal(0, 1); // weakly informative prior for effect of date squared
  beta_julian_sq[5] ~ normal(0, 1); // weakly informative prior for effect of date squared
  
  // residual effect
  for(i in 1:N){ // a residual dispersion for each count
    epsilon[i] ~ normal(0, sigma[species_vector[i]]); // with the variance of the dispersion allowed to change among species
  }
  sigma ~ normal(0, 0.5); // with the variance for each species drawn from a half-normal distribution
  
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
