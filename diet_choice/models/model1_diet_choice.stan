// STAN MODEL FOR BUMBLE BEE PHENOLOGY
// started June, 2024, J Ulrich

// Model is intended to estimate abundance by day of year for the different species

data {
  
  int<lower=0> N; // number of species*site*surveys
  
  int<lower=0> y[N]; // binary outcomes (1 = interact with invasive plant)
  
  int<lower=1> n_species;  // number of species
  matrix[N,n_species] X;  // species matrix
  
  int<lower=1> n_sites;  // number of sites 
  int<lower=1, upper=n_sites> sites[N];  // vector of site names 
  
  real prop_nvsv[N]; // proportion of plants that are invasive (on a given survey)
  
}

parameters {
  
  vector[n_species] beta; // species specific rate of interaction w invasive plants when
    // invasive plants are at an average proportion in the plant community
    
  vector[n_species] beta_prop_nvsv; // species specific effect of increasing prop. invasvive plants
  
  // site random effect
  vector[n_sites] beta_site; // site specific intercept for count outcomes
  real<lower=0> sigma_site; // variance in site intercepts

  
}

transformed parameters{

  // the linear predictor for the individual observations
  // that underlies the rate of choosing invasive plants
  real theta[N];

  // Individual flower mean
  for(i in 1:N){
      
      theta[i] = inv_logit( // exponential link function
             X[i] * beta + // species-specific intercepts
             beta_site[sites[i]] + // a site-effect on intercepts
             (X[i] * beta_prop_nvsv) * prop_nvsv[i] // a species-specific effect of increasing invasive plants
             )
            ; // end lambda
  }
  
}

model {
  
  // PRIORS
  
  // species-specific intercepts
  beta ~ normal(0, 2);
  
  // site random effect
  beta_site ~ normal(0, sigma_site); // hierarchical prior, i.e., random effect
  sigma_site ~ normal(0, 1); // weakly informative prior
  
  // species-specific effects of increasing prop. of invasive plants
  beta_prop_nvsv ~ normal(0, 2);
  
  // LIKELIHOOD
  y ~ bernoulli(theta); // y (count) is the outcome of a bernoulli trial w/ probability p

}

generated quantities {
  
  // posterior predictive check
  // group by site and predict the number of pollen limited plants
  // compare the number of pollen limited plants at each site to the
  // number of pollen limited plants in the real data
  
  // define abundances to predict
  int<lower=0> y_rep[N]; // simulated abundance outcomes
  int<lower=0> sum_y_rep; // simulated abundance outcomes
  
  // generating posterior predictive distribution
  // Predict invasive interaction for each instance,
  for(i in 1:N) { // loop across all plants
    y_rep[i] = bernoulli_rng(theta[i]);
  } 
  
  sum_y_rep = sum(y_rep);
  
}
