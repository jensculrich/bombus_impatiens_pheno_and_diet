// STAN MODEL FOR BUMBLE BEE PHENOLOGY
// started June, 2024, J Ulrich

// Model is intended to estimate abundance by day of year for the different species

data {
  
  int<lower=0> N; // number of species*site*surveys
  
  int<lower=0> y[N]; // counts
  
  //int<lower=1> n_species;  // number of species
  //matrix[N,n_species]  X;  // species matrix
  
  real prop_nvsv[N];
  
}

parameters {
  
  real beta0;
  //vector[n_species] beta;
  
  real beta_prop_nvsv;
  
}

transformed parameters{

  // the linear predictor for the individual observations
  // poisson process "intensity" that underlies the counts
  real p[N];

  // Individual flower mean
  for(i in 1:N){
      
      p[i] = inv_logit( // exponential link function
             beta0 +
             //X[i] * beta // species-specific intercepts
             beta_prop_nvsv * prop_nvsv[i]
             )
            ; // end lambda
  }
  
}

model {
  
  // PRIORS
  
  // species-specific intercept
  beta0 ~ normal(0, 2);
  //beta ~ normal(0, 2); // weakly informative prior for global intercept
  
  beta_prop_nvsv ~ normal(0, 2);
  
  // LIKELIHOOD
  y ~ bernoulli_logit(p); // y (count) is the outcome of a bernoulli trial w/ probability p

}

generated quantities {
  
  // posterior predictive check
  // group by site and predict the number of pollen limited plants
  // compare the number of pollen limited plants at each site to the
  // number of pollen limited plants in the real data
  
  // define abundances to predict
  int<lower=0> y_rep[N]; // simulated abundance outcomes
  
  // generating posterior predictive distribution
  // Predict invasive interaction for each instance,
  for(i in 1:N) { // loop across all plants
    y_rep[i] = bernoulli_logit_rng(p[i]);
  } 
  
  sum_y_rep = sum(y_rep);
  
}
