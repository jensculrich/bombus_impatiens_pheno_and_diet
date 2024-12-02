# *Bombus impatiens* phenology and diet project

Data and code for the ms titled "", by Melissa Platsko, Jens Ulrich and Risa Sargent.

This repository including the analysis instructions and metadata file prepared by J. Ulrich in December 2024.

### Phenology (Figure 1)
**To run the phenology analysis**, navigate to ./phenology/run_model/run_model_phenology.R

With this analysis, we used a Bayesian GLM (overdispersed Poisson distribution) to test whether species occur at different abundances ("beta"), whether they have different peak day of year of abundance ("beta_julian"), and/or have different width of flight season ("beta_julian_sq"). We used two years of data and included a control covariate of year ("beta_year") and also include site as a random effect ("beta_site").  
 
This file will prep the "./data/melissa_2022_2023_phenology_data.csv" file for analysis. See metadata for more information about the phenology data file. The run model file will prep a model covariate matrix that places B. impatiens as the reference level. With this framework, each species will get it's own abundance intercept ("beta") and species-specific covariate effects. However, the beta[1] is the intercept for B. impatiens. All other intercept terms are the difference that needs to be added to the B. impatiens reference level. E.g., species number 2 is B. flavifrons. So the abundance intercept for B. flavifrons is beta[1] + beta[2]. Species number 3 is bombus melanopygus. The abundance intercept for B. melanopygus is beta[1] + beta[3]. Etc. etc.
 
To run the model, enter all lines of code, finally ending with specifying the final model we used for the analysis
```stan_model <- "./phenology/models/final_model_phenology.stan"```

and then calling stan
```stan_out <- stan()```

The code to visualize model diagnostics (traceplots, pairs plots) and conduct goodness of fit tests (posterior predictive checks) are included in the run_model_phenology.R file, at the bottom of the script after the code used to run the model.

To make **Figure 1** (the results of the phenology analysis), navigate to ./phenology/figures/make_figures.R

### Diet Choice (Figure 2)
**To run the diet choice analysis**, navigate to ./diet_choice/run_model/run_model_diet_choice.R
 
With this analysis, we used a Bayesian GLM (Binomial distribution, i.e., logistic regression) to test whether species have different odds of choosing an invasive plant when invasive plants are at an average abundance ("beta"), and whether there are species-specific differences in how this rate changes as invasive plants become more common in the community ("beta_prop_nvsv"). We also included controlling random effect of site ("beta_site").
 
This file will prep the "./data/melissa_2022_2023_logistic_analyses_file.csv" file for analysis. See metadata for more information about this diet choice data file. The run model file will prep a model covariate matrix that places B. impatiens as the reference level. With this framework, each species will get it's own rate intercept ("beta") and species-specific covariate effects. However, the beta[1] is the intercept for B. impatiens. All other intercept terms are the difference that needs to be added to the B. impatiens reference level. E.g., species number 2 is B. flavifrons. So the rate intercept for B. flavifrons is beta[1] + beta[2]. Species number 3 is bombus melanopygus. The rate intercept for B. melanopygus is beta[1] + beta[3]. Etc. etc.
 
To run the model, enter all lines of code, finally ending with specifying the final model we used for the analysis
```stan_model <- "./phenology/models/model1_diet_choice.stan"```

and then calling stan
```stan_out <- stan()```

The code to visualize model diagnostics (traceplots, pairs plots) and conduct goodness of fit tests (posterior predictive checks) are included in the run_model_diet_choice.R file, at the bottom of the script after the code used to run the model.

To make **Figure 2** (the results of the diet choice analysis), navigate to ./diet_choice/figures/make_figures.R

### Plant-pollinator interaction NMDS (Figure 3)
**To run the interaction NMDS analysis**, navigate to ./nmds/Melissa_NMDS_code_subset.R. This file will construct an NMDS just for our 5 focal bumble bee species.
 
This file will use the R package "*vegan*" to construct an NMDS from the ./nmds/data/Melissa_commun_matrix_no0.csv" interaction data file. See metadata for more information about this diet choice data file. The Bray-Curtis dissimilarity test indicates that there were no significant differences among species. In other words, species in our bumble bee community did not have distinct clusters of plant species that they uniquely interacted with.

The code used to make **Figure 3** (the results of the NMDS analysis) are included in the ./nmds/Melissa_NMDS_code_subset.R file. 

### Metadata
View the "metadata.docx" file in the main directory to review the explanations of the data used for the analyses. The data used for the three core analyses are:
- ./data/melissa_2022_2023_phenology_data.csv
- ./data/melissa_2022_2023_logistic_analyses_file.csv
- ./nmds/data/ (multiple files)