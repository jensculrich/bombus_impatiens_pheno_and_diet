library(rstan)
library(tidyverse)

## predictor center scaling function
center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

#-----------------------------------------------------------------------
# To visualize species-specific phenology differences:
# Plot bumble bee counts (by species) across date with the upper and lower 95 BCI. 
# We could also make a caterpillar plot with all parameter estimates.

# First we need our original data to get the prediction interval
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

# IMPORTANT: enter names in alphabetical order
species_names <- c("bom_fla", "bom_imp", "bom_mel", "bom_mix", "bom_vos")
# how many species are we going to consider? Just the main 5 for now
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
         #species_factor = relevel(factor(species), ref = 'bom_imp'))
         species_factor = as.factor(species))



#-----------------------------------------------------------------------
# Use the estimates from the fitted poisson model with overdispersion in counts
# read in the saved RDS file containing the model outputs
fit_summary <- rstan::summary(readRDS("./phenology/model_outputs/stan_out_phenology.RDS"))
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

# get original julian dates back
pred_length <- 100 # divide the calendar year by some number of days
lower_date <- 70
upper_date <- 330
original_julian <- seq(lower_date, upper_date, length.out = pred_length) # reasonable prediction range
julian_pred <- (original_julian - mean(df$julian)) / sd(df$julian) # unscale the dates

count_mean <- array(dim = c(n_species, pred_length))
count_lower95 <- array(dim = c(n_species, pred_length))
count_upper95 <- array(dim = c(n_species, pred_length))
count_lower50 <- array(dim = c(n_species, pred_length))
count_upper50 <- array(dim = c(n_species, pred_length))

# expected counts
for(i in 1:n_species){
  for(j in 1:pred_length){
    
    if(i != 1){ # i != b. impatiens, ie. the reference level
      # compute for a comparison species
      count_mean[i,j] =
        exp(
          # intercept +
          fit_summary$summary[1,1] +
            # species differences in overall abundance
            fit_summary$summary[1+(i-1),1] + 
            # julian scaled with species effect
            (fit_summary$summary[6,1] + fit_summary$summary[6+(i-1),1]) * julian_pred[j] +
            # julian scaled squared with species effect
            (fit_summary$summary[11,1] + fit_summary$summary[11+(i-1),1]) * julian_pred[j]^2 #+
        )
      
      count_lower95[i,j] =
        exp(
          # intercept +
          fit_summary$summary[1,4] +
            # species differences in overall abundance
            fit_summary$summary[1+(i-1),4] + 
            # julian scaled with species effect
            (fit_summary$summary[6,4] + fit_summary$summary[6+(i-1),4]) * julian_pred[j] +
            # julian scaled squared with species effect
            (fit_summary$summary[11,4] + fit_summary$summary[11+(i-1),4]) * julian_pred[j]^2 #+
        )
      
      count_upper95[i,j] =
        exp(
          # intercept +
          fit_summary$summary[1,8] +
            # species differences in overall abundance
            fit_summary$summary[1+(i-1),8] + 
            # julian scaled with species effect
            (fit_summary$summary[6,8] + fit_summary$summary[6+(i-1),8]) * julian_pred[j] +
            # julian scaled squared with species effect
            (fit_summary$summary[11,8] + fit_summary$summary[11+(i-1),8]) * julian_pred[j]^2 #+
        )
      
      count_lower50[i,j] =
        exp(
          # intercept +
          fit_summary$summary[1,5] +
            # species differences in overall abundance
            fit_summary$summary[1+(i-1),5] + 
            # julian scaled with species effect
            (fit_summary$summary[6,5] + fit_summary$summary[6+(i-1),5]) * julian_pred[j] +
            # julian scaled squared with species effect
            (fit_summary$summary[11,5] + fit_summary$summary[11+(i-1),5]) * julian_pred[j]^2 #+
        )
      
      count_upper50[i,j] =
        exp(
          # intercept +
          fit_summary$summary[1,7] +
            # species differences in overall abundance
            fit_summary$summary[1+(i-1),7] + 
            # julian scaled with species effect
            (fit_summary$summary[6,7] + fit_summary$summary[6+(i-1),7]) * julian_pred[j] +
            # julian scaled squared with species effect
            (fit_summary$summary[11,7] + fit_summary$summary[11+(i-1),7]) * julian_pred[j]^2 #+
        )
      
    } else { 
      
      # else compute for the reference level
      count_mean[i,j] =
        exp(
          # intercept +
          fit_summary$summary[1,1] +
            # julian scaled
            fit_summary$summary[6,1] * julian_pred[j] +
            # julian scaled squared
            fit_summary$summary[11,1] * julian_pred[j]^2 
        )
      
      count_lower95[i,j] =
        exp(
          # intercept +
          fit_summary$summary[1,4] +
            # julian scaled
            fit_summary$summary[6,4] * julian_pred[j] +
            # julian scaled squared
            fit_summary$summary[11,4] * julian_pred[j]^2 
        )
      
      count_upper95[i,j] =
        exp(
          # intercept +
          fit_summary$summary[1,8] +
            # julian scaled
            fit_summary$summary[6,8] * julian_pred[j] +
            # julian scaled squared
            fit_summary$summary[11,8] * julian_pred[j]^2 
        )
      
      count_lower50[i,j] =
        exp(
          # intercept +
          fit_summary$summary[1,5] +
            # julian scaled
            fit_summary$summary[6,5] * julian_pred[j] +
            # julian scaled squared
            fit_summary$summary[11,5] * julian_pred[j]^2 
        )
      
      count_upper50[i,j] =
        exp(
          # intercept +
          fit_summary$summary[1,7] +
            # julian scaled
            fit_summary$summary[6,7] * julian_pred[j] +
            # julian scaled squared
            fit_summary$summary[11,7] * julian_pred[j]^2 
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
lower_50 <- vector()
upper_50 <- vector()

# fill predictions by species and add sequentially to the vector
for(i in 1:n_species){
  mean <- c(mean, count_mean[i,]) 
  lower_95 <- c(lower_95, count_lower95[i,])
  upper_95 <- c(upper_95, count_upper95[i,])
  lower_50 <- c(lower_50, count_lower50[i,])
  upper_50 <- c(upper_50, count_upper50[i,])
}
# unite as a dataframe and rename the repped julian dates
new_df <- as.data.frame(cbind(species, julian_pred, rep(original_julian, n_species),
                              mean, lower_95, upper_95, lower_50, upper_50)) %>%
  rename_with(.cols = 3, ~"julian") %>%
  mutate(species_new = as.factor(species)) %>%
  mutate(species_new_ordered = factor(species_new, 
                                      levels = c("1","4", "2", "5", "3")))

library(RColorBrewer)

#my_palette <- brewer.pal(5, "Blues")
#my_palette[1] <- "#7C0000"
  
my_palette <- hcl.colors(5, "Zissou 1")

# plot just the expected means (plot p)
p <- ggplot(data = new_df, aes(julian, mean, fill=as.factor(species_new_ordered))) +
  
  #geom_ribbon(aes(
    #ymin=lower_50, ymax=upper_50), alpha=0.5) +
  geom_ribbon(aes(
    ymin=lower_95, ymax=upper_95), alpha=0.5) +
  #geom_line(size=2) +
  xlim(c(lower_date, upper_date)) +
  theme_bw() +
  ylab("Expected abundance") +
  xlab("Julian date") +
  scale_fill_manual(name = "Species",
                    labels=c("B. impatiens",
                             "B. mixtus",
                             "B. flavifrons",
                             "B. vosnesenskii", 
                             "B. melanopygus"),
                    values = my_palette) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0, 5, 10)) +
  theme(legend.position = c(0.15, 0.7),
        legend.title=element_text(size=16),
        legend.text=element_text(size=14),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16, angle=45, vjust=-0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

#-------------------------------------------------------------------------------
# make caterpillar plots
## --------------------------------------------------
## Plot ecological paramter means and variation

# parameter means
# let's plot species effects of julian date (for 5 species)
X_eco <- c(1, 2, 3, 4) # params of interest

# mean of eco params
Y_eco <- c(
           fit_summary$summary[7,1],
           fit_summary$summary[8,1],
           fit_summary$summary[9,1], 
           fit_summary$summary[10,1]
)

# confidence intervals
lower_95_eco <- c( 
                  fit_summary$summary[7,4],
                  fit_summary$summary[8,4],
                  fit_summary$summary[9,4], 
                  fit_summary$summary[10,4]
)

upper_95_eco <- c(
                  fit_summary$summary[7,8],
                  fit_summary$summary[8,8],
                  fit_summary$summary[9,8], 
                  fit_summary$summary[10,8]
)

# confidence intervals
lower_50_eco <- c(
                  fit_summary$summary[7,5],
                  fit_summary$summary[8,5],
                  fit_summary$summary[9,5], 
                  fit_summary$summary[10,5]
)

upper_50_eco <- c(
                  fit_summary$summary[7,7],
                  fit_summary$summary[8,7],
                  fit_summary$summary[9,7], 
                  fit_summary$summary[10,7]
)


df_estimates_eco <- as.data.frame(cbind(X_eco, Y_eco, 
                                        lower_95_eco, upper_95_eco,
                                        lower_50_eco, upper_50_eco))

df_estimates_eco$X_eco <- as.factor(df_estimates_eco$X_eco)

## --------------------------------------------------
## Draw ecological parameter plot

(s <- ggplot(df_estimates_eco) +
   theme_bw() +
   # scale_color_viridis(discrete=TRUE) +
   scale_x_discrete(name="", breaks = c(1, 2, 3, 4),
                    labels=c(
                             "B. flavifrons",
                             "B. melanopygus",
                             "B. mixtus", 
                             "B. vosnesenskii")
                    ) +
   scale_y_continuous("Effect of species identity on phenology peak\n(relative to B. impatiens) (log-scaled)",
                      limits = c(-10, 0), breaks = c(-10, -5, 0)) +
   guides(color = guide_legend(title = "")) +
   geom_hline(yintercept = 0, lty = "dashed") +
   theme(legend.text=element_text(size=10),
         axis.text.x = element_text(size = 16),
         axis.text.y = element_text(size = 16, angle=0, vjust=0.5),
         axis.title.x = element_text(size = 16),
         axis.title.y = element_text(size = 16),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   coord_flip()
)

s <- s +
  geom_errorbar(aes(x=X_eco, ymin=lower_95_eco, ymax=upper_95_eco),
                color="black",width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X_eco, ymin=lower_50_eco, ymax=upper_50_eco),
                color="black",width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X_eco, y=Y_eco),
             size = 5, alpha = 0.8) 


s

#gridExtra::grid.arrange(p, s, ncol=1)
cowplot::plot_grid(p, s, 
                   labels = c("a)", "b)"), label_size = 16,
                   hjust = 0,
                   nrow = 2, rel_heights = c(2/3, 1/3))

#-------------------------------------------------------------------------------
## make nice tables of param estimates

phenology_table_data1 <- data.frame(
  Term = c("(Intercept)", "B. flavifrons", "B. mixtus", "B. melanopygus", "B. vosnesenskii",
           "(Sigma)", "B. flavifrons", "B. mixtus", "B. melanopygus", "B. vosnesenskii", "Year\n(2023 v. 2022)"),
  'mean' = c(fit_summary$summary[1,1], fit_summary$summary[2,1], 
             fit_summary$summary[3,1], fit_summary$summary[4,1], fit_summary$summary[5,1],
             fit_summary$summary[27,1], fit_summary$summary[28,1], 
             fit_summary$summary[29,1], fit_summary$summary[30,1], fit_summary$summary[31,1],
             fit_summary$summary[26,1]),
  'lower 95 BCI' = c(fit_summary$summary[1,4], fit_summary$summary[2,4], 
                     fit_summary$summary[3,4], fit_summary$summary[4,4], fit_summary$summary[5,4],
                     fit_summary$summary[27,4], fit_summary$summary[28,4], 
                     fit_summary$summary[29,4], fit_summary$summary[30,4], fit_summary$summary[31,4],
                     fit_summary$summary[26,4]),
  'upper 95 BCI' = c(fit_summary$summary[1,8], fit_summary$summary[2,8], 
                     fit_summary$summary[3,8], fit_summary$summary[4,8], fit_summary$summary[5,8],
                     fit_summary$summary[27,8], fit_summary$summary[28,8], 
                     fit_summary$summary[29,8], fit_summary$summary[30,8], fit_summary$summary[31,8],
                     fit_summary$summary[26,8])
)
print(phenology_table_data1)
nice_table(phenology_table_data1)


phenology_table_data2 <- data.frame(
  Term = c("(Julian date)", "B. flavifrons", "B. mixtus", "B. melanopygus", "B. vosnesenskii",
           "(Julian date)^2", "B. flavifrons", "B. mixtus", "B. melanopygus", "B. vosnesenskii"),
  'mean' = c(fit_summary$summary[6,1], fit_summary$summary[7,1], 
             fit_summary$summary[8,1], fit_summary$summary[9,1], fit_summary$summary[10,1],
             fit_summary$summary[11,1], fit_summary$summary[12,1], 
             fit_summary$summary[13,1], fit_summary$summary[14,1], fit_summary$summary[15,1]),
  'lower 95 BCI' = c(fit_summary$summary[6,4], fit_summary$summary[7,4], 
                     fit_summary$summary[8,4], fit_summary$summary[9,4], fit_summary$summary[10,4],
                     fit_summary$summary[11,4], fit_summary$summary[12,4], 
                     fit_summary$summary[13,4], fit_summary$summary[14,4], fit_summary$summary[15,4]),
  'upper 95 BCI' = c(fit_summary$summary[6,8], fit_summary$summary[7,8], 
                     fit_summary$summary[8,8], fit_summary$summary[9,8], fit_summary$summary[10,8],
                     fit_summary$summary[11,8], fit_summary$summary[12,8], 
                     fit_summary$summary[13,8], fit_summary$summary[14,8], fit_summary$summary[15,8])
)
print(phenology_table_data2)
nice_table(phenology_table_data2)
