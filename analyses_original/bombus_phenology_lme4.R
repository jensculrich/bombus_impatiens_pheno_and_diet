# bombus phenology analysis and plot for melissa's data
# using lme4 to estimate parameter point estimates with frequentist implementation

library(tidyverse)
library(lme4)
library(RColorBrewer)

## predictor center scaling function
center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

#------------------------------------------------------------------------------
# initial data cleanup

df_orig <- as.data.frame(read.csv("./data/melissa_2022_2023_phenology_data.csv"))

table(df_orig$species)
table(df_orig$plant)

# see how many detections per species
unique_detections_by_species <- df_orig %>%
  group_by(date, site, species) %>%
  slice(1) %>%
  ungroup() %>%
  group_by(species) %>%
  add_tally() %>%
  slice(1) %>%
  select(species, n)

# how many species are we going to consider? Just the main 4 for now
n_species = 5
species_names <- c("bom_fla", "bom_mel", "bom_mix", "bom_vos", "bom_imp")

df_subset_with_counts <- subset(df_orig, species %in% species_names) %>%
  # let's try looking at just the workers
  filter(caste == "w") %>%
  group_by(date, site, species) %>%
  add_tally() %>%
  slice(1) %>%
  ungroup() %>%
  select(date, year, site, species, n) 

# format data types as factors / scale the dates
df_site_visits <- df_orig %>%
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
# refit a model to the data

# poisson model
fit1 <- lme4::glmer( 
  formula = 
    n ~ (1|site_factor) + year_factor +
    species_factor +
    julian_scaled + julian_scaled_sq + 
    julian_scaled*species_factor +
    julian_scaled_sq*species_factor,
  data = df,
  family = "poisson")

summary(fit1)
library(car)
fit1_anova <- Anova(fit1, type = 3)
fit1_anova

fit0 <- lme4::glmer( 
  formula = 
    n ~ (1|site_factor) + year_factor +
    species_factor +
    julian_scaled + julian_scaled_sq ,
  data = df,
  family = "poisson")

Anova(model_name,type=3)

summary(fit0)

anova(fit0, fit1)

#------------------------------------------------------------------------------
# Create tables with model data

library(rempsyc)

phen_table_data <- data.frame(
  Term = c("(Intercept)", "year_factor2023", "julian_scaled", "julian_scaled_sq"),
  'b' = c(2.05362, 0.40912, 1.89563, -3.05580),
  'SE' = c(0.26602, 0.04119, 0.10885, 0.12783),
  'z' = c(7.720, 9.933, 17.416, -23.904),
  'p' = c(0.001, 0.001, 0.001, 0.001)
)
print(phen_table_data)
phen_table_data$p <- as.numeric(gsub("[^0-9.]+", "", phen_table_data$p))
nice_table(phen_table_data)


#------------------------------------------------------------------------------
# plot predictions for a fit
# let's use fit1 for now (poisson model with zero counts)

# enter the name of the fit you want to use here
fit <- fit1

se <- sqrt(diag(vcov(fit)))
# table of estimates with 95% CI
(fit_summary <- cbind(Est = fixef(fit), LL = fixef(fit) - 1.96 * se, UL = fixef(fit) + 1.96 *
                se))

# get original julian dates back
pred_length <- 100 # divide the calendar year by some number of days
lower_date <- 80
upper_date <- 325
original_julian <- seq(80, 325, length.out = pred_length) # reasonable prediction range
julian_pred <- (original_julian - mean(df$julian)) / sd(df$julian) # unscale the dates

#orig_jul_extended <- seq(70, 365, length.out = pred_length)
#jul_pred_extended <- (orig_jul_extended - mean(df$julian)) / sd(df$julian)

# prep arrays to fill with predictions
count_mean <- array(dim = c(n_species, pred_length))
count_lower95 <- array(dim = c(n_species, pred_length))
count_upper95 <- array(dim = c(n_species, pred_length))

## now predict abundance across the range of dates using the model fit
for(i in 1:n_species){
  for(j in 1:pred_length){
    
    if(i != 5){
      
      # compute expected mean abundance for a comparison species
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
      
      # compute lower 95 for the expected mean abundance
      count_lower95[i,j] =
        exp(
          # intercept +
          fit_summary[1,2] +
            # year * 0 for 2022 +
            fit_summary[2,2] * 0 + 
            # species differences in overall abundance
            fit_summary[3+(i-1),2] + 
            # julian scaled with species effect
            (fit_summary[7,2] + fit_summary[9+(i-1),2]) * julian_pred[j] +
            # julian scaled squared with species effect
            (fit_summary[8,2] + fit_summary[13+(i-1),2]) * julian_pred[j]^2 #+
        )
      
      # compute upper 95 for the expected mean abundance
      count_upper95[i,j] =
        exp(
          # intercept +
          fit_summary[1,3] +
            # year * 0 for 2022 +
            fit_summary[2,3] * 0 + 
            # species differences in overall abundance
            fit_summary[3+(i-1),3] + 
            # julian scaled with species effect
            (fit_summary[7,3] + fit_summary[9+(i-1),3]) * julian_pred[j] +
            # julian scaled squared with species effect
            (fit_summary[8,3] + fit_summary[13+(i-1),3]) * julian_pred[j]^2 #+
        )
      
    } else { 
      # there are no species effects for species 4 (B. impatiens)
      # because it was used as the reference level.
      # so we have to do the predictions WITHOUT any species effects.
      
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
          fit_summary[1,2] +
            # year * 0 for 2022 +
            fit_summary[2,2] * 0 + 
            # julian scaled
            fit_summary[7,2] * julian_pred[j] +
            # julian scaled squared
            fit_summary[8,2] * julian_pred[j]^2 
        )
      
      count_upper95[i,j] =
        exp(
          # intercept +
          fit_summary[1,3] +
            # year * 0 for 2022 +
            fit_summary[2,3] * 0 + 
            # julian scaled
            fit_summary[7,3] * julian_pred[j] +
            # julian scaled squared
            fit_summary[8,3] * julian_pred[j]^2 
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


#------------------------------------------------------------------------------
# Now we are ready to draw the plot

# choose a colour palette for the plot, I like "Blues"
my_palette_reduced <- brewer.pal(n_species, "Paired")

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

## plot the expected means AND the actual data

# need to get survey data 
df <- df %>%
  mutate(
         mean = as.numeric(n), # easiest if counts have same column name
         # need same factor level order used for the pred lines and ribbons
         species_factor = factor(species, levels = c(
           'bom_fla', 'bom_mel', 'bom_mix',
           'bom_vos', 'bom_imp'))
         ) 

# plot the expected means and the survey counts (plot q)
q <- ggplot(data = new_df, aes(julian, mean)) +
  geom_line(aes(group=as.factor(species)), size=2) +
  geom_ribbon(aes(ymin=lower_95, ymax=upper_95, 
                  fill = as.factor(species)), alpha=0.8) +
  geom_point(data = df, 
             aes(x=julian, y=mean, colour=species_factor)) +
  xlim(c(100, 300)) +
  ylim(c(0, 105)) +
  theme_bw() +
  ylab("abundance per survey \n(with 95% CI)") +
  xlab("julian date") +
  scale_fill_manual(name = "Species",
                    labels=c("B. flavifrons",
                             "B. melanopygus",
                             "B. mixtus",
                             "B. vosnesenskii", 
                             "B. impatiens"),
                    values=my_palette_reduced) +
  scale_colour_manual(name = "Species",
                      labels=c("B. flavifrons",
                               "B. melanopygus",
                               "B. mixtus",
                               "B. vosnesenskii", 
                               "B. impatiens"),
                      values=my_palette_reduced) +
  theme(legend.position = c(0.15, 0.8),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size = 20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

q

