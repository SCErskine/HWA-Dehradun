data <- Dehradun_data_final
#data=read.csv(file.choose())
data$DSter=as.numeric(data$Change_Ster)
data$DDensity=as.numeric(data$Change_Density)

# Remove NAs upfront
data <- data[complete.cases(data$DDensity, data$DSter), ]

# install.packages("brms")
library(brms)
library(dplyr)

# Check what as.integer gives you
table(as.integer(factor(data$Zone)))
# Explicitly recode to clean integers
data$Zone_id <- as.integer(factor(data$Zone))
# Verify the mapping
table(data$Zone, data$Zone_id)

# Split by season
data_s1 <- data[data$Season == "Breeding", ]
data_s2 <- data[data$Season == "Non-breeding", ]

# Breeding Season model
DeltaB <- brm(
  formula = bf(
    scale(DDensity) ~ 1 + scale(DSter) + (1 | Zone_id)
  ),
  data = data_s1,
  prior = c(
    prior(normal(0, 1), class = Intercept),
    prior(normal(0, 1), class = b,       coef = scaleDSter),
    prior(exponential(1), class = sd),
    prior(exponential(1), class = sigma)
  ),
  chains = 4, cores = 4, iter = 8000,
  save_pars = save_pars(all = TRUE),
  seed = 123
)

summary(DeltaB)
posterior_summary(DeltaB)
plot(DeltaB)

# Non-Breeding Season model
set.seed(123)
DeltaNB <- brm(
  formula = bf(
    scale(DDensity) ~ 1 + scale(DSter) + (1 | Zone_id)
  ),
  data = data_s2,
  prior = c(
    prior(normal(0, 1), class = Intercept),
    prior(normal(0, 1), class = b,       coef = scaleDSter),
    prior(exponential(1), class = sd),
    prior(exponential(1), class = sigma)
  ),
  chains = 4, cores = 4, iter = 8000,
  save_pars = save_pars(all = TRUE),
  seed = 123
)

summary(DeltaNB)
posterior_summary(DeltaNB)
plot(DeltaNB)


## Unstandardize main effects for first differences models
# Get the SD of density change for each season
B_DDensity_sd <- sd(data_s1$DDensity)
NB_DDensity_sd <- sd(data_s2$DDensity)

# Extract posterior summary for the DSter coefficient
B_Dcoef <- fixef(DeltaB)["scaleDSter", ]
NB_Dcoef <- fixef(DeltaNB)["scaleDSter", ]

# Get SD of sterilisation change for each season
B_DSter_sd <- sd(data_s1$DSter)
NB_DSter_sd <- sd(data_s2$DSter)

# Effect in dogs/km per 10 percentage point (0.10) change in sterilisation change
(B_Dcoef * B_DDensity_sd) / B_DSter_sd * 0.10
(NB_Dcoef * NB_DDensity_sd) / NB_DSter_sd * 0.10


## PLOT COUNTERFACTUALS
library(ggplot2)
library(brms)

# Generate prediction grids on ORIGINAL scale
s1_grid <- data.frame(
  DSter = seq(min(data_s1$DSter), max(data_s1$DSter), length.out = 100)
)
s2_grid <- data.frame(
  DSter = seq(min(data_s2$DSter), max(data_s2$DSter), length.out = 100)
)

# Get posterior predictions - brms will scale internally to match the formula
s1_preds <- fitted(DeltaB, newdata = s1_grid, re_formula = NA)
s2_preds <- fitted(DeltaNB, newdata = s2_grid, re_formula = NA)

# Back-transform predictions from scaled outcome to original DDensity scale
s1_DDensity_mean <- mean(data_s1$DDensity)
s1_DDensity_sd   <- sd(data_s1$DDensity)
s2_DDensity_mean <- mean(data_s2$DDensity)
s2_DDensity_sd   <- sd(data_s2$DDensity)

# Build prediction dataframes
s1_plot <- data.frame(
  DSter  = s1_grid$DSter,
  mean   = s1_preds[,"Estimate"] * s1_DDensity_sd + s1_DDensity_mean,
  lower  = s1_preds[,"Q2.5"]    * s1_DDensity_sd + s1_DDensity_mean,
  upper  = s1_preds[,"Q97.5"]   * s1_DDensity_sd + s1_DDensity_mean,
  Season = "Breeding"
)
s2_plot <- data.frame(
  DSter  = s2_grid$DSter,
  mean   = s2_preds[,"Estimate"] * s2_DDensity_sd + s2_DDensity_mean,
  lower  = s2_preds[,"Q2.5"]    * s2_DDensity_sd + s2_DDensity_mean,
  upper  = s2_preds[,"Q97.5"]   * s2_DDensity_sd + s2_DDensity_mean,
  Season = "Non-breeding"
)

pred_df <- rbind(s1_plot, s2_plot)

# Raw data points
raw_df <- rbind(
  data.frame(DSter = data_s1$DSter, DDensity = data_s1$DDensity, Season = "Breeding"),
  data.frame(DSter = data_s2$DSter, DDensity = data_s2$DDensity, Season = "Non-breeding")
)

# Plot
ggplot() +
  geom_ribbon(data = pred_df,
              aes(x = DSter, ymin = lower, ymax = upper, fill = Season),
              alpha = 0.2) +
  geom_line(data = pred_df,
            aes(x = DSter, y = mean, colour = Season),
            linewidth = 1) +
  geom_point(data = raw_df,
             aes(x = DSter, y = DDensity, colour = Season),
             alpha = 0.6, size = 2) +
  scale_colour_manual(values = c("Breeding" = "black", "Non-breeding" = "steelblue")) +
  scale_fill_manual(values   = c("Breeding" = "black", "Non-breeding" = "steelblue")) +
  labs(
    x      = "Change in Sterilisation Proportion",
    y      = "Change in Dog Density (dogs/km)",
    colour = "Season",
    fill   = "Season"
  ) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.85))

