data <- Dehradun_data_final
#data=read.csv(file.choose())
data$Ster=as.numeric(data$`% Total Sterilized`)
data$Skin=as.numeric(data$`% Skin Problem`)

# install.packages("brms")
library(brms)
library(dplyr)
# install.packages("tidybayes")

# Check what as.integer gives you
table(as.integer(factor(data$Zone)))
# Explicitly recode to clean integers
data$Zone_id <- as.integer(factor(data$Zone))
# Verify the mapping
table(data$Zone, data$Zone_id)

# Split by season
data_s1 <- data[data$Season == "Breeding", ]
data_s2 <- data[data$Season == "Non-breeding", ]

# Create unique time index within each season
data_s1 <- data_s1[order(data_s1$Zone_id, data_s1$Survey), ]
data_s2 <- data_s2[order(data_s2$Zone_id, data_s2$Survey), ]

data_s1$Time_id <- ave(seq_len(nrow(data_s1)), data_s1$Zone_id, FUN = seq_along)
data_s2$Time_id <- ave(seq_len(nrow(data_s2)), data_s2$Zone_id, FUN = seq_along)

# Breeding Season model
SkinB <- brm(
  formula = bf(
    scale(Skin) ~ 1 + scale(Ster) + (1 | Zone_id) +
      ar(time = Time_id, gr = Zone_id, p = 1)
  ),
  data = data_s1,
  prior = c(
    prior(normal(0, 1), class = Intercept),
    prior(normal(0, 1), class = b, coef = scaleSter),
    prior(normal(0, 1), class = ar),
    prior(exponential(1), class = sd),
    prior(exponential(1), class = sigma)
  ),
  chains = 4, cores = 4, iter = 8000,
  save_pars = save_pars(all = TRUE),
  seed = 123
)

summary(SkinB)
posterior_summary(SkinB)
plot(SkinB)

# Non-Breeding Season model
SkinNB <- brm(
  formula = bf(
    scale(Skin) ~ 1 + scale(Ster) + (1 | Zone_id) +
      ar(time = Time_id, gr = Zone_id, p = 1)
  ),
  data = data_s2,
  prior = c(
    prior(normal(0, 1), class = Intercept),
    prior(normal(0, 1), class = b, coef = scaleSter),
    prior(normal(0, 1), class = ar),
    prior(exponential(1), class = sd),
    prior(exponential(1), class = sigma)
  ),
  chains = 4, cores = 4, iter = 8000,
  save_pars = save_pars(all = TRUE),
  seed = 123
)

summary(SkinNB)
posterior_summary(SkinNB)
plot(SkinNB)


## Unstandardize main effects

# Get the SD of density for each season
B_Skin_sd <- sd(data_s1$Skin)
NB_Skin_sd <- sd(data_s2$Skin)

# Extract posterior summary for the Ster coefficient
B_coef <- fixef(SkinB)["scaleSter", ]
NB_coef <- fixef(SkinNB)["scaleSter", ]

# Get SD of sterilisation for each season
B_Ster_sd <- sd(data_s1$Ster)
NB_Ster_sd <- sd(data_s2$Ster)

# Effect in skin % per 10 percentage point (0.10) change in sterilisation
# Divide by Ster_sd to get per unit of Ster, multiply by 0.10 for 10pp change
(B_coef * B_Skin_sd) / B_Ster_sd * 0.10 * 100
(NB_coef * NB_Skin_sd) / NB_Ster_sd * 0.10 * 100



## Plot Counterfactuals

# Means and SDs for back-transformation
B_Skin_mean  <- mean(data_s1$Skin)
B_Skin_sd    <- sd(data_s1$Skin)
NB_Skin_mean <- mean(data_s2$Skin)
NB_Skin_sd   <- sd(data_s2$Skin)

B_Ster_mean  <- mean(data_s1$Ster)
B_Ster_sd    <- sd(data_s1$Ster)
NB_Ster_mean <- mean(data_s2$Ster)
NB_Ster_sd   <- sd(data_s2$Ster)

# Extract fixed effects
B_fe  <- fixef(SkinB)
NB_fe <- fixef(SkinNB)

# Sterilisation sequence on original scale
B_Ster_seq  <- seq(min(data_s1$Ster), max(data_s1$Ster), length.out = 100)
NB_Ster_seq <- seq(min(data_s2$Ster), max(data_s2$Ster), length.out = 100)

# Scale the predictor the same way brms did internally
B_Ster_scaled  <- (B_Ster_seq  - mean(data_s1$Ster)) / sd(data_s1$Ster)
NB_Ster_scaled <- (NB_Ster_seq - mean(data_s2$Ster)) / sd(data_s2$Ster)

library(tidybayes)

# Draw from posterior
B_draws  <- as_draws_df(SkinB)
NB_draws <- as_draws_df(SkinNB)

# For each sterilisation value, compute predicted density across all draws
B_Ster_seq  <- seq(min(data_s1$Ster), max(data_s1$Ster), length.out = 100)
NB_Ster_seq <- seq(min(data_s2$Ster), max(data_s2$Ster), length.out = 100)

B_Ster_scaled  <- (B_Ster_seq  - mean(data_s1$Ster)) / sd(data_s1$Ster)
NB_Ster_scaled <- (NB_Ster_seq - mean(data_s2$Ster)) / sd(data_s2$Ster)

# Compute predictions across draws then summarise
B_pred_matrix <- outer(B_draws$b_scaleSter, B_Ster_scaled) + B_draws$b_Intercept
NB_pred_matrix <- outer(NB_draws$b_scaleSter, NB_Ster_scaled) + NB_draws$b_Intercept

# Back-transform and summarise
B_plot <- data.frame(
  Ster  = B_Ster_seq,
  mean  = (apply(B_pred_matrix,  2, mean))  * B_Skin_sd  + B_Skin_mean,
  lower = (apply(B_pred_matrix,  2, quantile, 0.025)) * B_Skin_sd  + B_Skin_mean,
  upper = (apply(B_pred_matrix,  2, quantile, 0.975)) * B_Skin_sd  + B_Skin_mean,
  Season = "Breeding"
)

NB_plot <- data.frame(
  Ster  = NB_Ster_seq,
  mean  = (apply(NB_pred_matrix, 2, mean))  * NB_Skin_sd + NB_Skin_mean,
  lower = (apply(NB_pred_matrix, 2, quantile, 0.025)) * NB_Skin_sd + NB_Skin_mean,
  upper = (apply(NB_pred_matrix, 2, quantile, 0.975)) * NB_Skin_sd + NB_Skin_mean,
  Season = "Non-breeding"
)

pred_df <- rbind(B_plot, NB_plot)

# Raw data points
raw_df <- rbind(
  data.frame(Ster = data_s1$Ster, Skin = data_s1$Skin, Season = "Breeding"),
  data.frame(Ster = data_s2$Ster, Skin = data_s2$Skin, Season = "Non-breeding")
)

# Plot
ggplot() +
  geom_ribbon(data = pred_df,
              aes(x = Ster, ymin = lower, ymax = upper, fill = Season),
              alpha = 0.2) +
  geom_line(data = pred_df,
            aes(x = Ster, y = mean, colour = Season),
            linewidth = 1) +
  geom_point(data = raw_df,
             aes(x = Ster, y = Skin, colour = Season),
             alpha = 0.6, size = 2) +
  scale_colour_manual(values = c("Breeding" = "black", "Non-breeding" = "steelblue")) +
  scale_fill_manual(values   = c("Breeding" = "black", "Non-breeding" = "steelblue")) +
  labs(
    x      = "Sterilisation Proportion",
    y      = "Proportion of Dogs with Skin Conditions",
    colour = "Season",
    fill   = "Season"
  ) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.85))
