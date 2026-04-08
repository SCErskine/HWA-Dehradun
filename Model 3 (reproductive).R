data <- Dehradun_data_final
#data=read.csv(file.choose())
data$Ster=as.numeric(data$`% Total Sterilized`)
data$Pups=as.numeric(data$`% Pups`)
data$Lact=as.numeric(data$`% Lactating`)
data$iSeason <- factor(data$Season, levels = c("Non-breeding", "Breeding"))

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

# Pups model
Pup <- brm(
  formula = bf(
    scale(Pups) ~ 1 + scale(Ster) + (1 | Zone_id) + iSeason
  ),
  data = data,
  prior = c(
    prior(normal(0, 1), class = Intercept),
    prior(normal(0, 1), class = b,       coef = scaleSter),
    prior(normal(0, 1), class = b, coef = iSeasonBreeding),
    prior(exponential(1), class = sd),
    prior(exponential(1), class = sigma)
  ),
  chains = 4, cores = 4, iter = 8000,
  save_pars = save_pars(all = TRUE),
  seed = 123
)

summary(Pup)
posterior_summary(Pup)
plot(Pup)

# Lactation model
Lac <- brm(
  formula = bf(
    scale(Lact) ~ 1 + scale(Ster) + (1 | Zone_id) + iSeason
  ),
  data = data,
  prior = c(
    prior(normal(0, 1), class = Intercept),
    prior(normal(0, 1), class = b,       coef = scaleSter),
    prior(normal(0, 1), class = b, coef = iSeasonBreeding),
    prior(exponential(1), class = sd),
    prior(exponential(1), class = sigma)
  ),
  chains = 4, cores = 4, iter = 8000,
  save_pars = save_pars(all = TRUE),
  seed = 123
)

summary(Lac)
posterior_summary(Lac)
plot(Lac)


## Unstandardize main effects

# Get the SD of density
Pup_sd <- sd(data$Pups)
Lac_sd <- sd(data$Lact)

# Extract posterior summary for the Ster coefficient
Pup_coef <- fixef(Pup)["scaleSter", ]
Lac_coef <- fixef(Lac)["scaleSter", ]

# Get SD of sterilisation
Ster_sd <- sd(data$Ster)

# Effect in skin % per 10 percentage point (0.10) change in sterilisation
# Divide by Ster_sd to get per unit of Ster, multiply by 0.10 for 10pp change
(Pup_coef * Pup_sd) / Ster_sd * 0.10 * 100
(Lac_coef * Lac_sd) / Ster_sd * 0.10 * 100


## COUNTERFACTUAL PLOTS

library(ggplot2)

# Posterior draws
Pup_draws <- as_draws_df(Pup)
Lac_draws <- as_draws_df(Lac)

# Sterilisation sequence on original scale
Ster_seq <- seq(min(data$Ster), max(data$Ster), length.out = 100)
Ster_scaled <- (Ster_seq - mean(data$Ster)) / sd(data$Ster)

# Means and SDs for back-transformation
Pup_mean <- mean(data$Pups)
Lac_mean <- mean(data$Lact)

# Prediction matrices for each season
# Non-breeding (reference, iSeasonBreeding = 0)
Pup_pred_NB <- outer(Pup_draws$b_scaleSter, Ster_scaled) + Pup_draws$b_Intercept
Lac_pred_NB <- outer(Lac_draws$b_scaleSter, Ster_scaled) + Lac_draws$b_Intercept

# Breeding (iSeasonBreeding = 1)
Pup_pred_B <- outer(Pup_draws$b_scaleSter, Ster_scaled) + Pup_draws$b_Intercept + Pup_draws$b_iSeasonBreeding
Lac_pred_B <- outer(Lac_draws$b_scaleSter, Ster_scaled) + Lac_draws$b_Intercept + Lac_draws$b_iSeasonBreeding

# Build prediction dataframes
make_pred_df <- function(pred_NB, pred_B, Ster_seq, outcome_sd, outcome_mean, outcome_col) {
  rbind(
    data.frame(
      Ster   = Ster_seq,
      mean   = apply(pred_NB, 2, mean)             * outcome_sd + outcome_mean,
      lower  = apply(pred_NB, 2, quantile, 0.025)  * outcome_sd + outcome_mean,
      upper  = apply(pred_NB, 2, quantile, 0.975)  * outcome_sd + outcome_mean,
      Season = "Non-breeding"
    ),
    data.frame(
      Ster   = Ster_seq,
      mean   = apply(pred_B, 2, mean)             * outcome_sd + outcome_mean,
      lower  = apply(pred_B, 2, quantile, 0.025)  * outcome_sd + outcome_mean,
      upper  = apply(pred_B, 2, quantile, 0.975)  * outcome_sd + outcome_mean,
      Season = "Breeding"
    )
  )
}

Pup_pred_df <- make_pred_df(Pup_pred_NB, Pup_pred_B, Ster_seq, Pup_sd, Pup_mean)
Lac_pred_df <- make_pred_df(Lac_pred_NB, Lac_pred_B, Ster_seq, Lac_sd, Lac_mean)

# Raw data
Pup_raw <- data.frame(Ster = data$Ster, Outcome = data$Pups, Season = as.character(data$iSeason))
Lac_raw <- data.frame(Ster = data$Ster, Outcome = data$Lact, Season = as.character(data$iSeason))

# Plot function
make_plot <- function(pred_df, raw_df, y_label) {
  ggplot() +
    geom_ribbon(data = pred_df,
                aes(x = Ster, ymin = lower, ymax = upper, fill = Season),
                alpha = 0.2) +
    geom_line(data = pred_df,
              aes(x = Ster, y = mean, colour = Season),
              linewidth = 1) +
    geom_point(data = raw_df,
               aes(x = Ster, y = Outcome, colour = Season),
               alpha = 0.6, size = 2) +
    scale_colour_manual(values = c("Breeding" = "black", "Non-breeding" = "steelblue")) +
    scale_fill_manual(values   = c("Breeding" = "black", "Non-breeding" = "steelblue")) +
    labs(
      x      = "Sterilisation Proportion",
      y      = y_label,
      colour = "Season",
      fill   = "Season"
    ) +
    theme_classic() +
    theme(legend.position = c(0.85, 0.85))
}

# Render plots in consecutive windows
make_plot(Pup_pred_df, Pup_raw, "Proportion Pups")

make_plot(Lac_pred_df, Lac_raw, "Proportion Lactating Females")

