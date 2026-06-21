# ==============================================================================
# ADVANCED RESEARCH SOLUTION: Nonparametric Bayesian Intensity Estimation
# FOCUS: Well-Calibrated Uncertainty Quantification for Complex Point Processes
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(rjags)
library(coda)

cat("🚀 Initializing Nonparametric Bayesian Intensity Sampler...\n")

# 1. Generate a regular evaluation grid across our timeline
# (Using the event_times data structure calculated in Script 1)
grid_size <- 50
time_grid <- seq(0, max(event_times), length.out = grid_size)

# Compute empirical event counts falling into each grid interval
bin_counts <- hist(event_times, breaks = time_grid, plot = FALSE)$counts
N_bins <- length(bin_counts)
mid_points <- time_grid[-1]

# 2. Define a Nonparametric Bayesian Model in BUGS/JAGS syntax
# Instead of forcing an exponential curve, we place a flexible Random Walk (Prior)
# over the log-intensity to allow the data to dictate the shape dynamically.
jags_model_string <- "
model {
  # Likelihood: Counts in each interval follow a Poisson process driven by intensity lambda
  for (t in 1:N_bins) {
    bin_counts[t] ~ dpois(lambda[t])
    log(lambda[t]) <- intensity_lat[t]
  }
  
  # Nonparametric Prior: First-order Random Walk (Intrinsic Autoregressive Prior)
  # This acts as a smoothing prior, allowing flexible non-parametric shapes.
  intensity_lat[1] ~ dnorm(0, 0.01)
  for (t in 2:N_bins) {
    intensity_lat[t] ~ dnorm(intensity_lat[t-1], tau)
  }
  
  # Hyperprior on the smoothing precision (controls model complexity)
  tau ~ dgamma(1, 0.05)
  sigma <- 1 / sqrt(tau)
}
"

# 3. Package data for the JAGS engine
jags_data <- list(
  bin_counts = bin_counts,
  N_bins = N_bins
)

# 4. Compile and run the Markov Chain Monte Carlo (MCMC) sampler
cat("⏳ Compiling JAGS model and running MCMC chains...\n")
model_compile <- jags.model(textConnection(jags_model_string), 
                             data = jags_data, 
                             n.chains = 2, 
                             n.adapt = 500)

update(model_compile, 500) # Burn-in period

posterior_samples <- coda.samples(model_compile, 
                                  variable.names = c("lambda", "sigma"), 
                                  n.iter = 2000)

# 5. Extract posterior distributions for well-calibrated Uncertainty Quantification
samples_matrix <- as.matrix(posterior_samples)
lambda_cols <- grep("lambda", colnames(samples_matrix))

# Calculate the posterior mean and the exact 95% Bayesian Credible Intervals
post_mean <- colMeans(samples_matrix[, lambda_cols])
post_lower <- apply(samples_matrix[, lambda_cols], 2, quantile, probs = 0.025)
post_upper <- apply(samples_matrix[, lambda_cols], 2, quantile, probs = 0.975)

# Create a clean reporting dataframe
plot_df <- data.frame(
  Time = mid_points,
  Observed = bin_counts,
  PosteriorMean = post_mean,
  LowerCI = post_lower,
  UpperCI = post_upper
)

# 6. Visualize the Well-Calibrated Uncertainty Bounds
cat("\n📊 Plotting Nonparametric Bayesian Intensity and Uncertainty Bounds...\n")
ggplot(plot_df, aes(x = Time)) +
  geom_point(aes(y = Observed, color = "Observed Event Counts"), alpha = 0.6, size = 2) +
  geom_line(aes(y = PosteriorMean, color = "Posterior Mean Intensity"), size = 1.2) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), fill = "purple", alpha = 0.15) +
  scale_color_manual(values = c("Observed Event Counts" = "black", "Posterior Mean Intensity" = "purple")) +
  theme_minimal() +
  labs(
    title = "Nonparametric Bayesian Intensity Estimation",
    subtitle = "Shaded area tracks the explicit 95% Bayesian Credible Interval (Uncertainty Quantification)",
    x = "Timeline (Days)",
    y = "Event Intensity / Counts",
    color = "Legend"
  ) +
  theme(legend.position = "bottom")
