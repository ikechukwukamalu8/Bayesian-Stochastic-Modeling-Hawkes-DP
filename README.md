# Bayesian Hawkes Process and Nonparametric Modeling for Complex Systems

## 🎯 Project Objective
This project demonstrates the application of stochastic process modeling (Hawkes Process) and Nonparametric Bayesian methods (Dirichlet Process Mixture & Latent Intensity Estimation) to analyze self-exciting events in real-world time-series data. The goal is to make complex, temporal dependencies tractable and provide robust uncertainty quantification under rigid model constraints.

## 🛠️ Key Technical Skills Demonstrated
* **Stochastic Processes:** Implementation and parameter estimation (MLE) of the Hawkes Process.
* **Nonparametric Bayesian Methods:** Application of Dirichlet Process Mixture Models for cluster identification and Latent Intrinsic Autoregressive priors for intensity estimation.
* **Bayesian Inference & MCMC:** Utilizing the JAGS engine to compile graphical models and simulate full posterior distributions.
* **Model Diagnostics & Risk Assessment:** Executing Transformed Residual Analysis via the Time-Rescaling Theorem (QQ-Plots, Kolmogorov-Smirnov Test) to explicitly diagnose model limitations.

## 🔄 Methodological Workflow (Narrative)
1. **Parametric Modeling:** A standard Hawkes Process is fitted to discrete epidemiological event times using Maximum Likelihood Estimation (MLE).
2. **Diagnostic Validation:** The Time-Rescaling Theorem transforms the process residuals to evaluate model fit. Diagnostics show systemic flatlining and reject exponentiality ($p < 2.2 \times 10^{-16}$), proving a rigid parametric model fails to absorb localized temporal structural uncertainty.
3. **Nonparametric Resolution:** A Dirichlet Process Mixture Model isolates a rare, heavy-tailed anomaly regime (accounting for exactly 2.1% of the data). Finally, a Bayesian Nonparametric Intensity Sampler reconstructs the fluid velocity of the system alongside explicit 95% Bayesian Credible Intervals.

## 📄 Files
* `Hawkes_Process_Bayesian_Modeling.R`: Core script for data ingestion, Hawkes MLE optimization, and Dirichlet Process Mixture clustering.
* `Hawkes_Process_Diagnostic_Analysis.R`: Dedicated validation script running transformed residual analysis, QQ-plotting, and the Kolmogorov-Smirnov test.
* `Bayesian_Nonparametric_Intensity_Estimation.R`: Advanced research solution utilizing JAGS to fit a flexible random-walk smoothing prior over log-intensity for robust uncertainty quantification.
