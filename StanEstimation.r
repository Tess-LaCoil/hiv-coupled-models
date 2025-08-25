#Implementing HIV model in RStan and odin2

library(tidyverse)
library(deSolve)
library(rstan)
library(GGally)

#See up system with some data
# Timestep
t <- seq(0,120, by = 1)

# Variables
cells <- c(R = 200, 
           L = 0,
           V = 4e-07, 
           E = 0)

# Model parameters
params <- c(gamma = 1.36, 
            nu  = 0.00136, 
            tau  = 0.2, 
            beta = 0.00027,
            rho  = 0.1, 
            alpha = 0.036, 
            sigma = 2,
            delta  = 0.33,
            phi = 100)

# Equations
source("R/phillips_1996_viral_load.R")

# Solve the equations
s_CD4_HIV_dynamics_solution <-
  ode(
    y = cells,
    times = t,
    func = viral_load,
    parms = params,
    method = "rk4")

# load the plotting functions
source("R/phillips_1996_viz.R")

#Stan model adapted from https://blog.djnavarro.net/posts/2023-05-16_stan-ode/#multi-compartment-models
#Add measurement error
obs_data <- as_tibble(s_CD4_HIV_dynamics_solution)[seq(from = 1, 
                                                       to = nrow(s_CD4_HIV_dynamics_solution),
                                                       by = 5),] %>%
  mutate(R_obs = R + rnorm(1,0,0.001),
         L_obs = L*(1 + rnorm(1,0,0.001)), #1% error
         V_obs = V*(1 + rnorm(1,0,0.001)),
         E_obs = E*(1 + rnorm(1,0,0.001))
         )

stan_data <- list(
  n_obs = nrow(obs_data),
  n_fit = nrow(obs_data),
  R_obs = as.numeric(obs_data$R_obs), 
  L_obs = as.numeric(obs_data$L_obs), 
  E_obs = as.numeric(obs_data$E_obs), 
  V_obs = as.numeric(obs_data$V_obs), 
  t_obs = as.numeric(obs_data$time),
  t_fit = as.numeric(obs_data$time),
  a0 = as.numeric(obs_data[1,6:9]),
  t0 = -0.01
)

#Load model and fit
model <- stan_model(file="stan/HIV.stan")
fit <- sampling(model, data=stan_data, 
                iter=2000,
                chains=4,
                cores=4)
saveRDS(fit, "output/fits/HIV_Stan.rds")
rm(fit)

#Extract estimates
HIV_Stan <- readRDS("output/fits/HIV_Stan.rds")
samples <- rstan::extract(HIV_Stan, permuted=TRUE)
par_names <- c(
  "gamma", "tau", "nu", "beta", "rho", "alpha", 
  "delta", "phi", "sigma"
)

#Pairs plots and diagnostics
traceplot(HIV_Stan, pars=par_names, inc_warmup=TRUE)

pairs_plot_data <- bind_rows(samples[1:10]) %>%
  mutate(log_gamma = log(gamma),
         log_nu = log(nu),
         log_beta = log(beta),
         log_alpha = log(alpha),
         log_delta = log(delta),
         log_phi = log(phi),
         log_sigma = log(sigma)
         )
ggpairs(pairs_plot_data[, c(2,5, 10:17)])

true_vals <- tibble(
  true_val = params,
  par_name = names(params)
)

par_est_tibble <- tibble()
for(i in 1:length(par_names)){
  temp <- tibble(
    par_name = par_names[i],
    par_mean = mean(samples[[par_names[i]]]),
    par_median = median(samples[[par_names[i]]]),
    par_CI_lower = quantile(samples[[par_names[i]]],0.025),
    par_CI_upper = quantile(samples[[par_names[i]]],0.975),
  )
  par_est_tibble <- rbind(par_est_tibble, temp)
}

par_est_tibble <- left_join(par_est_tibble, true_vals, by = "par_name")
par_est_tibble <- par_est_tibble %>%
  mutate(true_in_ci = (true_val >= par_CI_lower) & (true_val <= par_CI_upper))
par_est_tibble
