#Implementing HIV model in RStan and odin2

library(tidyverse)
library(deSolve)
library(rstan)
#library(odin2)
#library(dust2)
#library(monty)

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
            mu  = 0.00136, 
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
obs_data <- as_tibble(s_CD4_HIV_dynamics_solution) %>%
  mutate(R_obs = R + rnorm(1,0,2),
         L_obs = L*(1 + rnorm(1,0,0.01)), #1% error
         V_obs = V*(1 + rnorm(1,0,0.01)),
         E_obs = E*(1 + rnorm(1,0,0.01))
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
                chains=2,
                cores=2)
saveRDS(fit, "output/fits/HIV_Stan.rds")