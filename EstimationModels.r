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
obs <- tibble::tibble(
  t_obs = 1:20,
  c_obs = c(
    13.48, 23.66, 28.17, 29.34, 28.02, 25.56, 23.25, 
    21.52, 21.21, 19.30, 14.54, 14.62, 11.71, 12.21, 
    10.21, 10.82, 8.43, 8.79, 7.60, 6.03
  )
)

t_fit <- seq(0, 20, .2)
two_cpt_data <- list(
  n_obs = nrow(obs),
  c_obs = obs$c_obs,
  t_obs = obs$t_obs,
  a0 = c(1000, 0, 0),
  t0 = -.01,
  t_fit = t_fit,
  n_fit = length(t_fit)
)
