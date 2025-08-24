if(FALSE){
  install.packages(
    "dust2",
    repos = c("https://mrc-ide.r-universe.dev", 
              "https://cloud.r-project.org"))
  
  install.packages(
    "monty",
    repos = c("https://mrc-ide.r-universe.dev", 
              "https://cloud.r-project.org"))
  
  install.packages(
    "odin2",
    repos = c("https://mrc-ide.r-universe.dev", 
              "https://cloud.r-project.org"))
}

library(dust2)
library(monty)
library(odin2)
library(deSolve)
library(tidyverse)

#Build data
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




