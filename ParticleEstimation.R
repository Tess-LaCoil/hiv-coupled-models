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

#Model adapted from https://blog.djnavarro.net/posts/2023-05-16_stan-ode/#multi-compartment-models
#Add measurement error
obs_data <- as_tibble(s_CD4_HIV_dynamics_solution)[seq(from = 1, 
                                                       to = nrow(s_CD4_HIV_dynamics_solution),
                                                       by = 5),] %>%
  mutate(R_obs = R + rnorm(1,0,0.001),
         L_obs = L*(1 + rnorm(1,0,0.001)), #0.1% error
         V_obs = V*(1 + rnorm(1,0,0.001)),
         E_obs = E*(1 + rnorm(1,0,0.001))
  )

# Model implementation
HIV <- odin({
  initial(R) <- 0
  initial(L) <- 0
  initial(E) <- 0
  initial(V) <- 0
  initial(estimated, zero_every = 1) <- 0
  
  update(R) <- gamma*tau - nu*R - beta*V*R
  update(L) <- rho*beta*V*R - nu*L - alpha*L
  update(E) <- (1 - rho)*beta*V*R + alpha*L - delta*E
  update(V) <- phi*E - sigma*V
  update(incidence) <- incidence + n_SI
  
  n_SI <- Binomial(S, p_SI)
  n_IR <- Binomial(I, p_IR)
  p_SI <- 1 - exp(-beta * I / N * dt)
  p_IR <- 1 - exp(-gamma * dt)
  
  gamma <- parameter(differentiate = TRUE)
  tau <- parameter(differentiate = TRUE)
  nu <- parameter(differentiate = TRUE)
  beta <- parameter(differentiate = TRUE)
  rho <- parameter(differentiate = TRUE)
  alpha <- parameter(differentiate = TRUE)
  delta <- parameter(differentiate = TRUE)
  phi <- parameter(differentiate = TRUE)
  sigma <- parameter(differentiate = TRUE)
  
  observations <- data()
  observations ~ Normal(estimated)
}, quiet = TRUE)


