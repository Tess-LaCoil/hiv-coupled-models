//HIV model from Phillips 1996
functions {
  vector hiv_des(real time,
                 vector state,
                 real gamma,
                 real tau,
                 real nu,
                 real beta,
                 real rho,
                 real alpha,
                 real delta,
                 real phi,
                 real sigma) {

    // convenience...
    real R = state[1]; // Susceptible CD4+ cells
    real L = state[2]; // Latently infected CD4+ cells
    real E = state[3]; // Actively infected CD4+ cells
    real V = state[4]; // Amount of virus

    // derivative of state vector with respect to time
    vector[4] dadt;

    // compute derivatives
    dadt[1] = gamma*tau - nu*R - beta*V*R;
    dadt[2] = rho*beta*V*R - nu*L - alpha*L;
    dadt[3] = (1 - rho)*beta*V*R + alpha*L - delta*E;
    dadt[4] = phi*E - sigma*V;

    return dadt;
  }
}

data {
  int<lower=1> n_obs;
  int<lower=1> n_fit;
  vector[n_obs] R_obs; 
  vector[n_obs] L_obs; 
  vector[n_obs] E_obs; 
  vector[n_obs] V_obs; 
  array[n_obs] real t_obs; //Array of times
  array[n_fit] real t_fit; //Array of times
  vector[4] a0; //Initial conditions
  real t0;
}

parameters {
  real<lower=.0001> gamma;
  real<lower=.0001,upper=1> tau;
  real<lower=.0001> nu;
  real<lower=.0001> beta;
  real<lower=.0001,upper=1> rho;
  real<lower=.0001> alpha;
  real<lower=.0001> delta;
  real<lower=.0001> phi;
  real<lower=.0001> sigma;
  
  real<lower=.0001> sigma_error;
}

transformed parameters {
  // use ode solver to find all amounts at all event times
  array[n_obs] vector[4] amount = ode_rk45(hiv_des,
                                           a0,
                                           t0,
                                           t_obs,
                                           gamma,
                                           tau,
                                           nu,
                                           beta,
                                           rho,
                                           alpha,
                                           delta,
                                           phi,
                                           sigma);

  // vector of predicted counts
  vector[n_obs] mu;
  for (j in 1:n_obs) {
    mu[j] = amount[j, 2] ;
  }
}

model {
  // priors
  gamma ~ lognormal(log(1), 1); //Production rate of CD4+
  tau ~ beta(2, 2); //Fraction susceptible
  nu ~ lognormal(log(1), 1); //Removal rate of CD4+
  beta ~ lognormal(log(1), 1); //Rate of T-cell infection
  rho ~ beta(2, 2); //Fraction of infected cells that are latently infected
  alpha ~ lognormal(log(1), 1); //Rate of latent cells becoming activated
  delta ~ lognormal(log(1), 1); //Death/removal rate for actively infected cells
  phi ~ lognormal(log(1), 1); //Virus production rate
  sigma ~ lognormal(log(1), 1); //Virus removal rate
  
  sigma_error ~ lognormal(log(0.1), 1); //Measurement error

  // likelihood
  R_obs ~ normal(mu[1], sigma_error);
  L_obs ~ normal(mu[2], sigma_error);
  E_obs ~ normal(mu[3], sigma_error);
  V_obs ~ normal(mu[4], sigma_error);
}

generated quantities {
  array[n_fit] vector[4] hiv_fit = ode_rk45(hiv_des,
                                           a0,
                                           t0,
                                           t_fit,
                                           gamma,
                                           tau,
                                           nu,
                                           beta,
                                           rho,
                                           alpha,
                                           delta,
                                           phi,
                                           sigma);
  
  vector[n_fit] fitted_vals;
  for (j in 1:n_fit) {
    fitted_vals[j] = hiv_fit[j, 2];
  }
}
