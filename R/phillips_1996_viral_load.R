viral_load <- function(t, cells, params) {
  with(as.list(c(cells, params)), {
    dR_dt = (gamma * tau) - (nu * R) - (beta * V * R)
    
    dL_dt =  (rho * beta * V * R) - (nu * L) - (alpha * L)
    
    dV_dt = (phi * E) - (sigma * V)
    
    dE_dt = ((1 - rho) * beta * V * R) + (alpha * L) - (delta * E)
    
    return(list(c(dR_dt, dL_dt, dV_dt, dE_dt)))
    
  })
}
