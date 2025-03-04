# utility functions used in simulation study

# generate a dataset based on an approximate design
# N: total sample size
# theta: model parameters
# x: dose levels
# w: sample size allocations
sim_data = function(N, theta, x, w, model_type) {
  
  # construct data
  n = floor(N*w)
  
  # dose response function
  # convert to BMDS parameterization
  # copying from
  # see https://github.com/NIEHS/ToxicR/blob/main/R/dicho_functions.R

  if (model_type == 'qlinear') {
    g = 1 / (1 + exp(-theta[1]))
    a = theta[2]
    p = g + (1 - g) * 1 - exp(-a * x)
  }
  else if (model_type == 'hill') {
    g = 1 / (1 + exp(-theta[1]))
    nu = 1 / (1 + exp(-theta[2]))
    a = theta[3]
    b = theta[4]
    p = g + (1 - g) * nu * (1 / (1 + exp(-a - b * log(x))))
  }
  else if (model_type == 'log-logistic') {
    g = 1/(1 + exp(-theta[1]))
    a = theta[2]
    b = theta[3]
    p = g + (1-g)/(1+exp(-a-b*log(x)))
  }
  else
    stop('model type not supported')
  
  # sample
  events = rbinom(length(n), n, p)
  
  data.frame(
    dose = x,
    n = n,
    events = events,
    p = p,
    phat = events/n
  )
  
}

# compute length of BMD interval from simulated data
# for the ToxicR package, this is a profile likelihood interval
compute_BMD_CI_length = function(bmd) {
  
  #bmd_CI_length = as.numeric(mod$bmd[3] - mod$bmd[2])
  bmd_CI_length = as.numeric(bmd[3] - bmd[2])
  return(bmd_CI_length)
}

# compute difference between the true and estimated bmd
# will later be used to compute the bias
compute_bmd_diff = function(mod, bmd) {
  #as.numeric(mod$bmd[1] - bmd)
  as.numeric(mod[1] - bmd)
}

bmd_coverage = function(mod, bmd) {
  #bmdl = as.numeric(mod$bmd[2])
  #bmdu = as.numeric(mod$bmd[3])
  bmdl = as.numeric(mod[2])
  bmdu = as.numeric(mod[3])
  if (is.na(bmdl) | is.na(bmdu)) {
    return(0)
  }
  else if (bmdl <= bmd & bmdu >= bmd)
    return(1)
  else
    return(0)
}

# extract BMDL values
compute_BMDL = function(mod) {
  #bmdl = as.numeric(mod$bmd[2])
  bmdl = as.numeric(mod[2])
  return(bmdl)
}

# compute mean squared error
compute_theta_mse = function(mod, theta) {
  
  theta2 = theta
  #theta2[1] = 1 / (1 + exp(-theta2[1])) # rescale for better mse
  theta3 = mod$parameters
  #theta3[1] = 1 / (1 + exp(-theta3[1]))
  sum((theta3 - theta2)^2) / length(theta)
}

# compute bias of theta
compute_theta_bias = function(mod, theta) {
  theta_hat = mod$parameters
  return(theta_hat - theta)
}
