# naive designs
# for a given max dose and number of doses, give the naive design

unif_design = function(max_dose, num_doses) {
  x = seq(0, max_dose, length.out = num_doses)
  w = rep(1/num_doses, num_doses)
  return(list(x = x, w = w))
}


log_unif_design = function(max_dose, num_dose) {
  x = c(0, exp(seq(1, log(max_dose), by = (log(max_dose)-1)/(num_doses - 2))))
  w = rep(1/num_doses, num_doses)
  return(list(x = x, w = w))
}