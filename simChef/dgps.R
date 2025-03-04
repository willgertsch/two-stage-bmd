# data generating functions
# return information about the data to be generated in a later step
# models selected by using highest posterior weight in MA: ma_mod = ma_dichotomous_fit(dat$dose, dat$events, dat$n)

kociba1978_dgp = create_dgp(
  .dgp_fun = function(N1 = 117, N = 234) {
    
    dat = data.frame(
      dose = c(0, 1, 10, 100),
      events = c(9, 3, 18, 34),
      n = c(86, 50, 50, 48)
    )
    max_dose = max(dat$dose)
    
    model_type = 'qlinear' 
    # mod = single_dichotomous_fit(
    #   D = dat$dose,
    #   Y = dat$events,
    #   N = dat$n,
    #   model_type = 'qlinear',
    #   fit_type = 'mle',
    #   BMR = 0.1,
    #   alpha = 0.025
    # )
    theta = c(-2.0897337,  0.0129241)
    bmd = 8.152252
    return(
      list(
        theta = theta,
        N1 = N1,
        N = N,
        dat = dat,
        max_dose = max_dose,
        model_type = model_type,
        bmd = bmd
      )
    )
    
  },
  .name = 'Kociba et al. 1978'
)

kociba_loglogistic_dgp = create_dgp(
  .dgp_fun = function(N1 = 117, N = 234) {
    
    dat = data.frame(
      dose = c(0.000001, 1, 10, 100),
      events = c(9, 3, 18, 34),
      n = c(86, 50, 50, 48)
    )
    max_dose = max(dat$dose)
    
    model_type = 'log-logistic' 
    # mod = single_dichotomous_fit(
    #   D = dat$dose,
    #   Y = dat$events,
    #   N = dat$n,
    #   model_type = 'log-logistic',
    #   fit_type = 'mle',
    #   BMR = 0.1,
    #   alpha = 0.025
    # )
    theta = c(-2.3619758, -3.1608187,  0.8760667)
    bmd = 3.0038944
    return(
      list(
        theta = theta,
        N1 = N1,
        N = N,
        dat = dat,
        max_dose = max_dose,
        model_type = model_type,
        bmd = bmd
      )
    )
    
  },
  .name = 'Kociba log-logistic'
)

NTP1982a_dgp = create_dgp(
  .dgp_fun = function(N1 = 111, N = 223) {
    
    dat = data.frame(
      dose = c(0, 1.4, 7.1, 71),
      events = c(5, 1, 3, 12),
      n = c(75, 49, 50, 49)
    )
    max_dose = max(dat$dose)
    model_type = 'hill'
    # mod = single_dichotomous_fit(
    #   D = dat$dose,
    #   Y = dat$events,
    #   N = dat$n,
    #   model_type = 'hill',
    #   fit_type = 'mle',
    #   BMR = 0.1,
    #   alpha = 0.025
    # )
    theta = c(-2.978925,  -1.346135, -19.314638,   8.441910)
    bmd = 9.781407
    return(
      list(
        theta = theta,
        N1 = N1,
        N = N,
        dat = dat,
        max_dose = max_dose,
        model_type = model_type,
        bmd = bmd
      )
    )
  },
  .name = 'NTP 1982'
)


JBRC_dgp = create_dgp(
  .dgp_fun = function(N1=99, N=198) {
    dat = ToxicR::dichotomousDR %>%
      dplyr::filter(ID == 131)
    
    max_dose = max(dat$dose)
    
    model_type = 'hill'
    # mod = single_dichotomous_fit(
    #   D = dat$dose,
    #   Y = dat$obs,
    #   N = dat$n,
    #   model_type = 'hill',
    #   fit_type = 'mle',
    #   BMR = 0.1,
    #   alpha = 0.025
    # )
    theta = c(-2.442347,  3.884780, -6.146044,  2.550849)
    bmd = 4.744152
    return(
      list(
        theta = theta,
        N1 = N1,
        N = N,
        dat = dat,
        max_dose = max_dose,
        model_type = model_type,
        bmd = bmd
      )
    )
  },
  .name = 'JBRC'
)

# deguelin data
deguelin_dgp = create_dgp(
  .dgp_fun = function(N1=146, N=292) {
    dat = drc::deguelin
    max_dose = max(dat$dose)
    model_type = 'log-logistic'
    # mod = ToxicR::single_dichotomous_fit(
    #   D = dat$dose,
    #   Y = dat$r,
    #   N = dat$n,
    #   model_type = 'log-logistic',
    #   fit_type = 'mle',
    #   BMR = 0.1,
    #   alpha = 0.025
    # )
    theta = c(-0.7363599, -12.0933855,4.1368046)
    bmd = 10.937786
    return(
      list(
        theta = theta,
        N1 = N1,
        N = N,
        dat = dat,
        max_dose = max_dose,
        model_type = model_type,
        bmd = bmd
      )
    )
  },
  .name = 'Deguelin'
)