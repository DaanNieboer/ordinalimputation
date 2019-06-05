#' Compute the profiled likelihood for the standard deviation for the random
#' effect in a fitted cumulative mixed model.
#'
profile.clmm <- function(fitted, alpha = 0.01, n_steps = 200){
  n_re <- length(fitted$ST)
  wald_se <- sqrt(diag(vcov(fitted)))
  wald_se <- wald_se[(length(wald_se) - n_re + 1):length(wald_se)]

  mu <- log(unlist(fitted$ST))

  x_profile <- matrix(nrow = n_steps, ncol = n_re)
  for(i in 1:n_re){
    x_profile[, i] <- seq(mu[i] - qnorm(1 - alpha/2) * wald_se[i],
                          mu[i] + qnorm(1 - alpha/2) * wald_se[i],
                          length.out = n_steps)
  }

  y_profile <- matrix(nrow = n_steps, ncol = n_re)
  for(i in 1:n_re){
    mu_start    <- mu
    for(j in 1:n_steps){
      mu_start[i]   <- x_profile[j, i]
      start         <- list(fitted$coefficients, exp(mu_start))
      fit_prof      <- update(fitted, start = start, eval.max = 1)
      y_profile[j, i] <- -2 * (fit_prof$logLik - fitted$logLik)

    }
  }
  res <- list(x_profile, y_profile)
}

prof_ci <- function(fits, alpha = 0.05, index, data = data){
  # Function to calculate the profile confidence interval around the
  # variance of the random effect,
  x_est   <- unlist(fits[[1]]$ST)[index]
  xmax    <- 5 * x_est
  ci.lb   <- optimize(f = score_prog_ll, interval = c(0, x_est),
                              alpha = alpha/2, fits = fits, data = data)
  ci.ub   <- optimize(f = score_prog_ll, interval = c(x_est, xmax),
                              alpha = 1 - alpha/2, fits = fits, data = data)

  return(c(ci.lb$minimum, ci.ub$minimum))
}


score_prog_ll <- function (x, fits, index, alpha, data)
{
  n_x      <- length(x)
  score    <- rep(NA, n_x)

  for(j in 1:n_x){
    ll_value <- rep(NA , length(fits))
    for(i in 1:length(fits)){
      fitted <- fits[[i]]

      start             <- list(fitted$coefficients, unlist(fitted$ST))
      gamma             <- start[[2]][index]
      start[[2]][index] <- x[j]

      fit_prof <- update(fitted, start = start, eval.max = 1,
                         data = complete(data, i))
      ll_value[i] <- pnorm(sign(x[j] - gamma[index]) * sqrt(-2 * (fit_prof$logLik -
                                                                    fitted$logLik)))
    }

    avg_ll   <- mean(ll_value)
    score[j] <- (avg_ll - alpha)^2
  }
  return(score)
}
