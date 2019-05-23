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

prof_ci <- function(fits, alpha = 0.05, index){
  # Function to calculate the profile confidence interval around the
  # variance of the random effect,

  xmax    <- 5 * unlist(fits[[1]]$ST)[index]
  ci.lb   <- optimize(f = score_prog_ll, interval = c(0, xmax),
                              alpha = alpha/2, fits = fits)
  ci.ub   <- optimize(f = score_prog_ll, interval = c(0, xmax),
                              alpha = 1 - alpha/2, fits = fits)

  return(c(ci.lb$par, ci.ub$par))
}

score_prog_ll <- function(x, fits, index, alpha){
  ll_value <- sapply(X = fits, FUN = prof_ll, x = x, index = index)
  avg_ll   <- mean(ll_value)

  score <- (avg_ll - alpha)^2

  return(score)
}
prof_ll <- function(fitted, x, index){
  start <- list(fitted$coefficients, unlist(fitted$ST))
  gamma <- start[[2]][index]
  start[[2]][index] <- x

  fit_prof <- update(fitted, start = start, eval.max = 1)
  res <- pnorm(sign(x - gamma) * sqrt(-2 * (fit_prof$logLik - fitted$logLik)))

  return(res)
}
