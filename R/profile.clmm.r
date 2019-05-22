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
      y_profile[j, i] <- sign(mu[i] - x_profile[j, i]) * sqrt(-2 * (fit_prof$logLik - fitted$logLik))

    }
  }
  x_min <- x_profile[which.min((pnorm(y_profile) - 1 - 0.05/2)^2)]
  x_max <- x_profile[which.min((pnorm(y_profile) - 0.05/2)^2)]

  limits <- exp(c(x_min, x_max))
}

ll_prof <- function(mu_prof, start, j, fitted){

}
