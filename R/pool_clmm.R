#' Calculate pooled estimates of fixed effects, standard deviation of the
#' random effect and conditional modes of the random effect terms.
#'
#' @param x an mira object returned from the with.mids function from the mice
#'          package.
#'
#' @return A pooled.clmm object. With the following elements:
#' @slot fixed_effects Pooled fixed effect estimates and associated standard errors.
#' @slot random_effects Pooled modes of random effects with associated conditional variances.
#' @slot random_dist Standard deviation of the random effect distribution and associated median odds ratio
#'
pooling.clmmmira <- function(x){
  # Pool fixed effects and standard deviation random effect.
  coefs        <- sapply(X = x, FUN = coefficients)
  std_re       <- sapply(X = lapply(X = x, FUN = get_re_std),
                         FUN = unlist)
  vcov_fits    <- lapply(X = x, FUN = vcov)
  coefs_fit    <- t(rbind(coefs, std_re))

  pool_fixed <- pool_rubin(coefs = coefs_fit, variance = vcov_fits)

  mu_fixed   <- pool_fixed$estimate[-ncol(coefs_fit)]
  se_fixed   <- sqrt(diag(pool_fixed$variance)[-ncol(coefs_fit)])
  ci_l_fixed <- mu_fixed - 1.96 * se_fixed
  ci_u_fixed <- mu_fixed + 1.96 * se_fixed
  fixef      <- cbind(mu_fixed, ci_l_fixed, ci_u_fixed)
  rownames(fixef) <- names(pool_fixed$estimate[-ncol(coefs_fit)])
  colnames(fixef) <- c("Estimate", "Lower 95% CI", "Upper 95% CI")

  std_dev <- pool_fixed$estimate[ncol(coefs_fit)]
  mor     <- exp(sqrt(2 * std_dev^2) * qnorm(0.75))
  std_re  <- data.frame(std_dev = std_dev, mor = mor)

  ranef_fits   <- t(sapply(X = lapply(X = x, FUN = ranef), FUN = unlist))
  condVar_fits <- t(sapply(X = lapply(X = x, FUN = condVar), FUN = unlist))

  random_effect <- pool_re(ranef_fits, condVar_fits)

  res <- list(fixed_effects = fixef, random_dist = std_re, random_effects = random_effect)
  class(res) <- c("pooled.clmm")
  return(res)
}
pooling <- function(x,...){
  UseMethod("pooling",x)
}


get_re_std <- function(x){
  res <- x$ST
  return(res)
}

pool_rubin <- function(coefs, variance){
  qbar <- colMeans(coefs)

  m            <- length(variance)

  bw_imp_var   <- var(coefs)
  wi_imp_var   <- Reduce("+", variance)/m

  totalVar     <- wi_imp_var + (1 + 1/m) * bw_imp_var
  qbar         <- colMeans(coefs)

  res <- list(estimate = qbar, variance = totalVar)
}

pool_re <- function(re_mode, condvar){
  n_cluster <- ncol(re_mode)
  m         <- nrow(re_mode)

  res <- data.frame(cluster = 1:n_cluster,
                    mode = rep(NA, n_cluster),
                    cond_var = rep(NA, n_cluster))
  for(i in 1:n_cluster){


    res$mode[i] <- mean(re_mode[, i])

    bw_imp_var <- var(re_mode[, i])
    wi_imp_var <- mean(condvar[, i])

    res$cond_var[i] <- wi_imp_var + (1 + 1/m) * bw_imp_var
  }
  return(res)
}

print.pooled.clmm <- function(x){
  cat("Random effect parameters:\n")
  print(x$random_dist)
  cat("\n\nFixed effect estimates:\n")
  print(x$fixed_effects)

}
