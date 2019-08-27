#' Calculate pooled estimates of fixed effects, standard deviation of the
#' random effect and conditional modes of the random effect terms.
#'
#' @param x A list containing clmm fists. For instance the analyses slot from
#'          a mira object of the mice package.
#' @param conf.int.re Type of 95% confidence interval for the variance of the
#'                    random effect to be calculated. Default is none.
#' @param data mira object used to fit the model. Needed for calculating the profile likelihood.
#'
#' @return A pooled.clmm object. With the following elements:
#' @slot fixed_effects Pooled fixed effect estimates and associated standard errors.
#' @slot random_effects Pooled modes of random effects with associated conditional variances.
#' @slot random_dist Standard deviation of the random effect distribution and associated median odds ratio
#' @slot conf_in_re Optional 95% confidence interval of the variance of the random effects.
pooling.clmm <- function(x, conf.int.re = c("none", "profile"), data = NULL){
  conf.int.re <- match.arg(conf.int.re)

  # Pool fixed effects and standard deviation random effect.
  coefs        <- sapply(X = x, FUN = coefficients)
  std_re       <- sapply(X = lapply(X = x, FUN = get_re_std),
                         FUN = unlist)
  coefs_fit    <- t(coefs)

  vcov_fits  <- lapply(X = x, FUN = get_vcov, coefs = rownames(coefs))
  pool_fixed <- pool_rubin(coefs = coefs_fit, variance = vcov_fits)

  if(class(std_re)=="numeric"){
    n_re <- 1
  }else{
    n_re <- nrow(std_re)
  }

  mu_fixed   <- pool_fixed$estimate
  se_fixed   <- sqrt(diag(pool_fixed$variance))
  ci_l_fixed <- mu_fixed - 1.96 * se_fixed
  ci_u_fixed <- mu_fixed + 1.96 * se_fixed

  fixef      <- cbind(mu_fixed, ci_l_fixed, ci_u_fixed)

  rownames(fixef) <- names(pool_fixed$estimate[1:(ncol(coefs_fit) - n_re)])
  colnames(fixef) <- c("Estimate", "Lower 95% CI", "Upper 95% CI")

  std_dev <- pool_fixed$estimate[(ncol(coefs_fit) - n_re + 1):ncol(coefs_fit)]
  mor     <- exp(sqrt(2 * std_dev^2) * qnorm(0.75))
  std_re  <- data.frame(std_dev = std_dev, mor = mor)
  if(conf.int.re=="profile"){
    if(is.null(data)){
      stop("To calculate the profile likelihood please supply the data used when fitting the model.")
    }
    cat("Calculating profile likelihood, may take a very long time.\n")
    conf_int <- matrix(nrow = n_re, ncol = 2)
    for(i in 1:n_re){
      conf_int[i, ] <- prof_ci(fits = x, index = i, data = data)
      cat("Random effect variabce number ", i, "done\n")
    }
  }else{
    conf_int <- NULL
  }

  ranef_fits   <- t(sapply(X = lapply(X = x, FUN = ranef), FUN = unlist))
  condVar_fits <- t(sapply(X = lapply(X = x, FUN = condVar), FUN = unlist))

  random_effect <- pool_re(ranef_fits, condVar_fits)

  mu_random <- random_effect$mode
  se_random <- sqrt(random_effect$cond_var)

  ci_l_random <- mu_random - 1.96 * se_random
  ci_u_random <- mu_random + 1.96 * se_random

  random <- cbind(mu_random, ci_l_random, ci_u_random)


  rownames(random) <- 1:nrow(random)
  colnames(random) <- c("Estimate", "Lower 95% CI", "Upper 95% CI")

  res <- list(fixed_effects = fixef, random_dist = std_re,
              random_effects = random, se_fixed = se_fixed,
              se_random = se_random, conf_int_re = conf_int)
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

get_vcov <- function(x, coefs){
  res <- vcov(x)
  res <- res[coefs, coefs]
  return(res)
}
