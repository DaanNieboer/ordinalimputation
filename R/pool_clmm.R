#' Calculate pooled estimates of fixed effects, standard deviation of the
#' random effect and conditional modes of the random effect terms.
#'
#' @param x an mira object returned from the with.mids function from the mice
#'          package.
#'
#' @return A clmm.pool object.
#'
pool.clmm <- function(x){
  # Pool fixed effects and standard deviation random effect.
  coefs        <- sapply(X = x$analyses, FUN = coefficients)
  std_re       <- sapply(X = lapply(X = x$analyses, FUN = get_re_std),
                         FUN = unlist)
  vcov_fits    <- lapply(X = x$analyses, FUN = vcov)
  coefs_fit    <- t(rbind(coefs, std_re))

  pool_fixed <- pool_rubin(coefs = coefs_fit, variance = vcov_fits)


  fixef  <- data.frame(variable = colnames(coefs_fit)[-ncol(coefs_fit)],
                       coefficient = pool_fixed$estimate[-ncol(coefs_fit)],
                       se = sqrt(diag(pool_fixed$variance)[-ncol(coefs_fit)]))

  std_re <- pool_fixed$estimate[ncol(coefs_fit)]

  ranef_fits   <- t(sapply(X = lapply(X = x$analyses, FUN = ranef), FUN = unlist))
  condVar_fits <- t(sapply(X = lapply(X = x$analyses, FUN = condVar), FUN = unlist))

  random_effect <- pool_re(ranef_fits, condVar_fits)

  res <- list(fixed_effects = fixef, stDev = std_re, random_effects = random_effect)
  class(res) <- c("pooled.clmm")
  return(res)
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
  qbar         <- colMeans(coefs_fit)

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
