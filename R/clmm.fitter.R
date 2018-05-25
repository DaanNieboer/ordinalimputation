#' A function to fit and pool results from a random effects with ordinal outcome
#' using multiple imputation. Uncertainty in model parameters is calculated
#' based on the model or through bootrstrapping.
#'
#' @param formula a two-sided linear formula object describing the fixed-effects
#'                part of the model, with the response on the left of a ~
#'                operator and the terms, separated by + operators, on the
#'                right. The vertical bar character "|" separates an expression
#'                for a model matrix and a grouping factor.
#' @param data A mids object containing the imputed data
#' @param bootstrap an optional parameter indicitating if a 95% CI should be
#'                  estimated using bootstrapping or model results.
#'                  Default is bootrstrapping.
#' @param n_boot optional parameter denoting the number of bootstrap samples
#'               that need to be drawn.
#' @return A pooled.clmm object containing the pooled parameter estimates.
#'
clmm.fitter <- function(formula, data, bootstrap = TRUE, n_boot = 1000,...){
  data <- complete(data, action = "long", include = FALSE)
  data <- split(x = data, f = data$.imp)

  if(boot){
    cluster <- attributes(terms(formula))$term.labels
    cluster <- cluster[grepl("\\|", cluster)]
    cluster <- trimws(gsub("\\|", "", gsub("1", "", cluster)))

    fit0      <- clmm_wrap(data = data[[1]], formula = formula)
    n_fixed   <- length(fit0$coefficients)
    n_random  <- length(fit0$ranef)

    res <- lapply(data, FUN = bootstrap.clmm, formula = formula,
                  cluster = cluster, n_fixed = n_fixed, n_random = n_random,
                  n_boot = n_boot,...)

    fixed_effect  <- do.call("rbind", sapply(res, "[", "fixed_effects"))
    random_effect <- do.call("rbind", sapply(res, "[", "random_effects"))
    std           <- unlist(sapply(res, "[", "std"))

    mu_fixed      <- colMeans(fixed_effect)
    mu_random     <- colMeans(random_effect)

    ci_l_fixed    <- apply(fixed_effect, 2, quantile, probs = 0.025)
    ci_h_fixed    <- apply(fixed_effect, 2, quantile, probs = 0.975)

    fixed           <- cbind(mu_fixed, ci_l_fixed, ci_h_fixed)
    rownames(fixed) <- names(fit0$coefficients)
    colnames(fixed) <- c("Estimate", "Lower 95% CI", "Upper 95% CI")

    ci_l_random     <- apply(random_effect, 2, quantile, probs = 0.025)
    ci_h_random     <- apply(random_effect, 2, quantile, probs = 0.975)

    random           <- cbind(mu_random, ci_l_random, ci_h_random)
    rownames(random) <- 1:nrow(random)
    colnames(random) <- c("Estimate", "Lower 95% CI", "Upper 95% CI")

    std <- c(mean(std), quantile(std, probs = c(0.025, 0.975)))
    names(std) <- c("Estimate", "Lower 95% CI", "Upper 95% CI")
    mor     <- exp(sqrt(2 * std^2) * qnorm(0.75))

    std <- rbind(std, mor)
    rownames(std) <- c("Standard deviation", "MOR")


    res <- list(fixed_effects = fixed, random_effects = random, random_dist = std)
    class(res) <- c("pooled.clmm")
    res$bootstrap <- TRUE
  }else{
    res          <- lapply(data, clmm_wrap, formula = formula,...)
    res_pooled   <- pooled(res)
    res_pooled$analyses  <- res
    res_pooled$bootstrap <- FALSE
    class(res) <- "pooled.clmm"
  }
  return(res)
}

bootstrap.clmm <- function(x, n_boot, formula, cluster, n_fixed, n_random,...){
  data_split <- split(x, x[, cluster])

  fixed_effects  <- matrix(nrow = n_boot, ncol = n_fixed)
  random_effects <- matrix(nrow = n_boot, ncol = n_random)
  std            <- rep(0, n_boot)
  for(i in 1:n_boot){
    dta_boot <- do.call("rbind", lapply(data_split, sample_pats))
    fit      <- clmm_wrap(data = dta_boot, formula,...)

    fixed_effects[i, ]  <- fit$coefficients
    random_effects[i, ] <- fit$ranef
    std[i]              <- fit$ST
  }
  res <- list(fixed_effects = fixed_effects, random_effects = random_effects, std = std)
  return(res)
}

sample_pats <- function(x){
  x <- x[sample(1:nrow(x), replace = TRUE), ]
  return(x)
}

clmm_wrap <- function(data,...){
  fit <- clmm(data = data,...)
  return(fit)
}
