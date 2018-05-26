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
  n    <- nrow(complete(data, 1))
  data <- complete(data, action = "long", include = FALSE)
  data <- split(x = data, f = data$.imp)

  res                  <- lapply(data, clmm_wrap, formula = formula,...)
  res_pooled           <- pooling.clmm(res)
  res_pooled$analyses  <- res
  if(bootstrap){
    cluster <- attributes(terms(formula))$term.labels
    cluster <- cluster[grepl("\\|", cluster)]
    cluster <- trimws(gsub("\\|", "", gsub("1", "", cluster)))

    fit0      <- clmm_wrap(data = data[[1]], formula = formula)
    n_fixed   <- length(fit0$coefficients)
    n_random  <- length(fit0$ranef)

    fixed_effects  <- matrix(nrow = n_boot, ncol = n_fixed)
    random_effects <- matrix(nrow = n_boot, ncol = n_random)
    std_dev        <- matrix(nrow = n_boot, ncol = 2)

    for(i in 1:n_boot){
      index    <- sample(1:n, replace = TRUE, size = n)
      dta_boot <- lapply(data, "[", index, )
      res_boot <- try(lapply(dta_boot, clmm_wrap, formula = formula,...))
      res_boot <- pooling.clmm(res_boot)
      if(class(res_boot)[1]!="try-error"&res_boot$fixed_effects==n_fixed){
        fixed_effects[i, ]  <- res_boot$fixed_effects
        random_effects[i, ] <- res_boot$random_effects
        std_dev[i, ]        <- res_boot$random_dist
      }else{
        warning("Unable to fit clmm in a bootstrapped sample.")
      }
    }

    ci_l_fixed    <- apply(fixed_effect, 2, quantile, probs = 0.025, na.rm = TRUE)
    ci_h_fixed    <- apply(fixed_effect, 2, quantile, probs = 0.975, na.rm = TRUE)

    res_pooled$fixed_effects <- cbind(res_pooled$fixed_effects, ci_l_fixed,
                                      ci_u_fixed)
    colnames(res_pooled$fixed_effects) <- c("Estimate", "Lower 95% CI",
                                            "Upper 95% CI",
                                            "Lower 95% CI (Bootstrap)",
                                            "Upper 95% CI (Bootstrap)")

    ci_l_random     <- apply(random_effect, 2, quantile, probs = 0.025, na.rm = TRUE)
    ci_h_random     <- apply(random_effect, 2, quantile, probs = 0.975, na.rm = TRUE)

    res_pooled$random_effects <- cbind(res_pooled$random_effects, ci_l_random,
                                       ci_u_random)
    colnames(res_pooled$random_effects) <- c("Estimate", "Lower 95% CI",
                                             "Upper 95% CI",
                                             "Lower 95% CI (Bootstrap)",
                                             "Upper 95% CI (Bootstrap)")

    ci_l_std <- apply(std_dev, 2, quantile, probs = 0.025, na.rm = TRUE)
    ci_u_std <- apply(std_dev, 2, quantile, probs = 0.975, na.rm = TRUE)

    res_pooled$random_dist <- cbind(res_pooled$random_dist, ci_l_std, ci_u_std)
    rownames(res_pooled$random_dist) <- c("Standard deviation", "MOR")
    colnames(res_pooled$random_dist) <- c("Estimate", "Lower 95% CI",
                                          "Upper 95% CI",
                                          "Lower 95% CI (Bootstrap)",
                                          "Upper 95% CI (Bootstrap)")
  }

  class(res_pooled)    <- "pooled.clmm"
  return(res_pooled)
}

clmm_wrap <- function(data,...){
  fit <- clmm(data = data,...)
  return(fit)
}
