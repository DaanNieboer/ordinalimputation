#' Function to visualize the effect of the clusters on the outcome
#'
#' @param x a clmm.pool object returned from the pool.clmm function.
#' @param xlim Limits of the x-axis. If this is not supplied a sensible default is chosen.
#' @param ylim Limits of the y-axis. If this is not supplied a sensible default is chosen.
#' @param xlab Label at the x-axis.
#' @param ylab Label at the y-axis
#' @param refline location of reference line in caterpillar plot. Default is 0.
#' @param ... Other arguments passed to the default plot function.
#' @return A plot showing the conditional modes of the clusters with the 95%
#' prediction intervals.

plot.pooled.clmm <- function(x, xlim = NULL, ylim = NULL, xlab = "Cluster effect", ylab = "Cluster", refline = 0,...){
  index <- order(x$random_effects$mode)
  x$random_effects <- x$random_effects[index, ]
  y      <- 1:nrow(x$random_effects)


  x_plot <- x$random_effects$mode

  if(is.null(xlim)){
    xmax   <- with(x$random_effects, mode + qnorm(0.975) * sqrt(cond_var))
    xmin   <- with(x$random_effects, mode - qnorm(0.975) * sqrt(cond_var))

    xlim <- c(floor(min(xmin) + 0.1), ceiling(max(xmax) + 0.1))
  }
  if(is.null(ylim)){
    ylim <- c(0, nrow(x$random_effects) + 1)
  }

  plot(x_plot, y, xlim = xlim, ylim = ylim, pch = 19, xlab = xlab, ylab = ylab,...)
  segments(x0 = xmin, x1 = xmax, y0 = y, y1 = y)
  abline(v = refline, lty = 3)
}
