#' Function to visualize the effect of the clusters on the outcome
#'
#' @param x a clmm.pool object returned from the pool.clmm function.
#'
#' @return A plot showing the conditional modes of the clusters with the 95%
#' prediction intervals.

plot.pooled.clmm <- function(x,...){
  index <- order(x$random_effects$mode)
  x$random_effects <- x$random_effects[index, ]
  y      <- 1:nrow(x$random_effects)


  x_plot <- x$random_effects$mode
  xmax   <- with(x$random_effects, mode + qnorm(0.975) * sqrt(cond_var))
  xmin   <- with(x$random_effects, mode - qnorm(0.975) * sqrt(cond_var))

  xlim <- c(floor(min(xmin) + 0.1), ceiling(max(xmax) + 0.1))
  ylim <- c(0, nrow(x$random_effects) + 1)

  plot(x_plot, y, xlim = xlim, ylim = ylim, pch = 19, xlab = "Center effect", ylab = "Center",
       axes = FALSE)
  axis(1)
  axis(2, at = y, las = 2, label = x$random_effects$cluster)
  segments(x0 = xmin, x1 = xmax, y0 = y, y1 = y)
  box()
  abline(v = 0, lty = 3)
}
