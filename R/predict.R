# MSA prediction layer — prediction from a fitted model.
# Input = a fit object; each predictor evaluates a forward model function
# (R/model.R) at the fit's (a1, a2) and propagates the fit's covariance to a
# delta-method error band. @family prediction.
#
# Goodness of fit is NOT here: the fitter computes it at the optimum and attaches it as
# `fit$gof` (see gof_from_primitives() in R/fitting.R), so there is no accessor to call.

#' Validate that a fit carries the delta-method inputs
#'
#' A fit must carry the point estimate + `theta`-scale covariance the delta method
#' needs. Fails loud on a wrong-type object rather than propagating an `NA` band.
#'
#' @param fit The object to validate; must be a list carrying `a1`, `a2`, and a 2x2
#'   `cov` matrix.
#' @param producer Name of the producing function, used in the error message.
#' @return Invisibly returns `fit`; stops with a clear message if an input is missing
#'   or malformed.
#' @noRd
validate_ml_fit <- function(fit, producer) {
  needed <- c("a1", "a2", "cov")
  if (!is.list(fit) || !all(needed %in% names(fit)) ||
      !is.matrix(fit$cov) || !all(dim(fit$cov) == c(2L, 2L))) {
    stop("fit must be an ML fit carrying a1, a2, and a 2x2 cov matrix; ",
         "got an object without them. Pass a fit from ", producer, ".")
  }
  invisible(fit)
}
