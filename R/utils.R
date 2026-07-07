# Generic numeric helpers shared across the package (no @family). @noRd internal.

# Weighted quantiles via midpoint cumulative weight, sorted interpolation. Small
# node sets, so this simple form is enough. Used by fit_lrmsd_i_msa_agq (fitting)
# and the predict_*_agq functions (prediction).
#' @noRd
weighted_quantile <- function(x, w, probs) {
  o  <- order(x)
  xo <- x[o]; wo <- w[o] / sum(w)
  cw <- cumsum(wo) - 0.5 * wo
  stats::approx(cw, xo, xout = probs, rule = 2)$y
}
