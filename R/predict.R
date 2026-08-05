# MSA prediction layer — prediction from a fitted model.
# Input = a fit object (an _ml fit); each predictor evaluates a forward model
# function (calculate_*, in R/model.R) at the fit's (a1, a2) and propagates the
# fit's covariance to a delta-method error band. @family prediction.

# ---- goodness of fit -------------------------------------------------------------
# Absolute (single-model) goodness of fit for the _ml fits, glm/broom-patterned: the
# fitter stores raw primitives (logLik, deviance, null_deviance, nobs, k), these derive
# the reportable row. Comparative AIC/BIC across nested models needs those models fit
# at their own maxima (a future fix_a1/fix_a2 fitter) -- not built here.

#' Derive the glance-style goodness-of-fit row from stored primitives
#'
#' Pure arithmetic, shared by the site and mode accessors so the formulas live in one
#' place. `D2` is deviance-explained (`= 1 - Var(resid)/Var(obs)` since both deviances
#' are sums of squares on the same `n`); `AIC`/`BIC` use the full profiled Gaussian
#' `logLik`.
#'
#' @param logLik The profiled Gaussian log-likelihood at the fit.
#' @param deviance The model deviance (residual sum of squares).
#' @param null_deviance The flat/mean-only null deviance.
#' @param nobs Number of observations.
#' @param k Number of estimated parameters.
#' @return A one-row tibble: `D2`, `AIC`, `BIC`, and the primitives (`logLik`,
#'   `deviance`, `null_deviance`, `nobs`, `k`).
#' @noRd
gof_from_primitives <- function(logLik, deviance, null_deviance, nobs, k) {
  tibble::tibble(
    D2            = 1 - deviance / null_deviance,
    AIC           = -2 * logLik + 2 * k,
    BIC           = -2 * logLik + k * log(nobs),
    logLik        = logLik,
    deviance      = deviance,
    null_deviance = null_deviance,
    nobs          = nobs,
    k             = k
  )
}

#' Validate that a fit carries the goodness-of-fit primitives
#'
#' Shared validator: an `_ml` fit must carry every GoF primitive. Fails loud on a raw
#' (pre-GoF) or wrong-type object rather than returning a silent `NA`.
#'
#' @param fit The object to validate; must be a list carrying `logLik`, `deviance`,
#'   `null_deviance`, `nobs`, `k`.
#' @param producer Name of the producing function, used in the error message.
#' @return Invisibly returns `fit`; stops with a clear message if a primitive is missing.
#' @noRd
validate_gof_fit <- function(fit, producer) {
  needed <- c("logLik", "deviance", "null_deviance", "nobs", "k")
  if (!is.list(fit) || !all(needed %in% names(fit))) {
    stop("fit must be a list carrying the goodness-of-fit primitives (",
         paste(needed, collapse = ", "), "); got an object without them. ",
         "Pass a fit from ", producer, ".")
  }
  invisible(fit)
}

#' Goodness of fit of an lrmsd_i ML fit
#'
#' Absolute goodness-of-fit summary for a maximum-likelihood site fit
#' ([fit_lrmsd_i_msa_ml()]), in the style of [broom::glance()] for a `glm`. Reports
#' the deviance-explained `D2` against a flat / mean-only null, together with `AIC`
#' and `BIC` and the primitives they derive from. Because the fit mean-centres the
#' profile, `D2 = 1 - deviance/null_deviance` is exactly the fraction of profile
#' variance explained, `1 - Var(residuals)/Var(observed)`.
#'
#' `D2` is at most `1` but has **no lower bound**: it is negative whenever the
#' prediction is worse than the flat null (`Var(residuals) > Var(observed)`), which
#' happens for parameter values that inflate the predicted amplitude. It is returned
#' unclamped -- a negative `D2` is a real signal that the fit is poor, not an error.
#'
#' `AIC`/`BIC` here are per-fit numbers; they are only meaningful *compared* against
#' another model fit at its own maximum. Comparing nested variants (MS/MA) requires
#' fitting them independently, which this package does not yet do.
#'
#' @param fit A list from [fit_lrmsd_i_msa_ml()] (carrying `logLik`, `deviance`,
#'   `null_deviance`, `nobs`, `k`).
#' @return A one-row tibble: `D2` (deviance explained vs the flat null,
#'   `1 - deviance/null_deviance`), `AIC` (`-2*logLik + 2*k`), `BIC`
#'   (`-2*logLik + k*log(nobs)`), and the primitives `logLik`, `deviance`,
#'   `null_deviance`, `nobs`, `k`.
#' @seealso [fit_lrmsd_i_msa_ml()] (produces the fit), [gof_lrmsd_n_msa_ml()] (mode
#'   counterpart).
#' @family prediction
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' ml <- fit_lrmsd_i_msa_ml(spm, znb_profile$pdb_site, znb_profile$lrmsd_i_obs)
#' gof_lrmsd_i_msa_ml(ml)
#' }
#' @export
gof_lrmsd_i_msa_ml <- function(fit) {
  validate_gof_fit(fit, "fit_lrmsd_i_msa_ml()")
  gof_from_primitives(fit$logLik, fit$deviance, fit$null_deviance, fit$nobs, fit$k)
}

#' Goodness of fit of an lrmsd_n ML fit (mode form)
#'
#' Mode counterpart of [gof_lrmsd_i_msa_ml()]: absolute goodness-of-fit summary for a
#' maximum-likelihood mode fit ([fit_lrmsd_n_msa_ml()]). Identical machinery and
#' output; see [gof_lrmsd_i_msa_ml()] for the interpretation of `D2` (unbounded below,
#' capped at `1`) and the caveat on per-fit `AIC`/`BIC`.
#'
#' @param fit A list from [fit_lrmsd_n_msa_ml()] (carrying `logLik`, `deviance`,
#'   `null_deviance`, `nobs`, `k`).
#' @return A one-row tibble with the same columns as [gof_lrmsd_i_msa_ml()]: `D2`,
#'   `AIC`, `BIC`, `logLik`, `deviance`, `null_deviance`, `nobs`, `k`.
#' @seealso [fit_lrmsd_n_msa_ml()] (produces the fit), [gof_lrmsd_i_msa_ml()] (site
#'   counterpart).
#' @family prediction
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' ml <- fit_lrmsd_n_msa_ml(spm, znb_profile_n$n, znb_profile_n$lrmsd_n_obs)
#' gof_lrmsd_n_msa_ml(ml)
#' }
#' @export
gof_lrmsd_n_msa_ml <- function(fit) {
  validate_gof_fit(fit, "fit_lrmsd_n_msa_ml()")
  gof_from_primitives(fit$logLik, fit$deviance, fit$null_deviance, fit$nobs, fit$k)
}

#' Validate that a fit carries the delta-method inputs
#'
#' An `_ml` fit must carry the point estimate + `theta`-scale covariance the delta method
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
