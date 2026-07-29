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
#' ml <- fit_lrmsd_i_msa_ml(spm, znb_profile)
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
#' ml <- fit_lrmsd_n_msa_ml(spm, znb_profile_n)
#' gof_lrmsd_n_msa_ml(ml)
#' }
#' @export
gof_lrmsd_n_msa_ml <- function(fit) {
  validate_gof_fit(fit, "fit_lrmsd_n_msa_ml()")
  gof_from_primitives(fit$logLik, fit$deviance, fit$null_deviance, fit$nobs, fit$k)
}

#' Validate that a fit carries the delta-method inputs
#'
#' An `_ml` fit must carry the point estimate + `t`-scale covariance the delta method
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

#' Predicted lrmsd_i profile with delta-method error bands from an ML fit
#'
#' Propagates the uncertainty of a point ML fit ([fit_lrmsd_i_msa_ml()]) to the
#' predicted per-site divergence profile, returning the point profile and a symmetric
#' error band at every residue. The band is by the delta method -- deterministic and
#' Gaussian (symmetric), a frequentist confidence band consistent with the fit's own
#' standard errors. This is the band on the *uncentred* `lrmsd`, holding the model's
#' own predicted level. For the mean-centred quantity that compares to observations,
#' use [predict_nlrmsd_i_msa_ml()].
#'
#' The band combines two uncertainty sources: the `(a1, a2)` parameter uncertainty
#' (delta method, above) and the SPM sampling error (each `lrmsd` is a mean over a
#' *finite* mutant scan, not the infinite limit). `uncertainty` selects which are
#' included; the total variance is their sum under an independence assumption (the
#' cross-term is dropped -- the two arms share the SPM table, so their covariance is
#' plausibly positive and the band may then very slightly under-cover). The column
#' schema is the same in every mode: `"none"` returns a zero-width band, not fewer
#' columns.
#'
#' @param fit A list from [fit_lrmsd_i_msa_ml()] (carrying `a1`, `a2`, and the 2x2
#'   `cov` on the `(a1, log2(a2+1))` scale).
#' @param spm The `spm` object from [generate_spm_data()] (the same one used for the fit).
#' @param level Confidence-band coverage (default 0.95).
#' @param uncertainty Which error sources enter the band: `"both"` (default, parameter
#'   + SPM sampling), `"spm"` (SPM sampling only), `"parameter"` (`(a1, a2)` only), or
#'   `"none"` (zero-width). The output columns are identical in every mode.
#' @return A tibble with one row per model residue: `i`, `pdb_site`, and the profile
#'   summary `lrmsd_i_msa_mean`, `lrmsd_i_msa_lower`, `lrmsd_i_msa_upper`.
#' @seealso [predict_nlrmsd_i_msa_ml()] (the mean-centred counterpart);
#'   [fit_lrmsd_i_msa_ml()] (produces the fit); [calculate_lrmsd_i_msa()] (the forward
#'   map linearized).
#' @family prediction
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' ml <- fit_lrmsd_i_msa_ml(spm, znb_profile)
#' head(predict_lrmsd_i_msa_ml(ml, spm))
#' }
#' @export
predict_lrmsd_i_msa_ml <- function(fit, spm, level = 0.95,
                                   uncertainty = c("both", "spm", "parameter", "none")) {
  validate_ml_fit(fit, "fit_lrmsd_i_msa_ml()")
  use <- uncertainty_gates(uncertainty, level)

  keys <- calculate_lrmsd_i_msa(spm, fit$a1, fit$a2) %>%
    dplyr::left_join(spm$site_map, by = "i")
  mean_v <- keys$lrmsd_i_msa
  vp <- if (use$param) {
    f <- function(t) calculate_lrmsd_i_msa(spm, t[1], 2^t[2] - 1)$lrmsd_i_msa
    var_param_delta(f, ml_t_hat(fit), fit$cov, centred = FALSE)
  } else 0
  vs <- if (use$spm) {
    w <- weights_jm_spm(spm, fit$a1, fit$a2)
    var_spm_lrmsd_i_msa(spm$dr2_ijm, w)
  } else 0
  band <- delta_band(mean_v, vp + vs, level, "lrmsd_i_msa")
  dplyr::bind_cols(tibble(i = keys$i, pdb_site = keys$pdb_site), band)
}

#' Predicted mean-centred nlrmsd_i profile with delta-method error bands from an ML fit
#'
#' Mean-centred counterpart of [predict_lrmsd_i_msa_ml()]: the per-site profile
#' centred by its own mean over the **full model support** (all model residues),
#' `nlrmsd_i = lrmsd_i - mean_S(lrmsd)`, with a symmetric delta-method band. Because
#' `mean_S(lrmsd)` is itself a function of `(a1, a2)`, the band uses the column-centred
#' gradient `g_i - mean_S(g)`; it is therefore **narrower** than (not a vertical shift
#' of) the [predict_lrmsd_i_msa_ml()] band. This is the quantity the ML fit is on: the
#' likelihood centres both sides, so `nlrmsd` is what the fit has information about.
#'
#' Support note: prediction centres over all model residues, deliberately agnostic to
#' which residues a given dataset happens to observe. The fit itself centres over the
#' observation-matched residues (its `fitted.values`), so predicted and fitted `nlrmsd`
#' sit at slightly different levels **by design**. When overlaying observations, centre
#' them on their own matched support -- that is the comparison being made.
#'
#' Like [predict_lrmsd_i_msa_ml()], the band combines the `(a1, a2)` parameter arm and
#' the SPM sampling arm (summed under independence), selected by `uncertainty`; the
#' column schema is stable across modes.
#'
#' @param fit A list from [fit_lrmsd_i_msa_ml()] (carrying `a1`, `a2`, and the 2x2
#'   `cov` on the `(a1, log2(a2+1))` scale).
#' @param spm The `spm` object from [generate_spm_data()] (the same one used for the fit).
#' @param level Confidence-band coverage (default 0.95).
#' @param uncertainty Which error sources enter the band: `"both"` (default), `"spm"`,
#'   `"parameter"`, or `"none"`. See [predict_lrmsd_i_msa_ml()].
#' @return A tibble with one row per model residue: `i`, `pdb_site`, and the centred
#'   profile summary `nlrmsd_i_msa_mean`, `nlrmsd_i_msa_lower`, `nlrmsd_i_msa_upper`.
#' @seealso [predict_lrmsd_i_msa_ml()] (the uncentred counterpart);
#'   [fit_lrmsd_i_msa_ml()] (produces the fit); [calculate_lrmsd_i_msa()] (the forward
#'   map linearized).
#' @family prediction
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' ml <- fit_lrmsd_i_msa_ml(spm, znb_profile)
#' head(predict_nlrmsd_i_msa_ml(ml, spm))
#' }
#' @export
predict_nlrmsd_i_msa_ml <- function(fit, spm, level = 0.95,
                                    uncertainty = c("both", "spm", "parameter", "none")) {
  validate_ml_fit(fit, "fit_lrmsd_i_msa_ml()")
  use <- uncertainty_gates(uncertainty, level)

  keys <- calculate_nlrmsd_i_msa(spm, fit$a1, fit$a2) %>%
    dplyr::left_join(spm$site_map, by = "i")
  mean_v <- keys$nlrmsd_i_msa                          # centred over full model support
  vp <- if (use$param) {
    f <- function(t) calculate_lrmsd_i_msa(spm, t[1], 2^t[2] - 1)$lrmsd_i_msa
    var_param_delta(f, ml_t_hat(fit), fit$cov, centred = TRUE)
  } else 0
  vs <- if (use$spm) {
    w <- weights_jm_spm(spm, fit$a1, fit$a2)
    var_spm_nlrmsd_i_msa(spm$dr2_ijm, w)
  } else 0
  band <- delta_band(mean_v, vp + vs, level, "nlrmsd_i_msa")
  dplyr::bind_cols(tibble(i = keys$i, pdb_site = keys$pdb_site), band)
}

#' Predicted lrmsd_n profile with delta-method error bands from an ML fit (mode form)
#'
#' Mode counterpart of [predict_lrmsd_i_msa_ml()]: propagates a mode ML fit
#' ([fit_lrmsd_n_msa_ml()]) to the predicted per-mode divergence profile with a
#' symmetric delta-method band on the *uncentred* `lrmsd`. Modes are not
#' residue-anchored, so there is no `pdb_site` column. For the mean-centred quantity,
#' use [predict_nlrmsd_n_msa_ml()].
#'
#' Like the site form, the band combines the `(a1, a2)` parameter arm and the SPM
#' sampling arm (summed under independence), selected by `uncertainty`.
#'
#' @param fit A list from [fit_lrmsd_n_msa_ml()] (carrying `a1`, `a2`, `cov`).
#' @param spm The `spm` object from [generate_spm_data()] (the same one used for the fit).
#' @param level Confidence-band coverage (default 0.95).
#' @param uncertainty Which error sources enter the band: `"both"` (default), `"spm"`,
#'   `"parameter"`, or `"none"`. See [predict_lrmsd_i_msa_ml()].
#' @return A tibble with one row per mode: `n`, and the profile summary
#'   `lrmsd_n_msa_mean`, `lrmsd_n_msa_lower`, `lrmsd_n_msa_upper`.
#' @seealso [predict_nlrmsd_n_msa_ml()] (the mean-centred counterpart);
#'   [fit_lrmsd_n_msa_ml()] (produces the fit); [predict_lrmsd_i_msa_ml()]
#'   (site counterpart); [calculate_lrmsd_n_msa()] (the forward map linearized).
#' @family prediction
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' ml <- fit_lrmsd_n_msa_ml(spm, znb_profile_n)
#' head(predict_lrmsd_n_msa_ml(ml, spm))
#' }
#' @export
predict_lrmsd_n_msa_ml <- function(fit, spm, level = 0.95,
                                   uncertainty = c("both", "spm", "parameter", "none")) {
  validate_ml_fit(fit, "fit_lrmsd_n_msa_ml()")
  use <- uncertainty_gates(uncertainty, level)

  keys <- calculate_lrmsd_n_msa(spm, fit$a1, fit$a2)
  mean_v <- keys$lrmsd_n_msa
  vp <- if (use$param) {
    f <- function(t) calculate_lrmsd_n_msa(spm, t[1], 2^t[2] - 1)$lrmsd_n_msa
    var_param_delta(f, ml_t_hat(fit), fit$cov, centred = FALSE)
  } else 0
  vs <- if (use$spm) {
    w <- weights_jm_spm(spm, fit$a1, fit$a2)
    var_spm_lrmsd_n_msa(spm$dr2_njm, w)
  } else 0
  band <- delta_band(mean_v, vp + vs, level, "lrmsd_n_msa")
  dplyr::bind_cols(tibble(n = keys$n), band)
}

#' Predicted mean-centred nlrmsd_n profile with delta-method error bands (mode form)
#'
#' Mode counterpart of [predict_nlrmsd_i_msa_ml()]: the per-mode profile centred by its
#' own mean over the full model support (all modes), with a symmetric delta-method band
#' using the column-centred gradient. No `pdb_site` (modes are not residue-anchored).
#' See [predict_nlrmsd_i_msa_ml()] for the support note (predict centres over all modes;
#' the fit centres over the observation-matched modes, so the two levels differ by
#' design).
#'
#' Like the site form, the band combines the `(a1, a2)` parameter arm and the SPM
#' sampling arm (summed under independence), selected by `uncertainty`.
#'
#' @param fit A list from [fit_lrmsd_n_msa_ml()] (carrying `a1`, `a2`, `cov`).
#' @param spm The `spm` object from [generate_spm_data()] (the same one used for the fit).
#' @param level Confidence-band coverage (default 0.95).
#' @param uncertainty Which error sources enter the band: `"both"` (default), `"spm"`,
#'   `"parameter"`, or `"none"`. See [predict_lrmsd_i_msa_ml()].
#' @return A tibble with one row per mode: `n`, and the centred profile summary
#'   `nlrmsd_n_msa_mean`, `nlrmsd_n_msa_lower`, `nlrmsd_n_msa_upper`.
#' @seealso [predict_lrmsd_n_msa_ml()] (the uncentred counterpart);
#'   [predict_nlrmsd_i_msa_ml()] (site counterpart); [fit_lrmsd_n_msa_ml()] (produces
#'   the fit); [calculate_lrmsd_n_msa()] (the forward map linearized).
#' @family prediction
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' ml <- fit_lrmsd_n_msa_ml(spm, znb_profile_n)
#' head(predict_nlrmsd_n_msa_ml(ml, spm))
#' }
#' @export
predict_nlrmsd_n_msa_ml <- function(fit, spm, level = 0.95,
                                    uncertainty = c("both", "spm", "parameter", "none")) {
  validate_ml_fit(fit, "fit_lrmsd_n_msa_ml()")
  use <- uncertainty_gates(uncertainty, level)

  keys <- calculate_nlrmsd_n_msa(spm, fit$a1, fit$a2)
  mean_v <- keys$nlrmsd_n_msa
  vp <- if (use$param) {
    f <- function(t) calculate_lrmsd_n_msa(spm, t[1], 2^t[2] - 1)$lrmsd_n_msa
    var_param_delta(f, ml_t_hat(fit), fit$cov, centred = TRUE)
  } else 0
  vs <- if (use$spm) {
    w <- weights_jm_spm(spm, fit$a1, fit$a2)
    var_spm_nlrmsd_n_msa(spm$dr2_njm, w)
  } else 0
  band <- delta_band(mean_v, vp + vs, level, "nlrmsd_n_msa")
  dplyr::bind_cols(tibble(n = keys$n), band)
}

#' Nested-model divergence profiles with delta-method bands from an ML fit
#'
#' The four nested variants (MM, MS, MA, MSA) of the predicted per-site *uncentred*
#' profile, each with a symmetric delta-method band, from a point ML fit
#' ([fit_lrmsd_i_msa_ml()]). Each variant's band linearizes
#' [calculate_lrmsd_i_nested_models()] about the fit's `(a1, a2)`; MM `(0, 0)` does not
#' depend on the parameters, so its *parameter-uncertainty* band is zero-width -- but its
#' SPM-sampling band is not, so MM has a nonzero total band whenever `uncertainty`
#' includes the SPM arm. For the mean-centred variants that compare to observations, use
#' [predict_nlrmsd_i_nested_models_ml()].
#'
#' Each variant carries its own selection strengths (MM `(0,0)`, MS `(a1,0)`, MA
#' `(0,a2)`, MSA `(a1,a2)`), so its SPM sampling weights -- and therefore its SPM band --
#' are computed at that variant's own point.
#'
#' @param fit A list from [fit_lrmsd_i_msa_ml()] (carrying `a1`, `a2`, `cov`).
#' @param spm The `spm` object from [generate_spm_data()] (the same one used for the fit).
#' @param level Confidence-band coverage (default 0.95).
#' @param uncertainty Which error sources enter the band: `"both"` (default), `"spm"`,
#'   `"parameter"`, or `"none"`. See [predict_lrmsd_i_msa_ml()].
#' @return A tibble with one row per residue: `i`, `pdb_site`, and a mean/lower/upper
#'   triple for each variant -- `lrmsd_i_mm_{mean,lower,upper}`,
#'   `lrmsd_i_ms_{mean,lower,upper}`, `lrmsd_i_ma_{mean,lower,upper}`,
#'   `lrmsd_i_msa_{mean,lower,upper}` (12 band columns).
#' @seealso [predict_nlrmsd_i_nested_models_ml()] (mean-centred counterpart);
#'   [fit_lrmsd_i_msa_ml()] (produces the fit);
#'   [calculate_lrmsd_i_nested_models()] (the forward map).
#' @family prediction
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' ml <- fit_lrmsd_i_msa_ml(spm, znb_profile)
#' head(predict_lrmsd_i_nested_models_ml(ml, spm))
#' }
#' @export
predict_lrmsd_i_nested_models_ml <- function(fit, spm, level = 0.95,
                                             uncertainty = c("both", "spm", "parameter", "none")) {
  validate_ml_fit(fit, "fit_lrmsd_i_msa_ml()")
  use   <- uncertainty_gates(uncertainty, level)
  keys  <- calculate_lrmsd_i_nested_models(spm, fit$a1, fit$a2)
  t_hat <- ml_t_hat(fit)
  # variant -> (a1, a2) for the per-variant SPM weights.
  variant_a <- list(mm = c(0, 0), ms = c(fit$a1, 0), ma = c(0, fit$a2), msa = c(fit$a1, fit$a2))

  bands <- lapply(names(variant_a), function(v) {
    col <- paste0("lrmsd_i_", v)
    mean_v <- keys[[col]]
    vp <- if (use$param) {
      f <- function(t) calculate_lrmsd_i_nested_models(spm, t[1], 2^t[2] - 1)[[col]]
      var_param_delta(f, t_hat, fit$cov, centred = FALSE)
    } else 0
    vs <- if (use$spm) {
      a <- variant_a[[v]]
      w <- weights_jm_spm(spm, a[1], a[2])
      var_spm_lrmsd_i_msa(spm$dr2_ijm, w)
    } else 0
    delta_band(mean_v, vp + vs, level, col)
  })
  dplyr::bind_cols(keys[c("i", "pdb_site")], bands)
}

#' Mean-centred nested-model divergence profiles with delta-method bands from an ML fit
#'
#' Mean-centred counterpart of [predict_lrmsd_i_nested_models_ml()]: the four nested
#' variants (MM, MS, MA, MSA), each centred by *its own* mean over the full model
#' support and banded with the column-centred gradient. MM's gradient is constant in
#' `(a1, a2)`, so its centred *parameter* band is exactly zero -- but its SPM-sampling
#' band is not, so MM has a nonzero total band whenever `uncertainty` includes the SPM
#' arm. See [predict_nlrmsd_i_msa_ml()] for the support note.
#'
#' Each variant's SPM band is computed at its own selection strengths (MM `(0,0)`, MS
#' `(a1,0)`, MA `(0,a2)`, MSA `(a1,a2)`).
#'
#' @param fit A list from [fit_lrmsd_i_msa_ml()] (carrying `a1`, `a2`, `cov`).
#' @param spm The `spm` object from [generate_spm_data()] (the same one used for the fit).
#' @param level Confidence-band coverage (default 0.95).
#' @param uncertainty Which error sources enter the band: `"both"` (default), `"spm"`,
#'   `"parameter"`, or `"none"`. See [predict_lrmsd_i_msa_ml()].
#' @return A tibble with one row per residue: `i`, `pdb_site`, and a mean/lower/upper
#'   triple for each centred variant -- `nlrmsd_i_mm_{mean,lower,upper}`,
#'   `nlrmsd_i_ms_{mean,lower,upper}`, `nlrmsd_i_ma_{mean,lower,upper}`,
#'   `nlrmsd_i_msa_{mean,lower,upper}` (12 band columns).
#' @seealso [predict_lrmsd_i_nested_models_ml()] (uncentred counterpart);
#'   [fit_lrmsd_i_msa_ml()] (produces the fit);
#'   [calculate_lrmsd_i_nested_models()] (the forward map).
#' @family prediction
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' ml <- fit_lrmsd_i_msa_ml(spm, znb_profile)
#' head(predict_nlrmsd_i_nested_models_ml(ml, spm))
#' }
#' @export
predict_nlrmsd_i_nested_models_ml <- function(fit, spm, level = 0.95,
                                              uncertainty = c("both", "spm", "parameter", "none")) {
  validate_ml_fit(fit, "fit_lrmsd_i_msa_ml()")
  use   <- uncertainty_gates(uncertainty, level)
  keys  <- calculate_nlrmsd_i_nested_models(spm, fit$a1, fit$a2)
  t_hat <- ml_t_hat(fit)
  variant_a <- list(mm = c(0, 0), ms = c(fit$a1, 0), ma = c(0, fit$a2), msa = c(fit$a1, fit$a2))

  bands <- lapply(names(variant_a), function(v) {
    mean_v <- keys[[paste0("nlrmsd_i_", v)]]
    vp <- if (use$param) {
      f <- function(t) calculate_lrmsd_i_nested_models(spm, t[1], 2^t[2] - 1)[[paste0("lrmsd_i_", v)]]
      var_param_delta(f, t_hat, fit$cov, centred = TRUE)
    } else 0
    vs <- if (use$spm) {
      a <- variant_a[[v]]
      w <- weights_jm_spm(spm, a[1], a[2])
      var_spm_nlrmsd_i_msa(spm$dr2_ijm, w)
    } else 0
    delta_band(mean_v, vp + vs, level, paste0("nlrmsd_i_", v))
  })
  dplyr::bind_cols(keys[c("i", "pdb_site")], bands)
}

#' Nested-model divergence profiles with delta-method bands from an ML fit (mode form)
#'
#' Mode counterpart of [predict_lrmsd_i_nested_models_ml()]: the four nested variants
#' of the predicted per-mode *uncentred* profile, each with a symmetric delta-method
#' band, from a mode ML fit ([fit_lrmsd_n_msa_ml()]). No `pdb_site` (modes are not
#' residue-anchored). For the mean-centred variants use
#' [predict_nlrmsd_n_nested_models_ml()].
#'
#' @param fit A list from [fit_lrmsd_n_msa_ml()] (carrying `a1`, `a2`, `cov`).
#' @param spm The `spm` object from [generate_spm_data()] (the same one used for the fit).
#' @param level Confidence-band coverage (default 0.95).
#' @return A tibble with one row per mode: `n` and a mean/lower/upper triple for each
#'   variant -- `lrmsd_n_mm_{mean,lower,upper}`, `lrmsd_n_ms_{mean,lower,upper}`,
#'   `lrmsd_n_ma_{mean,lower,upper}`, `lrmsd_n_msa_{mean,lower,upper}` (12 band columns).
#' @seealso [predict_nlrmsd_n_nested_models_ml()] (mean-centred counterpart);
#'   [predict_lrmsd_i_nested_models_ml()] (site counterpart);
#'   [calculate_lrmsd_n_nested_models()] (the forward map).
#' @family prediction
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' ml <- fit_lrmsd_n_msa_ml(spm, znb_profile_n)
#' head(predict_lrmsd_n_nested_models_ml(ml, spm))
#' }
#' @param uncertainty Which error sources enter the band: `"both"` (default), `"spm"`,
#'   `"parameter"`, or `"none"`. See [predict_lrmsd_i_msa_ml()].
#' @export
predict_lrmsd_n_nested_models_ml <- function(fit, spm, level = 0.95,
                                             uncertainty = c("both", "spm", "parameter", "none")) {
  validate_ml_fit(fit, "fit_lrmsd_n_msa_ml()")
  use   <- uncertainty_gates(uncertainty, level)
  keys  <- calculate_lrmsd_n_nested_models(spm, fit$a1, fit$a2)
  t_hat <- ml_t_hat(fit)
  variant_a <- list(mm = c(0, 0), ms = c(fit$a1, 0), ma = c(0, fit$a2), msa = c(fit$a1, fit$a2))

  bands <- lapply(names(variant_a), function(v) {
    col    <- paste0("lrmsd_n_", v)
    mean_v <- keys[[col]]
    vp <- if (use$param) {
      f <- function(t) calculate_lrmsd_n_nested_models(spm, t[1], 2^t[2] - 1)[[col]]
      var_param_delta(f, t_hat, fit$cov, centred = FALSE)
    } else 0
    vs <- if (use$spm) {
      a <- variant_a[[v]]
      w <- weights_jm_spm(spm, a[1], a[2])
      var_spm_lrmsd_n_msa(spm$dr2_njm, w)
    } else 0
    delta_band(mean_v, vp + vs, level, col)
  })
  dplyr::bind_cols(keys["n"], bands)
}

#' Mean-centred nested-model divergence profiles with delta-method bands (mode form)
#'
#' Mode counterpart of [predict_nlrmsd_i_nested_models_ml()]: the four nested variants
#' of the per-mode profile, each centred by its own mean over the full model support
#' and banded with the column-centred gradient. MM's centred *parameter* band is
#' zero-width, but its SPM-sampling band is not (so MM's total band is nonzero when
#' `uncertainty` includes the SPM arm). No `pdb_site` (modes are not residue-anchored).
#'
#' @param fit A list from [fit_lrmsd_n_msa_ml()] (carrying `a1`, `a2`, `cov`).
#' @param spm The `spm` object from [generate_spm_data()] (the same one used for the fit).
#' @param level Confidence-band coverage (default 0.95).
#' @param uncertainty Which error sources enter the band: `"both"` (default), `"spm"`,
#'   `"parameter"`, or `"none"`. See [predict_lrmsd_i_msa_ml()].
#' @return A tibble with one row per mode: `n` and a mean/lower/upper triple for each
#'   centred variant -- `nlrmsd_n_mm_{mean,lower,upper}`, `nlrmsd_n_ms_{mean,lower,upper}`,
#'   `nlrmsd_n_ma_{mean,lower,upper}`, `nlrmsd_n_msa_{mean,lower,upper}` (12 band columns).
#' @seealso [predict_lrmsd_n_nested_models_ml()] (uncentred counterpart);
#'   [predict_nlrmsd_i_nested_models_ml()] (site counterpart);
#'   [calculate_lrmsd_n_nested_models()] (the forward map).
#' @family prediction
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' ml <- fit_lrmsd_n_msa_ml(spm, znb_profile_n)
#' head(predict_nlrmsd_n_nested_models_ml(ml, spm))
#' }
#' @export
predict_nlrmsd_n_nested_models_ml <- function(fit, spm, level = 0.95,
                                              uncertainty = c("both", "spm", "parameter", "none")) {
  validate_ml_fit(fit, "fit_lrmsd_n_msa_ml()")
  use   <- uncertainty_gates(uncertainty, level)
  keys  <- calculate_nlrmsd_n_nested_models(spm, fit$a1, fit$a2)
  t_hat <- ml_t_hat(fit)
  variant_a <- list(mm = c(0, 0), ms = c(fit$a1, 0), ma = c(0, fit$a2), msa = c(fit$a1, fit$a2))

  bands <- lapply(names(variant_a), function(v) {
    mean_v <- keys[[paste0("nlrmsd_n_", v)]]
    vp <- if (use$param) {
      f <- function(t) calculate_lrmsd_n_nested_models(spm, t[1], 2^t[2] - 1)[[paste0("lrmsd_n_", v)]]
      var_param_delta(f, t_hat, fit$cov, centred = TRUE)
    } else 0
    vs <- if (use$spm) {
      a <- variant_a[[v]]
      w <- weights_jm_spm(spm, a[1], a[2])
      var_spm_nlrmsd_n_msa(spm$dr2_njm, w)
    } else 0
    delta_band(mean_v, vp + vs, level, paste0("nlrmsd_n_", v))
  })
  dplyr::bind_cols(keys["n"], bands)
}

#' Decomposition of nlrmsd_i_msa with delta-method bands from an ML fit
#'
#' The three contributions (`nphi_mut`, `nphi_stab`, `nphi_act`) of the *mean-centred*
#' per-site profile `nlrmsd_i_msa`, each with a symmetric delta-method band, from a
#' point ML fit ([fit_lrmsd_i_msa_ml()]). Each contribution is centred by its own mean
#' over the full model support and banded with the column-centred gradient (see
#' [predict_nlrmsd_i_msa_ml()] for the support note). The decomposition is defined only
#' for the centred quantity; the three `*_mean` columns sum to `nlrmsd_i_msa` (centring
#' is linear), but the bands do not sum (SEs are not additive). `nphi_mut` (the MM term)
#' has a constant *parameter* gradient, so its parameter band is zero-width -- but its
#' SPM-sampling band is not (it equals the `nlrmsd_mm` SPM band, since `nphi_mut = nlrmsd_mm`).
#'
#' The SPM arm of each contribution is a *difference* of two nested models evaluated on
#' the **same** mutant scan (`nphi_stab = nlrmsd_ms - nlrmsd_mm`,
#' `nphi_act = nlrmsd_msa - nlrmsd_ms`), so its sampling variance is computed from the
#' per-mutant contributions of both models jointly -- retaining the between-model
#' covariance, not summing the two variances.
#'
#' @param fit A list from [fit_lrmsd_i_msa_ml()] (carrying `a1`, `a2`, `cov`).
#' @param spm The `spm` object from [generate_spm_data()] (the same one used for the fit).
#' @param level Confidence-band coverage (default 0.95).
#' @param uncertainty Which error sources enter the band: `"both"` (default), `"spm"`,
#'   `"parameter"`, or `"none"`. See [predict_lrmsd_i_msa_ml()].
#' @return A tibble with one row per residue: `i`, `pdb_site`, and a mean/lower/upper
#'   triple for each contribution -- `nphi_mut_{mean,lower,upper}`,
#'   `nphi_stab_{mean,lower,upper}`, `nphi_act_{mean,lower,upper}` (9 band columns).
#' @seealso [predict_nlrmsd_i_msa_ml()] (the profile it decomposes);
#'   [fit_lrmsd_i_msa_ml()] (produces the fit);
#'   [calculate_decomposition_i_msa()] (the forward map).
#' @family prediction
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' ml <- fit_lrmsd_i_msa_ml(spm, znb_profile)
#' head(predict_nlrmsd_i_msa_decomposition_ml(ml, spm))
#' }
#' @export
predict_nlrmsd_i_msa_decomposition_ml <- function(fit, spm, level = 0.95,
                                                  uncertainty = c("both", "spm", "parameter", "none")) {
  validate_ml_fit(fit, "fit_lrmsd_i_msa_ml()")
  use   <- uncertainty_gates(uncertainty, level)
  keys  <- calculate_nlrmsd_i_msa_decomposition(spm, fit$a1, fit$a2)
  t_hat <- ml_t_hat(fit)

  # SPM arm: each phi is a linear contrast of nested models on the SAME scan; the helper
  # builds the differenced per-mutant contributions (shared-sample covariance kept).
  spm_var <- if (use$spm) var_spm_nphi_i_msa(spm, fit$a1, fit$a2) else NULL

  bands <- lapply(c("mut", "stab", "act"), function(v) {
    mean_v <- keys[[paste0("nphi_", v)]]
    vp <- if (use$param) {
      f <- function(t) calculate_decomposition_i_msa(spm, t[1], 2^t[2] - 1)[[paste0("phi_", v)]]
      var_param_delta(f, t_hat, fit$cov, centred = TRUE)
    } else 0
    vs <- if (use$spm) spm_var[[v]] else 0
    delta_band(mean_v, vp + vs, level, paste0("nphi_", v))
  })
  dplyr::bind_cols(keys[c("i", "pdb_site")], bands)
}

#' Decomposition of nlrmsd_n_msa with delta-method bands from an ML fit (mode form)
#'
#' Mode counterpart of [predict_nlrmsd_i_msa_decomposition_ml()]: the three
#' contributions of the mean-centred per-mode profile `nlrmsd_n_msa`, each with a
#' symmetric delta-method band. No `pdb_site` (modes are not residue-anchored). The
#' three `*_mean` columns sum to `nlrmsd_n_msa`; the bands do not sum. `nphi_mut`'s
#' *parameter* band is zero-width, but its SPM-sampling band is not. See
#' [predict_nlrmsd_i_msa_decomposition_ml()] for the shared-scan SPM contrast note.
#'
#' @param fit A list from [fit_lrmsd_n_msa_ml()] (carrying `a1`, `a2`, `cov`).
#' @param spm The `spm` object from [generate_spm_data()] (the same one used for the fit).
#' @param level Confidence-band coverage (default 0.95).
#' @param uncertainty Which error sources enter the band: `"both"` (default), `"spm"`,
#'   `"parameter"`, or `"none"`. See [predict_lrmsd_i_msa_ml()].
#' @return A tibble with one row per mode: `n` and a mean/lower/upper triple for each
#'   contribution -- `nphi_mut_{mean,lower,upper}`, `nphi_stab_{mean,lower,upper}`,
#'   `nphi_act_{mean,lower,upper}` (9 band columns).
#' @seealso [predict_nlrmsd_i_msa_decomposition_ml()] (site counterpart);
#'   [predict_nlrmsd_n_msa_ml()] (the profile it decomposes);
#'   [calculate_decomposition_n_msa()] (the forward map).
#' @family prediction
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' ml <- fit_lrmsd_n_msa_ml(spm, znb_profile_n)
#' head(predict_nlrmsd_n_msa_decomposition_ml(ml, spm))
#' }
#' @export
predict_nlrmsd_n_msa_decomposition_ml <- function(fit, spm, level = 0.95,
                                                  uncertainty = c("both", "spm", "parameter", "none")) {
  validate_ml_fit(fit, "fit_lrmsd_n_msa_ml()")
  use   <- uncertainty_gates(uncertainty, level)
  keys  <- calculate_nlrmsd_n_msa_decomposition(spm, fit$a1, fit$a2)
  t_hat <- ml_t_hat(fit)

  spm_var <- if (use$spm) var_spm_nphi_n_msa(spm, fit$a1, fit$a2) else NULL

  bands <- lapply(c("mut", "stab", "act"), function(v) {
    mean_v <- keys[[paste0("nphi_", v)]]
    vp <- if (use$param) {
      f <- function(t) calculate_decomposition_n_msa(spm, t[1], 2^t[2] - 1)[[paste0("phi_", v)]]
      var_param_delta(f, t_hat, fit$cov, centred = TRUE)
    } else 0
    vs <- if (use$spm) spm_var[[v]] else 0
    delta_band(mean_v, vp + vs, level, paste0("nphi_", v))
  })
  dplyr::bind_cols(keys["n"], bands)
}
