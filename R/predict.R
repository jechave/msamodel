# MSA prediction layer — prediction from a fitted model.
# Input = a fit object (an _ml fit); each predictor evaluates a forward model
# function (calculate_*, in R/model.R) at the fit's (a1, a2) and propagates the
# fit's covariance to a delta-method error band. @family prediction.

# ---- goodness of fit -------------------------------------------------------------
# Absolute (single-model) goodness of fit for the _ml fits, glm/broom-patterned: the
# fitter stores raw primitives (logLik, deviance, null_deviance, nobs, k), these derive
# the reportable row. Comparative AIC/BIC across nested models needs those models fit
# at their own maxima (a future fix_a1/fix_a2 fitter) -- not built here.

# Derive the glance-style GoF row from stored primitives. Pure arithmetic, shared by
# the site and mode accessors so the formulas live in one place. D2 is deviance-
# explained (= 1 - Var(resid)/Var(obs) since both deviances are sums of squares on
# the same n); AIC/BIC use the full profiled Gaussian logLik.
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

# Shared validator: an _ml fit must carry every GoF primitive. Fail loud on a raw
# (pre-GoF) or wrong-type object rather than returning a silent NA.
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
#' pp <- preprocess_spm(znb_spm)
#' ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)
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
#' pp <- preprocess_spm_mode(znb_spm)
#' ml <- fit_lrmsd_n_msa_ml(pp, znb_profile_n)
#' gof_lrmsd_n_msa_ml(ml)
#' }
#' @export
gof_lrmsd_n_msa_ml <- function(fit) {
  validate_gof_fit(fit, "fit_lrmsd_n_msa_ml()")
  gof_from_primitives(fit$logLik, fit$deviance, fit$null_deviance, fit$nobs, fit$k)
}

# ---- ML delta-method prediction ---------------------------------------------------
# Predict from an _ml fit: evaluate the forward maps (calculate_*) at the fit's point
# estimate and propagate the fit's asymptotic covariance through the forward map by the
# delta method to a symmetric error band. The fit carries `cov` on the
# t = (a1, b = log2(a2+1)) scale, so we differentiate g(t) = calc(a1=t1, a2=2^t2-1)
# w.r.t. t directly and sandwich with `cov` -- no covariance transform. Bands are
# symmetric on the reported scale (g(t_hat) +/- z*se), mirroring the fit's own SEs.

# Central-difference gradient of a vector-valued f (returns a per-site numeric vector)
# w.r.t. the 2-vector t. Returns an [nrow(f) x 2] Jacobian: column j is the partial
# derivative w.r.t. t[j]. Two f-evaluations per dimension (4 forward-map calls).
#' @noRd
grad_t <- function(f, t, h = 1e-5) {
  n <- length(f(t))
  J <- vapply(seq_along(t), function(j) {
    tp <- t; tp[j] <- tp[j] + h
    tm <- t; tm[j] <- tm[j] - h
    (f(tp) - f(tm)) / (2 * h)
  }, numeric(n))
  # vapply collapses to a bare length-length(t) vector when n == 1; force [n x length(t)].
  matrix(J, nrow = n)
}

# Delta-method mean + symmetric band of one per-site quantity. f(t) returns the
# per-site vector at parameter t; cov_t is the 2x2 covariance on the t-scale. Per-site
# variance is the sandwich diag(J cov_t J^T), computed row-wise as
# rowSums((J %*% cov_t) * J). Returns a tibble with <name>_mean/_lower/_upper.
#' @noRd
delta_band <- function(f, t_hat, cov_t, level, name) {
  mean_v <- f(t_hat)
  J      <- grad_t(f, t_hat)              # [nsite x 2]
  var_v  <- rowSums((J %*% cov_t) * J)    # diag(J cov_t J^T)
  se_v   <- sqrt(pmax(var_v, 0))          # pmax guards tiny negative round-off only
  z      <- stats::qnorm(1 - (1 - level) / 2)
  cols <- list(mean_v, mean_v - z * se_v, mean_v + z * se_v)
  names(cols) <- paste0(name, c("_mean", "_lower", "_upper"))
  tibble::as_tibble(cols)
}

# Delta-method mean + symmetric band of the MEAN-CENTRED quantity nq = q - mean_S(q),
# where S is the full model support (all rows f returns). Because mean_S(q) is itself a
# function of the parameters, the gradient of nq is the COLUMN-CENTRED Jacobian
# g_i - mean_S(g), NOT the raw g_i -- so the centred band differs from delta_band()'s,
# it is not delta_band() shifted by a constant. `name` should already carry its centred
# (n-prefixed) form; returns <name>_mean/_lower/_upper. Same cov_t and level handling as
# delta_band(); the mean subtraction and the gradient centring are one operation on one
# quantity over one support.
#' @noRd
delta_band_centred <- function(f, t_hat, cov_t, level, name) {
  v      <- f(t_hat)
  mean_v <- v - mean(v)                       # centred over the full model support
  J      <- grad_t(f, t_hat)                  # [nsite x 2]
  Jc     <- sweep(J, 2, colMeans(J))          # g_i - mean_S(g), per column
  var_v  <- rowSums((Jc %*% cov_t) * Jc)      # diag(Jc cov_t Jc^T)
  se_v   <- sqrt(pmax(var_v, 0))              # pmax guards tiny negative round-off only
  z      <- stats::qnorm(1 - (1 - level) / 2)
  cols <- list(mean_v, mean_v - z * se_v, mean_v + z * se_v)
  names(cols) <- paste0(name, c("_mean", "_lower", "_upper"))
  tibble::as_tibble(cols)
}

# An _ml fit must carry the point estimate + t-scale covariance the delta method needs.
# Fail loud on a wrong-type object rather than propagating an NA band.
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

# t = (a1, log2(a2+1)) from an _ml fit's natural-scale point estimate.
#' @noRd
ml_t_hat <- function(fit) c(fit$a1, log2(fit$a2 + 1))

#' Predicted lrmsd_i profile with delta-method error bands from an ML fit
#'
#' Propagates the asymptotic covariance of a point ML fit ([fit_lrmsd_i_msa_ml()]) to
#' the predicted per-site divergence profile, returning the point profile and a
#' symmetric error band at every residue. The band is the delta method -- the forward
#' map [calculate_lrmsd_i_msa()] linearized about the fit's `(a1, a2)`, its gradient
#' sandwiched with the fit's covariance -- so it is deterministic and Gaussian
#' (symmetric), a frequentist confidence band consistent with the fit's own standard
#' errors. This is the band on the *uncentred* `lrmsd`: the uncertainty in the model's
#' predicted `lrmsd` induced by uncertainty in `(a1, a2)`, holding the model's own
#' predicted level. For the mean-centred quantity that compares to observations, use
#' [predict_nlrmsd_i_msa_ml()].
#'
#' @param fit A list from [fit_lrmsd_i_msa_ml()] (carrying `a1`, `a2`, and the 2x2
#'   `cov` on the `(a1, log2(a2+1))` scale).
#' @param spm_pp Preprocessed data from [preprocess_spm()] (the same used for the fit).
#' @param level Confidence-band coverage (default 0.95).
#' @return A tibble with one row per model residue: `i`, `pdb_site`, and the profile
#'   summary `lrmsd_i_msa_mean`, `lrmsd_i_msa_lower`, `lrmsd_i_msa_upper`.
#' @seealso [predict_nlrmsd_i_msa_ml()] (the mean-centred counterpart);
#'   [fit_lrmsd_i_msa_ml()] (produces the fit); [calculate_lrmsd_i_msa()] (the forward
#'   map linearized).
#' @family prediction
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)
#' head(predict_lrmsd_i_msa_ml(ml, pp))
#' }
#' @export
predict_lrmsd_i_msa_ml <- function(fit, spm_pp, level = 0.95) {
  validate_ml_fit(fit, "fit_lrmsd_i_msa_ml()")
  if (length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a single number in (0, 1)")
  }
  keys <- calculate_lrmsd_i_msa(spm_pp, fit$a1, fit$a2) %>%
    dplyr::left_join(spm_pp$site_map, by = "i")
  f <- function(t) calculate_lrmsd_i_msa(spm_pp, t[1], 2^t[2] - 1)$lrmsd_i_msa
  band <- delta_band(f, ml_t_hat(fit), fit$cov, level, "lrmsd_i_msa")
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
#' @param fit A list from [fit_lrmsd_i_msa_ml()] (carrying `a1`, `a2`, and the 2x2
#'   `cov` on the `(a1, log2(a2+1))` scale).
#' @param spm_pp Preprocessed data from [preprocess_spm()] (the same used for the fit).
#' @param level Confidence-band coverage (default 0.95).
#' @return A tibble with one row per model residue: `i`, `pdb_site`, and the centred
#'   profile summary `nlrmsd_i_msa_mean`, `nlrmsd_i_msa_lower`, `nlrmsd_i_msa_upper`.
#' @seealso [predict_lrmsd_i_msa_ml()] (the uncentred counterpart);
#'   [fit_lrmsd_i_msa_ml()] (produces the fit); [calculate_lrmsd_i_msa()] (the forward
#'   map linearized).
#' @family prediction
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)
#' head(predict_nlrmsd_i_msa_ml(ml, pp))
#' }
#' @export
predict_nlrmsd_i_msa_ml <- function(fit, spm_pp, level = 0.95) {
  validate_ml_fit(fit, "fit_lrmsd_i_msa_ml()")
  if (length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a single number in (0, 1)")
  }
  keys <- calculate_lrmsd_i_msa(spm_pp, fit$a1, fit$a2) %>%
    dplyr::left_join(spm_pp$site_map, by = "i")
  f <- function(t) calculate_lrmsd_i_msa(spm_pp, t[1], 2^t[2] - 1)$lrmsd_i_msa
  band <- delta_band_centred(f, ml_t_hat(fit), fit$cov, level, "nlrmsd_i_msa")
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
#' @param fit A list from [fit_lrmsd_n_msa_ml()] (carrying `a1`, `a2`, `cov`).
#' @param spm_pp_mode Preprocessed data from [preprocess_spm_mode()].
#' @param level Confidence-band coverage (default 0.95).
#' @return A tibble with one row per mode: `n`, and the profile summary
#'   `lrmsd_n_msa_mean`, `lrmsd_n_msa_lower`, `lrmsd_n_msa_upper`.
#' @seealso [predict_nlrmsd_n_msa_ml()] (the mean-centred counterpart);
#'   [fit_lrmsd_n_msa_ml()] (produces the fit); [predict_lrmsd_i_msa_ml()]
#'   (site counterpart); [calculate_lrmsd_n_msa()] (the forward map linearized).
#' @family prediction
#' @examples
#' \dontrun{
#' pp <- preprocess_spm_mode(znb_spm)
#' ml <- fit_lrmsd_n_msa_ml(pp, znb_profile_n)
#' head(predict_lrmsd_n_msa_ml(ml, pp))
#' }
#' @export
predict_lrmsd_n_msa_ml <- function(fit, spm_pp_mode, level = 0.95) {
  validate_ml_fit(fit, "fit_lrmsd_n_msa_ml()")
  if (length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a single number in (0, 1)")
  }
  keys <- calculate_lrmsd_n_msa(spm_pp_mode, fit$a1, fit$a2)
  f <- function(t) calculate_lrmsd_n_msa(spm_pp_mode, t[1], 2^t[2] - 1)$lrmsd_n_msa
  band <- delta_band(f, ml_t_hat(fit), fit$cov, level, "lrmsd_n_msa")
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
#' @param fit A list from [fit_lrmsd_n_msa_ml()] (carrying `a1`, `a2`, `cov`).
#' @param spm_pp_mode Preprocessed data from [preprocess_spm_mode()].
#' @param level Confidence-band coverage (default 0.95).
#' @return A tibble with one row per mode: `n`, and the centred profile summary
#'   `nlrmsd_n_msa_mean`, `nlrmsd_n_msa_lower`, `nlrmsd_n_msa_upper`.
#' @seealso [predict_lrmsd_n_msa_ml()] (the uncentred counterpart);
#'   [predict_nlrmsd_i_msa_ml()] (site counterpart); [fit_lrmsd_n_msa_ml()] (produces
#'   the fit); [calculate_lrmsd_n_msa()] (the forward map linearized).
#' @family prediction
#' @examples
#' \dontrun{
#' pp <- preprocess_spm_mode(znb_spm)
#' ml <- fit_lrmsd_n_msa_ml(pp, znb_profile_n)
#' head(predict_nlrmsd_n_msa_ml(ml, pp))
#' }
#' @export
predict_nlrmsd_n_msa_ml <- function(fit, spm_pp_mode, level = 0.95) {
  validate_ml_fit(fit, "fit_lrmsd_n_msa_ml()")
  if (length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a single number in (0, 1)")
  }
  keys <- calculate_lrmsd_n_msa(spm_pp_mode, fit$a1, fit$a2)
  f <- function(t) calculate_lrmsd_n_msa(spm_pp_mode, t[1], 2^t[2] - 1)$lrmsd_n_msa
  band <- delta_band_centred(f, ml_t_hat(fit), fit$cov, level, "nlrmsd_n_msa")
  dplyr::bind_cols(tibble(n = keys$n), band)
}

#' Nested-model divergence profiles with delta-method bands from an ML fit
#'
#' The four nested variants (MM, MS, MA, MSA) of the predicted per-site *uncentred*
#' profile, each with a symmetric delta-method band, from a point ML fit
#' ([fit_lrmsd_i_msa_ml()]). Each variant's band linearizes
#' [calculate_lrmsd_i_nested_models()] about the fit's `(a1, a2)`; MM `(0, 0)` does not
#' depend on the parameters, so its (parameter-uncertainty) band is zero-width. For the
#' mean-centred variants that compare to observations, use
#' [predict_nlrmsd_i_nested_models_ml()].
#'
#' @param fit A list from [fit_lrmsd_i_msa_ml()] (carrying `a1`, `a2`, `cov`).
#' @param spm_pp Preprocessed data from [preprocess_spm()].
#' @param level Confidence-band coverage (default 0.95).
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
#' pp <- preprocess_spm(znb_spm)
#' ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)
#' head(predict_lrmsd_i_nested_models_ml(ml, pp))
#' }
#' @export
predict_lrmsd_i_nested_models_ml <- function(fit, spm_pp, level = 0.95) {
  validate_ml_fit(fit, "fit_lrmsd_i_msa_ml()")
  if (length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a single number in (0, 1)")
  }
  keys     <- calculate_lrmsd_i_nested_models(spm_pp, fit$a1, fit$a2)
  variants <- c("lrmsd_i_mm", "lrmsd_i_ms", "lrmsd_i_ma", "lrmsd_i_msa")
  t_hat    <- ml_t_hat(fit)

  bands <- lapply(variants, function(v) {
    f <- function(t) calculate_lrmsd_i_nested_models(spm_pp, t[1], 2^t[2] - 1)[[v]]
    delta_band(f, t_hat, fit$cov, level, v)
  })
  dplyr::bind_cols(keys[c("i", "pdb_site")], bands)
}

#' Mean-centred nested-model divergence profiles with delta-method bands from an ML fit
#'
#' Mean-centred counterpart of [predict_lrmsd_i_nested_models_ml()]: the four nested
#' variants (MM, MS, MA, MSA), each centred by *its own* mean over the full model
#' support and banded with the column-centred gradient. MM's gradient is constant in
#' `(a1, a2)`, so its centred gradient is exactly zero and its band is zero-width (its
#' SPM-sampling band, added later, is not). See [predict_nlrmsd_i_msa_ml()] for the
#' support note.
#'
#' @param fit A list from [fit_lrmsd_i_msa_ml()] (carrying `a1`, `a2`, `cov`).
#' @param spm_pp Preprocessed data from [preprocess_spm()].
#' @param level Confidence-band coverage (default 0.95).
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
#' pp <- preprocess_spm(znb_spm)
#' ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)
#' head(predict_nlrmsd_i_nested_models_ml(ml, pp))
#' }
#' @export
predict_nlrmsd_i_nested_models_ml <- function(fit, spm_pp, level = 0.95) {
  validate_ml_fit(fit, "fit_lrmsd_i_msa_ml()")
  if (length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a single number in (0, 1)")
  }
  keys     <- calculate_lrmsd_i_nested_models(spm_pp, fit$a1, fit$a2)
  variants <- c("mm", "ms", "ma", "msa")
  t_hat    <- ml_t_hat(fit)

  bands <- lapply(variants, function(v) {
    src <- paste0("lrmsd_i_", v)
    f <- function(t) calculate_lrmsd_i_nested_models(spm_pp, t[1], 2^t[2] - 1)[[src]]
    delta_band_centred(f, t_hat, fit$cov, level, paste0("nlrmsd_i_", v))
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
#' @param spm_pp_mode Preprocessed data from [preprocess_spm_mode()].
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
#' pp <- preprocess_spm_mode(znb_spm)
#' ml <- fit_lrmsd_n_msa_ml(pp, znb_profile_n)
#' head(predict_lrmsd_n_nested_models_ml(ml, pp))
#' }
#' @export
predict_lrmsd_n_nested_models_ml <- function(fit, spm_pp_mode, level = 0.95) {
  validate_ml_fit(fit, "fit_lrmsd_n_msa_ml()")
  if (length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a single number in (0, 1)")
  }
  keys     <- calculate_lrmsd_n_nested_models(spm_pp_mode, fit$a1, fit$a2)
  variants <- c("lrmsd_n_mm", "lrmsd_n_ms", "lrmsd_n_ma", "lrmsd_n_msa")
  t_hat    <- ml_t_hat(fit)

  bands <- lapply(variants, function(v) {
    f <- function(t) calculate_lrmsd_n_nested_models(spm_pp_mode, t[1], 2^t[2] - 1)[[v]]
    delta_band(f, t_hat, fit$cov, level, v)
  })
  dplyr::bind_cols(keys["n"], bands)
}

#' Mean-centred nested-model divergence profiles with delta-method bands (mode form)
#'
#' Mode counterpart of [predict_nlrmsd_i_nested_models_ml()]: the four nested variants
#' of the per-mode profile, each centred by its own mean over the full model support
#' and banded with the column-centred gradient. MM's centred band is zero-width. No
#' `pdb_site` (modes are not residue-anchored).
#'
#' @param fit A list from [fit_lrmsd_n_msa_ml()] (carrying `a1`, `a2`, `cov`).
#' @param spm_pp_mode Preprocessed data from [preprocess_spm_mode()].
#' @param level Confidence-band coverage (default 0.95).
#' @return A tibble with one row per mode: `n` and a mean/lower/upper triple for each
#'   centred variant -- `nlrmsd_n_mm_{mean,lower,upper}`, `nlrmsd_n_ms_{mean,lower,upper}`,
#'   `nlrmsd_n_ma_{mean,lower,upper}`, `nlrmsd_n_msa_{mean,lower,upper}` (12 band columns).
#' @seealso [predict_lrmsd_n_nested_models_ml()] (uncentred counterpart);
#'   [predict_nlrmsd_i_nested_models_ml()] (site counterpart);
#'   [calculate_lrmsd_n_nested_models()] (the forward map).
#' @family prediction
#' @examples
#' \dontrun{
#' pp <- preprocess_spm_mode(znb_spm)
#' ml <- fit_lrmsd_n_msa_ml(pp, znb_profile_n)
#' head(predict_nlrmsd_n_nested_models_ml(ml, pp))
#' }
#' @export
predict_nlrmsd_n_nested_models_ml <- function(fit, spm_pp_mode, level = 0.95) {
  validate_ml_fit(fit, "fit_lrmsd_n_msa_ml()")
  if (length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a single number in (0, 1)")
  }
  keys     <- calculate_lrmsd_n_nested_models(spm_pp_mode, fit$a1, fit$a2)
  variants <- c("mm", "ms", "ma", "msa")
  t_hat    <- ml_t_hat(fit)

  bands <- lapply(variants, function(v) {
    src <- paste0("lrmsd_n_", v)
    f <- function(t) calculate_lrmsd_n_nested_models(spm_pp_mode, t[1], 2^t[2] - 1)[[src]]
    delta_band_centred(f, t_hat, fit$cov, level, paste0("nlrmsd_n_", v))
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
#' is linear), but the bands do not sum (SEs are not additive). `nphi_mut` (the MM
#' term) has a constant gradient, so its band is zero-width.
#'
#' @param fit A list from [fit_lrmsd_i_msa_ml()] (carrying `a1`, `a2`, `cov`).
#' @param spm_pp Preprocessed data from [preprocess_spm()].
#' @param level Confidence-band coverage (default 0.95).
#' @return A tibble with one row per residue: `i`, `pdb_site`, and a mean/lower/upper
#'   triple for each contribution -- `nphi_mut_{mean,lower,upper}`,
#'   `nphi_stab_{mean,lower,upper}`, `nphi_act_{mean,lower,upper}` (9 band columns).
#' @seealso [predict_nlrmsd_i_msa_ml()] (the profile it decomposes);
#'   [fit_lrmsd_i_msa_ml()] (produces the fit);
#'   [calculate_decomposition_i_msa()] (the forward map).
#' @family prediction
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)
#' head(predict_nlrmsd_i_msa_decomposition_ml(ml, pp))
#' }
#' @export
predict_nlrmsd_i_msa_decomposition_ml <- function(fit, spm_pp, level = 0.95) {
  validate_ml_fit(fit, "fit_lrmsd_i_msa_ml()")
  if (length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a single number in (0, 1)")
  }
  keys <- calculate_decomposition_i_msa(spm_pp, fit$a1, fit$a2)
  phis <- c("mut", "stab", "act")
  t_hat <- ml_t_hat(fit)

  bands <- lapply(phis, function(v) {
    src <- paste0("phi_", v)
    f <- function(t) calculate_decomposition_i_msa(spm_pp, t[1], 2^t[2] - 1)[[src]]
    delta_band_centred(f, t_hat, fit$cov, level, paste0("nphi_", v))
  })
  dplyr::bind_cols(keys[c("i", "pdb_site")], bands)
}

#' Decomposition of nlrmsd_n_msa with delta-method bands from an ML fit (mode form)
#'
#' Mode counterpart of [predict_nlrmsd_i_msa_decomposition_ml()]: the three
#' contributions of the mean-centred per-mode profile `nlrmsd_n_msa`, each with a
#' symmetric delta-method band. No `pdb_site` (modes are not residue-anchored). The
#' three `*_mean` columns sum to `nlrmsd_n_msa`; the bands do not sum. `nphi_mut`'s band
#' is zero-width.
#'
#' @param fit A list from [fit_lrmsd_n_msa_ml()] (carrying `a1`, `a2`, `cov`).
#' @param spm_pp_mode Preprocessed data from [preprocess_spm_mode()].
#' @param level Confidence-band coverage (default 0.95).
#' @return A tibble with one row per mode: `n` and a mean/lower/upper triple for each
#'   contribution -- `nphi_mut_{mean,lower,upper}`, `nphi_stab_{mean,lower,upper}`,
#'   `nphi_act_{mean,lower,upper}` (9 band columns).
#' @seealso [predict_nlrmsd_i_msa_decomposition_ml()] (site counterpart);
#'   [predict_nlrmsd_n_msa_ml()] (the profile it decomposes);
#'   [calculate_decomposition_n_msa()] (the forward map).
#' @family prediction
#' @examples
#' \dontrun{
#' pp <- preprocess_spm_mode(znb_spm)
#' ml <- fit_lrmsd_n_msa_ml(pp, znb_profile_n)
#' head(predict_nlrmsd_n_msa_decomposition_ml(ml, pp))
#' }
#' @export
predict_nlrmsd_n_msa_decomposition_ml <- function(fit, spm_pp_mode, level = 0.95) {
  validate_ml_fit(fit, "fit_lrmsd_n_msa_ml()")
  if (length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a single number in (0, 1)")
  }
  keys <- calculate_decomposition_n_msa(spm_pp_mode, fit$a1, fit$a2)
  phis <- c("mut", "stab", "act")
  t_hat <- ml_t_hat(fit)

  bands <- lapply(phis, function(v) {
    src <- paste0("phi_", v)
    f <- function(t) calculate_decomposition_n_msa(spm_pp_mode, t[1], 2^t[2] - 1)[[src]]
    delta_band_centred(f, t_hat, fit$cov, level, paste0("nphi_", v))
  })
  dplyr::bind_cols(keys["n"], bands)
}
