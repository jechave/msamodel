# MSA prediction layer — posterior prediction from a fitted model.
# Input = a fit object (class "msa_agq"); each predictor node-weights a forward
# model function (calculate_*, in R/model.R) over the stored quadrature nodes to
# return a posterior mean + credible band. @family prediction.

#' Predicted lrmsd_i profile with posterior credible bands from an AGQ fit
#'
#' Propagates the `(a1, a2)` posterior of an adaptive Gauss-Hermite quadrature fit
#' ([fit_lrmsd_i_msa_agq()]) to the predicted per-site divergence profile, returning
#' a posterior mean and a central credible interval at every residue. For each
#' quadrature node the full-model profile `lrmsd_i_msa` is evaluated with
#' [calculate_lrmsd_i_msa()]; the per-site posterior is then those node
#' values weighted by their (normalized) posterior masses `exp(log_weight)`.
#' Everything -- mean and band -- comes from the same quadrature, making no Gaussian
#' assumption (consistent with the fit's parameter summaries). Deterministic.
#'
#' The band is read off the node masses as weighted quantiles, so its tails are only
#' as fine as the node grid. To sharpen the band (especially at high `level`, e.g.
#' 0.99), refit [fit_lrmsd_i_msa_agq()] with a larger `n_nodes`; on the bundled
#' `znb_profile` the band is already stable at the default.
#'
#' @param object An `"msa_agq"` fit, the output of [fit_lrmsd_i_msa_agq()].
#' @param spm_pp Preprocessed data from [preprocess_spm()] (the same used for the fit).
#' @param level Credible-interval coverage (default 0.95).
#' @return A tibble with one row per residue: `i`, `pdb_site`, the posterior summary
#'   of the full-model profile (`lrmsd_i_msa_mean`, `lrmsd_i_msa_lower`,
#'   `lrmsd_i_msa_upper`), and their mean-centred counterparts (`nlrmsd_i_msa_mean`,
#'   `nlrmsd_i_msa_lower`, `nlrmsd_i_msa_upper`).
#' @seealso [fit_lrmsd_i_msa_agq()] (produces the posterior);
#'   [calculate_lrmsd_i_msa()] (evaluated at each node).
#' @family prediction
#' @examples
#' \dontrun{
#' pp  <- preprocess_spm(znb_spm)
#' agq <- fit_lrmsd_i_msa_agq(pp, znb_profile)
#' head(predict_lrmsd_i_msa_agq(agq, pp))
#' }
#' @export
predict_lrmsd_i_msa_agq <- function(object, spm_pp, level = 0.95) {
  if (!inherits(object, "msa_agq")) {
    stop("object must be an 'msa_agq' fit (from fit_lrmsd_i_msa_agq)")
  }
  if (length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a single number in (0, 1)")
  }
  nodes <- object$nodes
  w <- exp(nodes$log_weight)
  w <- w / sum(w)
  alpha <- (1 - level) / 2

  # lrmsd_i_msa at every node: [site x node] matrix, sites carried from node 1.
  # site_map recovers pdb_site (the forward rung returns only i); left_join keeps
  # the rung's row order (unique i), so first's rows stay aligned with P's rows.
  first <- calculate_lrmsd_i_msa(spm_pp, nodes$a1[1], nodes$a2[1]) %>%
    dplyr::left_join(spm_pp$site_map, by = "i")
  P <- vapply(seq_len(nrow(nodes)), function(k)
    calculate_lrmsd_i_msa(spm_pp, nodes$a1[k], nodes$a2[k])$lrmsd_i_msa,
    numeric(nrow(first)))

  mean_i  <- as.numeric(P %*% w)
  lower_i <- apply(P, 1L, function(row) weighted_quantile(row, w, alpha))
  upper_i <- apply(P, 1L, function(row) weighted_quantile(row, w, 1 - alpha))

  out <- tibble(
    i = first$i, pdb_site = first$pdb_site,
    lrmsd_i_msa_mean = mean_i,
    lrmsd_i_msa_lower = lower_i,
    lrmsd_i_msa_upper = upper_i
  )
  # Mean-centre by a single shift (the mean profile's mean) so the band stays a band,
  # for comparison with mean-centred observed profiles.
  shift <- mean(out$lrmsd_i_msa_mean)
  out %>%
    mutate(
      nlrmsd_i_msa_mean  = lrmsd_i_msa_mean  - shift,
      nlrmsd_i_msa_lower = lrmsd_i_msa_lower - shift,
      nlrmsd_i_msa_upper = lrmsd_i_msa_upper - shift
    )
}

# Node masses of an AGQ fit, normalized to sum 1. Shared by the banded predictors.
#' @noRd
agq_node_weights <- function(nodes) {
  w <- exp(nodes$log_weight)
  w / sum(w)
}

# Posterior mean + central credible band of one per-site quantity across AGQ nodes.
# `values` is a [row x node] matrix (one column per node); returns a tibble with
# <name>_mean/_lower/_upper columns. The band is weighted quantiles of the node masses.
#' @noRd
agq_band <- function(values, w, alpha, name) {
  cols <- list(
    as.numeric(values %*% w),
    apply(values, 1L, function(r) weighted_quantile(r, w, alpha)),
    apply(values, 1L, function(r) weighted_quantile(r, w, 1 - alpha))
  )
  names(cols) <- paste0(name, c("_mean", "_lower", "_upper"))
  tibble::as_tibble(cols)
}

#' Nested-model divergence profiles with posterior credible bands from an AGQ fit
#'
#' Propagates the `(a1, a2)` posterior of an adaptive Gauss-Hermite quadrature fit
#' ([fit_lrmsd_i_msa_agq()]) to all four nested model variants (MM, MS, MA, MSA),
#' returning a posterior mean and a central credible band for each, at every residue.
#' It is [predict_lrmsd_i_msa_agq()] widened from the single full model to the four
#' variants: at each quadrature node the four profiles are evaluated with
#' [calculate_lrmsd_i_nested_models()], then each is summarized across nodes by their
#' (normalized) posterior masses `exp(log_weight)`. All-quadrature, no Gaussian
#' assumption, deterministic.
#'
#' The four variants turn the two selection strengths on or off: MM `(0, 0)`,
#' MS `(a1, 0)`, MA `(0, a2)`, MSA `(a1, a2)`. Bands are read off the node masses as
#' weighted quantiles, so their tails are only as fine as the node grid; to sharpen
#' them, refit with a larger `n_nodes`.
#'
#' @param object An `"msa_agq"` fit, the output of [fit_lrmsd_i_msa_agq()].
#' @param spm_pp Preprocessed data from [preprocess_spm()] (the same used for the fit).
#' @param level Credible-interval coverage (default 0.95).
#' @return A tibble with one row per residue: `i`, `pdb_site`, and a mean/lower/upper
#'   triple for each of the four variants -- `lrmsd_i_mm_{mean,lower,upper}`,
#'   `lrmsd_i_ms_{mean,lower,upper}`, `lrmsd_i_ma_{mean,lower,upper}`,
#'   `lrmsd_i_msa_{mean,lower,upper}` (12 band columns).
#' @seealso [predict_lrmsd_i_msa_agq()] (the single full-model profile);
#'   [calculate_lrmsd_i_nested_models()] (evaluated at each node);
#'   [predict_decomposition_i_msa_agq()] (the phi decomposition, banded).
#' @family prediction
#' @examples
#' \dontrun{
#' pp  <- preprocess_spm(znb_spm)
#' agq <- fit_lrmsd_i_msa_agq(pp, znb_profile)
#' head(predict_lrmsd_i_nested_models_agq(agq, pp))
#' }
#' @export
predict_lrmsd_i_nested_models_agq <- function(object, spm_pp, level = 0.95) {
  if (!inherits(object, "msa_agq")) {
    stop("object must be an 'msa_agq' fit (from fit_lrmsd_i_msa_agq)")
  }
  if (length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a single number in (0, 1)")
  }
  nodes <- object$nodes
  w     <- agq_node_weights(nodes)
  alpha <- (1 - level) / 2

  # The forward fn carries i + pdb_site itself; take keys from node 1.
  first    <- calculate_lrmsd_i_nested_models(spm_pp, nodes$a1[1], nodes$a2[1])
  variants <- c("lrmsd_i_mm", "lrmsd_i_ms", "lrmsd_i_ma", "lrmsd_i_msa")

  # [site x node] matrix of each variant across nodes, then band it.
  bands <- lapply(variants, function(v) {
    P <- vapply(seq_len(nrow(nodes)), function(k)
      calculate_lrmsd_i_nested_models(spm_pp, nodes$a1[k], nodes$a2[k])[[v]],
      numeric(nrow(first)))
    agq_band(P, w, alpha, v)
  })

  dplyr::bind_cols(first[c("i", "pdb_site")], bands)
}

#' Divergence-decomposition contributions with credible bands from an AGQ fit
#'
#' Propagates the `(a1, a2)` posterior of an adaptive Gauss-Hermite quadrature fit
#' ([fit_lrmsd_i_msa_agq()]) to the site-wise divergence decomposition, returning a
#' posterior mean and a central credible band for each of the three contributions
#' (`phi_mut`, `phi_stab`, `phi_act`), at every residue. At each quadrature node the
#' contributions are evaluated with [calculate_decomposition_i_msa()], then each is
#' summarized across nodes by their (normalized) posterior masses `exp(log_weight)`.
#' All-quadrature, no Gaussian assumption, deterministic -- the decomposition-band
#' path: credible bands on the phi contributions without any sampling.
#'
#' Bands are read off the node masses as weighted quantiles, so their tails are only
#' as fine as the node grid; to sharpen them, refit with a larger `n_nodes`. Note the
#' band on a sum is not the sum of the per-contribution bands (quantiles do not add),
#' but the posterior *means* do add: the three `phi_*_mean` columns sum to the
#' full-model mean profile of [predict_lrmsd_i_msa_agq()].
#'
#' @param object An `"msa_agq"` fit, the output of [fit_lrmsd_i_msa_agq()].
#' @param spm_pp Preprocessed data from [preprocess_spm()] (the same used for the fit).
#' @param level Credible-interval coverage (default 0.95).
#' @return A tibble with one row per residue: `i`, `pdb_site`, and a mean/lower/upper
#'   triple for each contribution -- `phi_mut_{mean,lower,upper}`,
#'   `phi_stab_{mean,lower,upper}`, `phi_act_{mean,lower,upper}` (9 band columns).
#' @seealso [calculate_decomposition_i_msa()] (evaluated at each node);
#'   [predict_lrmsd_i_msa_agq()] (the full-model profile the means sum to);
#'   [predict_lrmsd_i_nested_models_agq()] (the four nested variants, banded).
#' @family prediction
#' @examples
#' \dontrun{
#' pp  <- preprocess_spm(znb_spm)
#' agq <- fit_lrmsd_i_msa_agq(pp, znb_profile)
#' head(predict_decomposition_i_msa_agq(agq, pp))
#' }
#' @export
predict_decomposition_i_msa_agq <- function(object, spm_pp, level = 0.95) {
  if (!inherits(object, "msa_agq")) {
    stop("object must be an 'msa_agq' fit (from fit_lrmsd_i_msa_agq)")
  }
  if (length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a single number in (0, 1)")
  }
  nodes <- object$nodes
  w     <- agq_node_weights(nodes)
  alpha <- (1 - level) / 2

  first <- calculate_decomposition_i_msa(spm_pp, nodes$a1[1], nodes$a2[1])
  phis  <- c("phi_mut", "phi_stab", "phi_act")

  bands <- lapply(phis, function(v) {
    P <- vapply(seq_len(nrow(nodes)), function(k)
      calculate_decomposition_i_msa(spm_pp, nodes$a1[k], nodes$a2[k])[[v]],
      numeric(nrow(first)))
    agq_band(P, w, alpha, v)
  })

  dplyr::bind_cols(first[c("i", "pdb_site")], bands)
}

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
# Independent siblings of the AGQ predictors: same forward maps (calculate_*), but the
# band comes from propagating the fit's asymptotic covariance through the forward map
# by the delta method, not from a node cloud. The fit carries `cov` on the
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

# Delta-method mean + symmetric band of one per-site quantity (the delta twin of
# agq_band). f(t) returns the per-site vector at parameter t; cov_t is the 2x2
# covariance on the t-scale. Per-site variance is the sandwich diag(J cov_t J^T),
# computed row-wise as rowSums((J %*% cov_t) * J). Returns a tibble with
# <name>_mean/_lower/_upper, matching agq_band's column shape.
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
#' Maximum-likelihood counterpart of [predict_lrmsd_i_msa_agq()]: propagates the
#' asymptotic covariance of a point ML fit ([fit_lrmsd_i_msa_ml()]) to the predicted
#' per-site divergence profile, returning the point profile and a symmetric error band
#' at every residue. The band is the delta method -- the forward map
#' [calculate_lrmsd_i_msa()] linearized about the fit's `(a1, a2)`, its gradient
#' sandwiched with the fit's covariance -- so it is deterministic and Gaussian
#' (symmetric), consistent with the fit's own standard errors. Where the AGQ predictor
#' returns a posterior credible band, this returns a frequentist confidence band.
#'
#' @param fit A list from [fit_lrmsd_i_msa_ml()] (carrying `a1`, `a2`, and the 2x2
#'   `cov` on the `(a1, log2(a2+1))` scale).
#' @param spm_pp Preprocessed data from [preprocess_spm()] (the same used for the fit).
#' @param level Confidence-band coverage (default 0.95).
#' @return A tibble with one row per model residue: `i`, `pdb_site`, the profile
#'   summary (`lrmsd_i_msa_mean`, `lrmsd_i_msa_lower`, `lrmsd_i_msa_upper`), and their
#'   mean-centred counterparts (`nlrmsd_i_msa_mean`, `nlrmsd_i_msa_lower`,
#'   `nlrmsd_i_msa_upper`).
#' @seealso [predict_lrmsd_i_msa_agq()] (the Bayesian/posterior counterpart);
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

  out <- dplyr::bind_cols(tibble(i = keys$i, pdb_site = keys$pdb_site), band)
  # Mean-centre by a single shift so the band stays a band (mirrors the AGQ predictor).
  shift <- mean(out$lrmsd_i_msa_mean)
  out %>%
    mutate(
      nlrmsd_i_msa_mean  = lrmsd_i_msa_mean  - shift,
      nlrmsd_i_msa_lower = lrmsd_i_msa_lower - shift,
      nlrmsd_i_msa_upper = lrmsd_i_msa_upper - shift
    )
}

#' Predicted lrmsd_n profile with delta-method error bands from an ML fit (mode form)
#'
#' Mode counterpart of [predict_lrmsd_i_msa_ml()]: propagates a mode ML fit
#' ([fit_lrmsd_n_msa_ml()]) to the predicted per-mode divergence profile with a
#' symmetric delta-method band. Modes are not residue-anchored, so there is no
#' `pdb_site` column and no AGQ counterpart (AGQ is site-only).
#'
#' @param fit A list from [fit_lrmsd_n_msa_ml()] (carrying `a1`, `a2`, `cov`).
#' @param spm_pp_mode Preprocessed data from [preprocess_spm_mode()].
#' @param level Confidence-band coverage (default 0.95).
#' @return A tibble with one row per mode: `n`, the profile summary
#'   (`lrmsd_n_msa_mean`, `lrmsd_n_msa_lower`, `lrmsd_n_msa_upper`), and their
#'   mean-centred counterparts (`nlrmsd_n_msa_mean`, `nlrmsd_n_msa_lower`,
#'   `nlrmsd_n_msa_upper`).
#' @seealso [fit_lrmsd_n_msa_ml()] (produces the fit); [predict_lrmsd_i_msa_ml()]
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

  out <- dplyr::bind_cols(tibble(n = keys$n), band)
  shift <- mean(out$lrmsd_n_msa_mean)
  out %>%
    mutate(
      nlrmsd_n_msa_mean  = lrmsd_n_msa_mean  - shift,
      nlrmsd_n_msa_lower = lrmsd_n_msa_lower - shift,
      nlrmsd_n_msa_upper = lrmsd_n_msa_upper - shift
    )
}

#' Nested-model divergence profiles with delta-method bands from an ML fit
#'
#' Maximum-likelihood counterpart of [predict_lrmsd_i_nested_models_agq()]: the four
#' nested variants (MM, MS, MA, MSA) of the predicted per-site profile, each with a
#' symmetric delta-method band, from a point ML fit ([fit_lrmsd_i_msa_ml()]). Each
#' variant's band linearizes [calculate_lrmsd_i_nested_models()] about the fit's
#' `(a1, a2)`; MM `(0, 0)` does not depend on the parameters, so its band is
#' zero-width (as in the AGQ version).
#'
#' @param fit A list from [fit_lrmsd_i_msa_ml()] (carrying `a1`, `a2`, `cov`).
#' @param spm_pp Preprocessed data from [preprocess_spm()].
#' @param level Confidence-band coverage (default 0.95).
#' @return A tibble with one row per residue: `i`, `pdb_site`, and a mean/lower/upper
#'   triple for each variant -- `lrmsd_i_mm_{mean,lower,upper}`,
#'   `lrmsd_i_ms_{mean,lower,upper}`, `lrmsd_i_ma_{mean,lower,upper}`,
#'   `lrmsd_i_msa_{mean,lower,upper}` (12 band columns).
#' @seealso [predict_lrmsd_i_nested_models_agq()] (the Bayesian counterpart);
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

#' Nested-model divergence profiles with delta-method bands from an ML fit (mode form)
#'
#' Mode counterpart of [predict_lrmsd_i_nested_models_ml()]: the four nested variants
#' of the predicted per-mode profile, each with a symmetric delta-method band, from a
#' mode ML fit ([fit_lrmsd_n_msa_ml()]). No `pdb_site` (modes are not residue-anchored)
#' and no AGQ counterpart.
#'
#' @param fit A list from [fit_lrmsd_n_msa_ml()] (carrying `a1`, `a2`, `cov`).
#' @param spm_pp_mode Preprocessed data from [preprocess_spm_mode()].
#' @param level Confidence-band coverage (default 0.95).
#' @return A tibble with one row per mode: `n` and a mean/lower/upper triple for each
#'   variant -- `lrmsd_n_mm_{mean,lower,upper}`, `lrmsd_n_ms_{mean,lower,upper}`,
#'   `lrmsd_n_ma_{mean,lower,upper}`, `lrmsd_n_msa_{mean,lower,upper}` (12 band columns).
#' @seealso [predict_lrmsd_i_nested_models_ml()] (site counterpart);
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

#' Per-site divergence decomposition with delta-method bands from an ML fit
#'
#' Maximum-likelihood counterpart of [predict_decomposition_i_msa_agq()]: the three
#' contributions (`phi_mut`, `phi_stab`, `phi_act`) of the predicted per-site profile,
#' each with a symmetric delta-method band, from a point ML fit
#' ([fit_lrmsd_i_msa_ml()]). Each band linearizes [calculate_decomposition_i_msa()]
#' about the fit's `(a1, a2)`; `phi_mut` (the MM term) does not depend on the
#' parameters, so its band is zero-width. The three `*_mean` columns sum to the
#' full-model mean profile, but the bands do not sum (SEs are not additive).
#'
#' @param fit A list from [fit_lrmsd_i_msa_ml()] (carrying `a1`, `a2`, `cov`).
#' @param spm_pp Preprocessed data from [preprocess_spm()].
#' @param level Confidence-band coverage (default 0.95).
#' @return A tibble with one row per residue: `i`, `pdb_site`, and a mean/lower/upper
#'   triple for each contribution -- `phi_mut_{mean,lower,upper}`,
#'   `phi_stab_{mean,lower,upper}`, `phi_act_{mean,lower,upper}` (9 band columns).
#' @seealso [predict_decomposition_i_msa_agq()] (the Bayesian counterpart);
#'   [calculate_decomposition_i_msa()] (the forward map).
#' @family prediction
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)
#' head(predict_decomposition_i_msa_ml(ml, pp))
#' }
#' @export
predict_decomposition_i_msa_ml <- function(fit, spm_pp, level = 0.95) {
  validate_ml_fit(fit, "fit_lrmsd_i_msa_ml()")
  if (length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a single number in (0, 1)")
  }
  keys <- calculate_decomposition_i_msa(spm_pp, fit$a1, fit$a2)
  phis <- c("phi_mut", "phi_stab", "phi_act")
  t_hat <- ml_t_hat(fit)

  bands <- lapply(phis, function(v) {
    f <- function(t) calculate_decomposition_i_msa(spm_pp, t[1], 2^t[2] - 1)[[v]]
    delta_band(f, t_hat, fit$cov, level, v)
  })
  dplyr::bind_cols(keys[c("i", "pdb_site")], bands)
}

#' Per-mode divergence decomposition with delta-method bands from an ML fit (mode form)
#'
#' Mode counterpart of [predict_decomposition_i_msa_ml()]: the three contributions of
#' the predicted per-mode profile, each with a symmetric delta-method band, from a mode
#' ML fit ([fit_lrmsd_n_msa_ml()]). No `pdb_site` and no AGQ counterpart.
#'
#' @param fit A list from [fit_lrmsd_n_msa_ml()] (carrying `a1`, `a2`, `cov`).
#' @param spm_pp_mode Preprocessed data from [preprocess_spm_mode()].
#' @param level Confidence-band coverage (default 0.95).
#' @return A tibble with one row per mode: `n` and a mean/lower/upper triple for each
#'   contribution -- `phi_mut_{mean,lower,upper}`, `phi_stab_{mean,lower,upper}`,
#'   `phi_act_{mean,lower,upper}` (9 band columns).
#' @seealso [predict_decomposition_i_msa_ml()] (site counterpart);
#'   [calculate_decomposition_n_msa()] (the forward map).
#' @family prediction
#' @examples
#' \dontrun{
#' pp <- preprocess_spm_mode(znb_spm)
#' ml <- fit_lrmsd_n_msa_ml(pp, znb_profile_n)
#' head(predict_decomposition_n_msa_ml(ml, pp))
#' }
#' @export
predict_decomposition_n_msa_ml <- function(fit, spm_pp_mode, level = 0.95) {
  validate_ml_fit(fit, "fit_lrmsd_n_msa_ml()")
  if (length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a single number in (0, 1)")
  }
  keys <- calculate_decomposition_n_msa(spm_pp_mode, fit$a1, fit$a2)
  phis <- c("phi_mut", "phi_stab", "phi_act")
  t_hat <- ml_t_hat(fit)

  bands <- lapply(phis, function(v) {
    f <- function(t) calculate_decomposition_n_msa(spm_pp_mode, t[1], 2^t[2] - 1)[[v]]
    delta_band(f, t_hat, fit$cov, level, v)
  })
  dplyr::bind_cols(keys["n"], bands)
}
