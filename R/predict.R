# The fitted-model layer: predict from an ML fit.
#
# These verbs answer "what does the model predict, and how uncertain is it?" -- both
# halves, always. The predicted values come from the forward maps in R/model.R evaluated
# at the fit's estimate; the standard errors come from R/predict_se.R. The model layer
# (calculate_*, R/model.R) answers the same question at a GIVEN (a1, a2), where there is
# no fit and hence nothing to be uncertain about.
#
# Each verb returns list(site = <tibble>, mode = <tibble>) -- both response axes every
# call. The two axes are the same math on a different divergence matrix (`dr2mat_site` for
# sites, `dr2mat_mode` for modes), and are written out rather than looped: there are exactly
# two of them.
#
# Goodness of fit is NOT here: the fitter computes it at the optimum and attaches it as
# `fit$gof` (see R/fitting.R), so there is no accessor to call.

#' Validate that a fit carries the delta-method inputs
#'
#' A fit must carry the point estimate + `theta`-scale covariance the delta method
#' needs. Fails loud on a wrong-type object rather than propagating an `NA` standard error.
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

# ---- predict_profiles ------------------------------------------------------------

#' Predicted divergence profiles with standard errors from an ML fit (both axes)
#'
#' [calculate_profiles()] with uncertainty: the model's per-response log
#' structural-divergence profile evaluated at a fit's `(a1, a2)`, on **both** response
#' axes, with the value column followed by its `_se` sibling. `metric` selects
#' `"lrmsd"` (the absolute profile) or `"nlrmsd"` (the mean-centred profile the fit is on).
#'
#' The standard error sums two independent sources: the fit's parameter uncertainty
#' (propagated by the delta method) and the SPM finite-mutation sampling error of the
#' scan itself.
#'
#' For `metric = "nlrmsd"` the profile is centred by its own mean over the full model
#' support. Prediction centres over all model residues, agnostic to which residues a
#' given dataset observes; when overlaying observed data, centre it on its own matched
#' support.
#'
#' @param fit A list from [fit_lrmsd_msa_site()] (site) or [fit_lrmsd_msa_mode()] (mode),
#'   carrying `a1`, `a2`, and the 2x2 `cov` on the `(a1, log2(a2+1))` scale. One fit
#'   drives both axes.
#' @param spm The `spm` object from [generate_spm()] (the same one used for the fit).
#' @param metric `"lrmsd"` (absolute) or `"nlrmsd"` (mean-centred). Default `"lrmsd"`.
#' @return A list with two tibbles. `$site`: `site`, `pdb_site`, the profile column
#'   (`lrmsd_msa` or `nlrmsd_msa`), and its `_se`. `$mode`: `mode`, the same profile
#'   column, and its `_se`.
#' @seealso [calculate_profiles()] (point values at a given `(a1, a2)`, no fit);
#'   [predict_decomposition()] (the profile split into contributions, with standard errors).
#' @family api
#' @examples
#' if (requireNamespace("bio3d", quietly = TRUE)) {
#'   ex  <- function(f) system.file("extdata", f, package = "msamodel")
#'   wt  <- penm::set_enm(bio3d::read.pdb(ex("1d6o_A.pdb")), node = "ca",
#'                               model = "ming_wall", d_max = 10.5, frustrated = FALSE)
#'   act <- read.csv(ex("1d6o_A_active_site.csv"))
#'   spm <- generate_spm(wt, pdb_site_active = act$pdb_site, ensemble = 1L)
#'
#'   obs <- read.csv(ex("1d6o_A_lrmsd_obs_site.csv"))
#'   ml  <- fit_lrmsd_msa_site(spm, obs$pdb_site, obs$lrmsd_obs)
#'   predict_profiles(ml, spm, metric = "nlrmsd")$site
#' }
#' @export
predict_profiles <- function(fit, spm, metric = c("lrmsd", "nlrmsd")) {
  metric <- match.arg(metric)
  validate_ml_fit(fit, "fit_lrmsd_msa_site() / fit_lrmsd_msa_mode()")

  # --- site axis (responses are residues; dr2mat_site)
  if (metric == "nlrmsd") {
    site_profile <- tibble(
      nlrmsd_msa    = nlrmsd_msa(spm$dr2mat_site, spm$energy_data, fit$a1, fit$a2),
      nlrmsd_msa_se = se_profile_nlrmsd(spm$dr2mat_site, spm$energy_data, fit))
  } else if (metric == "lrmsd") {
    site_profile <- tibble(
      lrmsd_msa    = lrmsd_msa(spm$dr2mat_site, spm$energy_data, fit$a1, fit$a2),
      lrmsd_msa_se = se_profile_lrmsd(spm$dr2mat_site, spm$energy_data, fit))
  } else {
    stop(unimplemented_metric_message(metric))
  }
  site_profile <- prepend_site_key(spm$site_map, site_profile)

  # --- mode axis (responses are normal modes; dr2mat_mode)
  if (metric == "nlrmsd") {
    mode_profile <- tibble(
      nlrmsd_msa    = nlrmsd_msa(spm$dr2mat_mode, spm$energy_data, fit$a1, fit$a2),
      nlrmsd_msa_se = se_profile_nlrmsd(spm$dr2mat_mode, spm$energy_data, fit))
  } else if (metric == "lrmsd") {
    mode_profile <- tibble(
      lrmsd_msa    = lrmsd_msa(spm$dr2mat_mode, spm$energy_data, fit$a1, fit$a2),
      lrmsd_msa_se = se_profile_lrmsd(spm$dr2mat_mode, spm$energy_data, fit))
  } else {
    stop(unimplemented_metric_message(metric))
  }
  mode_profile <- prepend_mode_key(spm$mode_map, mode_profile)

  list(site = site_profile, mode = mode_profile)
}

# ---- predict_decomposition -------------------------------------------------------

#' Divergence decomposition with standard errors from an ML fit (both axes)
#'
#' [calculate_decomposition()] with uncertainty: the four nested-model profiles and the
#' three sequential contributions, then the seven matching `_se` columns, on **both**
#' response axes.
#'
#' Only `"nlrmsd"` is currently available, and it is the default so a bare call works.
#' `metric = "lrmsd"` stops: the uncentred standard error is not yet derived.
#'
#' The standard-error construction differs between the two objects, for mathematical
#' reasons set out in `R/predict_se.R`:
#' - **Nested models** -- the parameter arm is differentiated at the fit's estimate for
#'   all four variants, while the SPM arm is evaluated at each variant's own `(a1, a2)`.
#'   MM has a zero parameter arm but a nonzero SPM arm.
#' - **Components** -- the SPM arm differences the per-mutant contributions across nested
#'   models before squaring (retaining the between-model covariance); the parameter arm
#'   needs no differencing, the forward map returning the already-formed contrast.
#'
#' @param fit A list from [fit_lrmsd_msa_site()] (site) or [fit_lrmsd_msa_mode()] (mode),
#'   carrying `a1`, `a2`, `cov`.
#' @param spm The `spm` object from [generate_spm()] (the same one used for the fit).
#' @param metric `"nlrmsd"` (default, the mean-centred profile the fit is on) or `"lrmsd"`
#'   (accepted by the signature but not yet derived -- it stops).
#' @return A list with two tibbles (`$site`, `$mode`). Each holds the axis key (`site`,
#'   `pdb_site` for site; `mode` for mode), the four nested-model columns
#'   (`nlrmsd_mm`...), the three contribution columns (`nphi_*`), and then the seven
#'   corresponding `_se` columns in the same order.
#' @seealso [calculate_decomposition()] (point values at a given `(a1, a2)`, no fit);
#'   [predict_profiles()] (the profile these contributions sum to).
#' @family api
#' @examples
#' if (requireNamespace("bio3d", quietly = TRUE)) {
#'   ex  <- function(f) system.file("extdata", f, package = "msamodel")
#'   wt  <- penm::set_enm(bio3d::read.pdb(ex("1d6o_A.pdb")), node = "ca",
#'                               model = "ming_wall", d_max = 10.5, frustrated = FALSE)
#'   act <- read.csv(ex("1d6o_A_active_site.csv"))
#'   spm <- generate_spm(wt, pdb_site_active = act$pdb_site, ensemble = 1L)
#'
#'   obs <- read.csv(ex("1d6o_A_lrmsd_obs_site.csv"))
#'   ml  <- fit_lrmsd_msa_site(spm, obs$pdb_site, obs$lrmsd_obs)
#'   predict_decomposition(ml, spm)$site
#' }
#' @export
predict_decomposition <- function(fit, spm, metric = c("nlrmsd", "lrmsd")) {
  metric <- match.arg(metric)
  validate_ml_fit(fit, "fit_lrmsd_msa_site() / fit_lrmsd_msa_mode()")

  # Each block: the four nested-model columns then the three contributions, and only
  # then the seven matching `_se` columns. The source lists name their elements three
  # different ways -- `nested`/`nested_se` are mm/ms/ma/msa, `decomposition` is
  # nphi_mut/..., and `decomposition_se` is a bare mut/stab/act -- so every column is
  # named where it is built.

  # --- site axis (responses are residues; dr2mat_site)
  if (metric == "nlrmsd") {
    nested           <- nlrmsd_nested_models(spm$dr2mat_site, spm$energy_data, fit$a1, fit$a2)
    nested_se        <- se_nested_nlrmsd(spm$dr2mat_site, spm$energy_data, fit)
    decomposition    <- nlrmsd_msa_decomposition(spm$dr2mat_site, spm$energy_data, fit$a1, fit$a2)
    decomposition_se <- se_components_nlrmsd(spm$dr2mat_site, spm$energy_data, fit)
    site_decomposition <- tibble(
      nlrmsd_mm     = nested$mm,
      nlrmsd_ms     = nested$ms,
      nlrmsd_ma     = nested$ma,
      nlrmsd_msa    = nested$msa,
      nphi_mut      = decomposition$nphi_mut,
      nphi_stab     = decomposition$nphi_stab,
      nphi_act      = decomposition$nphi_act,
      nlrmsd_mm_se  = nested_se$mm,
      nlrmsd_ms_se  = nested_se$ms,
      nlrmsd_ma_se  = nested_se$ma,
      nlrmsd_msa_se = nested_se$msa,
      nphi_mut_se   = decomposition_se$mut,
      nphi_stab_se  = decomposition_se$stab,
      nphi_act_se   = decomposition_se$act)
  } else if (metric == "lrmsd") {
    nested           <- lrmsd_nested_models(spm$dr2mat_site, spm$energy_data, fit$a1, fit$a2)
    nested_se        <- se_nested_lrmsd(spm$dr2mat_site, spm$energy_data, fit)       # stops
    decomposition    <- lrmsd_msa_decomposition(spm$dr2mat_site, spm$energy_data, fit$a1, fit$a2)
    decomposition_se <- se_components_lrmsd(spm$dr2mat_site, spm$energy_data, fit)   # stops
    site_decomposition <- tibble(
      lrmsd_mm     = nested$mm,
      lrmsd_ms     = nested$ms,
      lrmsd_ma     = nested$ma,
      lrmsd_msa    = nested$msa,
      phi_mut      = decomposition$phi_mut,
      phi_stab     = decomposition$phi_stab,
      phi_act      = decomposition$phi_act,
      lrmsd_mm_se  = nested_se$mm,
      lrmsd_ms_se  = nested_se$ms,
      lrmsd_ma_se  = nested_se$ma,
      lrmsd_msa_se = nested_se$msa,
      phi_mut_se   = decomposition_se$mut,
      phi_stab_se  = decomposition_se$stab,
      phi_act_se   = decomposition_se$act)
  } else {
    stop(unimplemented_metric_message(metric))
  }
  site_decomposition <- prepend_site_key(spm$site_map, site_decomposition)

  # --- mode axis (responses are normal modes; dr2mat_mode)
  if (metric == "nlrmsd") {
    nested           <- nlrmsd_nested_models(spm$dr2mat_mode, spm$energy_data, fit$a1, fit$a2)
    nested_se        <- se_nested_nlrmsd(spm$dr2mat_mode, spm$energy_data, fit)
    decomposition    <- nlrmsd_msa_decomposition(spm$dr2mat_mode, spm$energy_data, fit$a1, fit$a2)
    decomposition_se <- se_components_nlrmsd(spm$dr2mat_mode, spm$energy_data, fit)
    mode_decomposition <- tibble(
      nlrmsd_mm     = nested$mm,
      nlrmsd_ms     = nested$ms,
      nlrmsd_ma     = nested$ma,
      nlrmsd_msa    = nested$msa,
      nphi_mut      = decomposition$nphi_mut,
      nphi_stab     = decomposition$nphi_stab,
      nphi_act      = decomposition$nphi_act,
      nlrmsd_mm_se  = nested_se$mm,
      nlrmsd_ms_se  = nested_se$ms,
      nlrmsd_ma_se  = nested_se$ma,
      nlrmsd_msa_se = nested_se$msa,
      nphi_mut_se   = decomposition_se$mut,
      nphi_stab_se  = decomposition_se$stab,
      nphi_act_se   = decomposition_se$act)
  } else if (metric == "lrmsd") {
    nested           <- lrmsd_nested_models(spm$dr2mat_mode, spm$energy_data, fit$a1, fit$a2)
    nested_se        <- se_nested_lrmsd(spm$dr2mat_mode, spm$energy_data, fit)       # stops
    decomposition    <- lrmsd_msa_decomposition(spm$dr2mat_mode, spm$energy_data, fit$a1, fit$a2)
    decomposition_se <- se_components_lrmsd(spm$dr2mat_mode, spm$energy_data, fit)   # stops
    mode_decomposition <- tibble(
      lrmsd_mm     = nested$mm,
      lrmsd_ms     = nested$ms,
      lrmsd_ma     = nested$ma,
      lrmsd_msa    = nested$msa,
      phi_mut      = decomposition$phi_mut,
      phi_stab     = decomposition$phi_stab,
      phi_act      = decomposition$phi_act,
      lrmsd_mm_se  = nested_se$mm,
      lrmsd_ms_se  = nested_se$ms,
      lrmsd_ma_se  = nested_se$ma,
      lrmsd_msa_se = nested_se$msa,
      phi_mut_se   = decomposition_se$mut,
      phi_stab_se  = decomposition_se$stab,
      phi_act_se   = decomposition_se$act)
  } else {
    stop(unimplemented_metric_message(metric))
  }
  mode_decomposition <- prepend_mode_key(spm$mode_map, mode_decomposition)

  list(site = site_decomposition, mode = mode_decomposition)
}
