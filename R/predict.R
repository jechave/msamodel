# The fitted-model layer: predict from an ML fit.
#
# These verbs answer "what does the model predict, and how uncertain is it?" -- both
# halves, always. The predicted values come from the forward maps in R/model.R evaluated
# at the fit's estimate; the standard errors come from R/predict_se.R. The model layer
# (calculate_*, R/model.R) answers the same question at a GIVEN (a1, a2), where there is
# no fit and hence nothing to be uncertain about.
#
# Each verb returns list(site = <tibble>, mode = <tibble>) -- both response axes every
# call. The two axes are the same math on a different divergence matrix (`dr2_ijm` for
# sites, `dr2_njm` for modes), and are written out rather than looped: there are exactly
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
#' axes, each value column followed by an adjacent `_se` sibling. `which` selects
#' `"lrmsd"` (the absolute profile) or `"nlrmsd"` (the mean-centred profile the fit is on).
#'
#' The standard error sums two independent sources: the fit's parameter uncertainty
#' (propagated by the delta method) and the SPM finite-mutation sampling error of the
#' scan itself.
#'
#' For `which = "nlrmsd"` the profile is centred by its own mean over the full model
#' support. Prediction centres over all model residues, agnostic to which residues a
#' given dataset observes; when overlaying observed data, centre it on its own matched
#' support.
#'
#' @param fit A list from [fit_lrmsd_msa_site()] (site) or [fit_lrmsd_msa_mode()] (mode),
#'   carrying `a1`, `a2`, and the 2x2 `cov` on the `(a1, log2(a2+1))` scale. One fit
#'   drives both axes.
#' @param spm The `spm` object from [generate_spm_data()] (the same one used for the fit).
#' @param which `"lrmsd"` (absolute) or `"nlrmsd"` (mean-centred). Default `"lrmsd"`.
#' @return A list with two tibbles. `$site`: `site`, `pdb_site`, the profile column
#'   (`lrmsd_msa` or `nlrmsd_msa`), and its `_se`. `$mode`: `mode`, the same profile
#'   column, and its `_se`.
#' @seealso [calculate_profiles()] (point values at a given `(a1, a2)`, no fit);
#'   [predict_decomposition()] (the profile split into contributions, with standard errors).
#' @family api
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' ml  <- fit_lrmsd_msa_site(spm, znb_profile$pdb_site, znb_profile$lrmsd_obs)
#' predict_profiles(ml, spm, which = "nlrmsd")$site
#' }
#' @export
predict_profiles <- function(fit, spm, which = c("lrmsd", "nlrmsd")) {
  which <- match.arg(which)
  validate_ml_fit(fit, "fit_lrmsd_msa_site() / fit_lrmsd_msa_mode()")
  column <- paste0(which, "_msa")

  if (which == "nlrmsd") {
    site_value <- nlrmsd_msa(spm$dr2_ijm, spm$energy_data, fit$a1, fit$a2)
    site_se    <- se_profile_nlrmsd(spm$dr2_ijm, spm$energy_data, fit)
    mode_value <- nlrmsd_msa(spm$dr2_njm, spm$energy_data, fit$a1, fit$a2)
    mode_se    <- se_profile_nlrmsd(spm$dr2_njm, spm$energy_data, fit)
  } else {
    site_value <- lrmsd_msa(spm$dr2_ijm, spm$energy_data, fit$a1, fit$a2)
    site_se    <- se_profile_lrmsd(spm$dr2_ijm, spm$energy_data, fit)
    mode_value <- lrmsd_msa(spm$dr2_njm, spm$energy_data, fit$a1, fit$a2)
    mode_se    <- se_profile_lrmsd(spm$dr2_njm, spm$energy_data, fit)
  }

  site_columns <- list(site_value, site_se)
  mode_columns <- list(mode_value, mode_se)
  names(site_columns) <- c(column, paste0(column, "_se"))
  names(mode_columns) <- c(column, paste0(column, "_se"))

  # key_profile prepends the axis key: site + pdb_site, or mode.
  list(site = key_profile(site_columns, spm, "site"),
       mode = key_profile(mode_columns, spm, "mode"))
}

# ---- predict_decomposition -------------------------------------------------------

#' Divergence decomposition with standard errors from an ML fit (both axes)
#'
#' [calculate_decomposition()] with uncertainty: the four nested-model profiles and the
#' three sequential contributions, each followed by its `_se`, on **both** response axes.
#'
#' Only the `"nlrmsd"` family is currently available, and it is the default so a bare
#' call works. `which = "lrmsd"` stops: the uncentred standard error is not yet derived.
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
#' @param spm The `spm` object from [generate_spm_data()] (the same one used for the fit).
#' @param which `"nlrmsd"` (default, the mean-centred family the fit is on) or `"lrmsd"`
#'   (accepted by the signature but not yet derived -- it stops).
#' @return A list with two tibbles (`$site`, `$mode`). Each holds the axis key (`site`,
#'   `pdb_site` for site; `mode` for mode), the four nested-model columns
#'   (`nlrmsd_mm`...) and the three contribution columns (`nphi_*`), each with its `_se`.
#' @seealso [calculate_decomposition()] (point values at a given `(a1, a2)`, no fit);
#'   [predict_profiles()] (the profile these contributions sum to).
#' @family api
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' ml  <- fit_lrmsd_msa_site(spm, znb_profile$pdb_site, znb_profile$lrmsd_obs)
#' predict_decomposition(ml, spm)$site
#' }
#' @export
predict_decomposition <- function(fit, spm, which = c("nlrmsd", "lrmsd")) {
  which <- match.arg(which)
  validate_ml_fit(fit, "fit_lrmsd_msa_site() / fit_lrmsd_msa_mode()")

  # --- site axis (responses are residues; dr2_ijm)
  if (which == "nlrmsd") {
    site_nested     <- nlrmsd_nested_models(spm$dr2_ijm, spm$energy_data, fit$a1, fit$a2)
    site_nested_se  <- se_nested_nlrmsd(spm$dr2_ijm, spm$energy_data, fit)
    site_phi        <- nlrmsd_msa_decomposition(spm$dr2_ijm, spm$energy_data, fit$a1, fit$a2)
    site_phi_se     <- se_components_nlrmsd(spm$dr2_ijm, spm$energy_data, fit)
  } else {
    site_nested     <- lrmsd_nested_models(spm$dr2_ijm, spm$energy_data, fit$a1, fit$a2)
    site_nested_se  <- se_nested_lrmsd(spm$dr2_ijm, spm$energy_data, fit)        # stops
    site_phi        <- lrmsd_msa_decomposition(spm$dr2_ijm, spm$energy_data, fit$a1, fit$a2)
    site_phi_se     <- se_components_lrmsd(spm$dr2_ijm, spm$energy_data, fit)    # stops
  }

  # --- mode axis (responses are normal modes; dr2_njm)
  if (which == "nlrmsd") {
    mode_nested     <- nlrmsd_nested_models(spm$dr2_njm, spm$energy_data, fit$a1, fit$a2)
    mode_nested_se  <- se_nested_nlrmsd(spm$dr2_njm, spm$energy_data, fit)
    mode_phi        <- nlrmsd_msa_decomposition(spm$dr2_njm, spm$energy_data, fit$a1, fit$a2)
    mode_phi_se     <- se_components_nlrmsd(spm$dr2_njm, spm$energy_data, fit)
  } else {
    mode_nested     <- lrmsd_nested_models(spm$dr2_njm, spm$energy_data, fit$a1, fit$a2)
    mode_nested_se  <- se_nested_lrmsd(spm$dr2_njm, spm$energy_data, fit)        # stops
    mode_phi        <- lrmsd_msa_decomposition(spm$dr2_njm, spm$energy_data, fit$a1, fit$a2)
    mode_phi_se     <- se_components_lrmsd(spm$dr2_njm, spm$energy_data, fit)    # stops
  }

  # Column order: each nested model then each component, value immediately followed by
  # its _se. `phi` names itself (phi_mut / nphi_mut ...); the nested list does not.
  phi <- if (which == "nlrmsd") "nphi" else "phi"

  site_columns <- list(site_nested$mm,  site_nested_se$mm,
                       site_nested$ms,  site_nested_se$ms,
                       site_nested$ma,  site_nested_se$ma,
                       site_nested$msa, site_nested_se$msa,
                       site_phi[[paste0(phi, "_mut")]],  site_phi_se$mut,
                       site_phi[[paste0(phi, "_stab")]], site_phi_se$stab,
                       site_phi[[paste0(phi, "_act")]],  site_phi_se$act)

  mode_columns <- list(mode_nested$mm,  mode_nested_se$mm,
                       mode_nested$ms,  mode_nested_se$ms,
                       mode_nested$ma,  mode_nested_se$ma,
                       mode_nested$msa, mode_nested_se$msa,
                       mode_phi[[paste0(phi, "_mut")]],  mode_phi_se$mut,
                       mode_phi[[paste0(phi, "_stab")]], mode_phi_se$stab,
                       mode_phi[[paste0(phi, "_act")]],  mode_phi_se$act)

  column_names <- c(paste0(which, "_mm"),  paste0(which, "_mm_se"),
                    paste0(which, "_ms"),  paste0(which, "_ms_se"),
                    paste0(which, "_ma"),  paste0(which, "_ma_se"),
                    paste0(which, "_msa"), paste0(which, "_msa_se"),
                    paste0(phi, "_mut"),   paste0(phi, "_mut_se"),
                    paste0(phi, "_stab"),  paste0(phi, "_stab_se"),
                    paste0(phi, "_act"),   paste0(phi, "_act_se"))
  names(site_columns) <- column_names
  names(mode_columns) <- column_names

  # key_profile prepends the axis key: site + pdb_site, or mode.
  list(site = key_profile(site_columns, spm, "site"),
       mode = key_profile(mode_columns, spm, "mode"))
}
