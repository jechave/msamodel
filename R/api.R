# Leaf API: four profile/decomposition verbs assembled from the non-exported
# primitives (R/model.R forward maps + R/predict_band_helpers.R band machinery).
#
# HARD CONSTRAINT: these verbs and their helpers call ONLY non-exported primitives.
# They MUST NOT call any exported calculate_* / predict_* / gof_* / pfix_msa function;
# they reach DOWN to the @noRd forward-map primitives (R/model.R) instead.
#
# Each verb returns list(site = <tibble>, mode = <tibble>) -- both axes every call.
# `which` selects the quantity family: "lrmsd" (absolute) or "nlrmsd" (mean-centred).

# ---- internal assembly helpers ---------------------------------------------------

#' Attach the axis key (and, for site, pdb_site) to bare value vectors
#'
#' The forward-map primitives return bare per-response vectors; this wraps a NAMED list
#' of equal-length bare vectors into the keyed tibble the verbs return. The list names
#' become the value-column names verbatim -- the verb owns naming (including any `_se`
#' siblings); `key_profile` only prepends the axis key and, for the site axis, joins
#' `pdb_site` via `spm$site_map`.
#'
#' On the site axis `site_map` IS the key table -- `(site, pdb_site)`, one row per site,
#' in `dr2_ijm` column order -- so it is bound on positionally rather than joined against
#' a manufactured index. The row-count equality that makes that valid is asserted, not
#' assumed: a positional bind fails loud on a mismatch where a join would have silently
#' filled `pdb_site` with `NA`.
#'
#' @param cols A NAMED list of equal-length bare numeric vectors. Names become columns.
#' @param spm The `spm` object (for `site_map` on the site axis).
#' @param axis `"site"` (keys `site`, `pdb_site`) or `"mode"` (key `mode`).
#' @return A tibble: site -> `site, pdb_site, <cols...>`; mode -> `mode, <cols...>`.
#' @noRd
key_profile <- function(cols, spm, axis) {
  n <- length(cols[[1]])
  body <- tibble::as_tibble(cols)
  if (axis == "site") {
    if (nrow(spm$site_map) != n) {
      stop("site_map has ", nrow(spm$site_map), " rows but the profile has ", n,
           " sites; the spm object is inconsistent.")
    }
    dplyr::bind_cols(spm$site_map, body)
  } else {
    dplyr::bind_cols(tibble(mode = seq_len(n)), body)
  }
}

#' Value + standard-error column pair for one banded quantity
#'
#' Given a point value vector and its TOTAL per-element variance (parameter + SPM arms
#' already summed by the caller), returns a two-element named list `<name>` (the value)
#' and `<name>_se` (the standard error), so the pair lands as adjacent columns.
#'
#' @param value_v Numeric vector: the point profile.
#' @param var_v Numeric vector (same length): the TOTAL per-element variance.
#' @param name Character stem for the value column; the SE column is `<name>_se`.
#' @return A named list `list(<name> = value_v, <name>_se = sqrt(pmax(var_v, 0)))`.
#' @noRd
se_cols <- function(value_v, var_v, name) {
  out <- list(value_v, sqrt(pmax(var_v, 0)))   # pmax guards tiny negative round-off only
  names(out) <- c(name, paste0(name, "_se"))
  out
}

#' Per-axis dispatch table: the divergence matrix for one branch
#'
#' The single place the site/mode branch is chosen: `dr2_ijm` + `"site"` or
#' `dr2_njm` + `"mode"`. The verbs iterate over the two branches and call the axis-blind
#' primitives with the branch's matrix.
#'
#' Carries no name tag: both branches emit the SAME value-column vocabulary
#' (`lrmsd_msa`, `lrmsd_mm`, ...), so only the key column differs, and `key_profile()`
#' derives that from `axis`.
#'
#' @param spm The `spm` object.
#' @return A named list of two branches, each `list(dr2_mat =, axis =)`.
#' @noRd
axis_branches <- function(spm) {
  list(site = list(dr2_mat = spm$dr2_ijm, axis = "site"),
       mode = list(dr2_mat = spm$dr2_njm, axis = "mode"))
}

# ---- calculate_profiles ----------------------------------------------------------

#' Predicted divergence profiles at one selection strength (both axes)
#'
#' The model's forward per-response log structural-divergence profile at a single pair
#' of selection strengths `(a1, a2)`, on **both** response axes at once. `which` selects
#' the quantity: `"lrmsd"` (the absolute profile `log(sqrt(dr2))`) or `"nlrmsd"` (the
#' mean-centred profile the fit is on). Returns point values only -- no error bands
#' (there is no fit and hence no parameter covariance); for bands from a fit use
#' [predict_profiles()].
#'
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()].
#' @param a1 Stability selection strength (non-negative). `0` disables it.
#' @param a2 Activity selection strength (non-negative). `0` disables it.
#' @param which `"lrmsd"` (absolute) or `"nlrmsd"` (mean-centred). Default `"lrmsd"`.
#' @return A list with two tibbles. `$site`: `site`, `pdb_site`, and the profile column
#'   (`lrmsd_msa` or `nlrmsd_msa`). `$mode`: `mode` and the same profile column. Both
#'   branches use identical value-column names; only the key column differs.
#' @seealso [predict_profiles()] (the same profiles with error bands, from a fit);
#'   [calculate_decomposition()] (the profile split into contributions).
#' @family api
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' calculate_profiles(spm, a1 = 1, a2 = 1, which = "nlrmsd")$site
#' }
#' @export
calculate_profiles <- function(spm, a1, a2, which = c("lrmsd", "nlrmsd")) {
  which <- match.arg(which)
  fwd   <- if (which == "nlrmsd") nlrmsd_msa else lrmsd_msa
  out <- lapply(axis_branches(spm), function(b) {
    v    <- fwd(b$dr2_mat, spm$energy_data, a1, a2)
    name <- paste0(which, "_msa")
    key_profile(stats::setNames(list(v), name), spm, b$axis)
  })
  out
}

# ---- predict_profiles ------------------------------------------------------------

#' Predicted divergence profiles with error bands from an ML fit (both axes)
#'
#' [calculate_profiles()] with uncertainty: propagates a fit's `(a1, a2)` uncertainty
#' (delta method) and the SPM finite-mutation sampling error to the predicted profile,
#' on **both** response axes. `which` selects `"lrmsd"` or `"nlrmsd"`. Each value column
#' gets an adjacent `_se` sibling holding the standard error (parameter and SPM arms
#' summed under independence).
#'
#' For `which = "nlrmsd"` the profile is centred by its own mean over the full model
#' support, and the parameter arm uses the column-centred gradient (the mean is itself a
#' function of the parameters). Prediction centres over all model residues, agnostic to
#' which residues a given dataset observes; when overlaying observed data, centre it on
#' its own matched support.
#'
#' @param fit A list from [fit_lrmsd_msa_site()] (site) or [fit_lrmsd_msa_mode()] (mode),
#'   carrying `a1`, `a2`, and the 2x2 `cov` on the `(a1, log2(a2+1))` scale. One fit
#'   drives both axes.
#' @param spm The `spm` object from [generate_spm_data()] (the same one used for the fit).
#' @param which `"lrmsd"` (absolute) or `"nlrmsd"` (mean-centred). Default `"lrmsd"`.
#' @return A list with two tibbles. `$site`: `site`, `pdb_site`, the profile column
#'   (`lrmsd_msa` or `nlrmsd_msa`), and its `_se`. `$mode`: `mode`, the same profile
#'   column, and its `_se`.
#' @seealso [calculate_profiles()] (point values, no bands);
#'   [predict_decomposition()] (the profile split into contributions, with bands).
#' @family api
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' ml  <- fit_lrmsd_msa_site(spm, znb_profile$pdb_site, znb_profile$lrmsd_i_obs)
#' predict_profiles(ml, spm, which = "nlrmsd")$site
#' }
#' @export
predict_profiles <- function(fit, spm, which = c("lrmsd", "nlrmsd")) {
  which <- match.arg(which)
  validate_ml_fit(fit, "fit_lrmsd_msa_site() / fit_lrmsd_msa_mode()")
  centred <- which == "nlrmsd"
  mean_fwd <- if (centred) nlrmsd_msa else lrmsd_msa     # value column primitive
  spm_var  <- if (centred) var_spm_nlrmsd else var_spm_lrmsd
  theta_hat <- as_theta(fit$a1, fit$a2)

  out <- lapply(axis_branches(spm), function(b) {
    mean_v <- mean_fwd(b$dr2_mat, spm$energy_data, fit$a1, fit$a2)
    # parameter arm: ALWAYS the uncentred primitive; centring handled inside var_param_delta.
    f  <- function(theta) lrmsd_msa(b$dr2_mat, spm$energy_data, theta[1], 2^theta[2] - 1)
    vp <- var_param_delta(f, theta_hat, fit$cov, centred = centred)
    w  <- weights_jm(spm$energy_data, fit$a1, fit$a2)
    vs <- spm_var(b$dr2_mat, w)
    name <- paste0(which, "_msa")
    key_profile(se_cols(mean_v, vp + vs, name), spm, b$axis)
  })
  out
}

# ---- calculate_decomposition -----------------------------------------------------

#' Divergence decomposition at one selection strength (both axes)
#'
#' The nested-model profiles AND the three sequential contributions of the divergence
#' profile at a single `(a1, a2)`, on **both** response axes. `which` selects the family
#' in lockstep: `"lrmsd"` gives the uncentred nested models (`lrmsd_mm`...) and the
#' uncentred contributions (`phi_mut`, `phi_stab`, `phi_act`); `"nlrmsd"` gives the
#' mean-centred nested models (`nlrmsd_mm`...) and the centred contributions
#' (`nphi_mut`, `nphi_stab`, `nphi_act`). Point values only; for bands from a fit use
#' [predict_decomposition()].
#'
#' The three contributions sum exactly to the full-model (`msa`) profile.
#'
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()].
#' @param a1 Stability selection strength (non-negative).
#' @param a2 Activity selection strength (non-negative).
#' @param which `"lrmsd"` (absolute) or `"nlrmsd"` (mean-centred). Default `"lrmsd"`.
#' @return A list with two tibbles (`$site`, `$mode`). Each holds the axis key
#'   (`site`, `pdb_site` for site; `mode` for mode), the four nested-model columns, and
#'   the three contribution columns, all in the family selected by `which`. Both branches
#'   use identical value-column names.
#' @seealso [predict_decomposition()] (the same, with error bands from a fit);
#'   [calculate_profiles()] (the profile these contributions sum to).
#' @family api
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' calculate_decomposition(spm, a1 = 1, a2 = 1, which = "nlrmsd")$site
#' }
#' @export
calculate_decomposition <- function(spm, a1, a2, which = c("lrmsd", "nlrmsd")) {
  which <- match.arg(which)
  centred     <- which == "nlrmsd"
  nested_fwd  <- if (centred) nlrmsd_nested_models else lrmsd_nested_models
  decomp_fwd  <- if (centred) nlrmsd_msa_decomposition else lrmsd_msa_decomposition
  phi_prefix  <- if (centred) "nphi" else "phi"

  out <- lapply(axis_branches(spm), function(b) {
    nm  <- nested_fwd(b$dr2_mat, spm$energy_data, a1, a2)          # list mm/ms/ma/msa
    phi <- decomp_fwd(b$dr2_mat, spm$energy_data, a1, a2)          # list *_mut/stab/act
    nested_cols <- stats::setNames(
      nm[c("mm", "ms", "ma", "msa")],
      paste0(which, "_", c("mm", "ms", "ma", "msa")))
    # decomp_fwd already names phi_mut/... (lrmsd) or nphi_mut/... (nlrmsd) -- take as is.
    key_profile(c(nested_cols, phi), spm, b$axis)
  })
  out
}

# ---- predict_decomposition -------------------------------------------------------

#' Divergence decomposition with error bands from an ML fit (both axes)
#'
#' [calculate_decomposition()] with uncertainty: the nested-model profiles and the three
#' contributions, each with a `_se` sibling, on **both** response axes, from a fit. Only
#' `which = "nlrmsd"` (the mean-centred family the fit is on) is available; the uncentred
#' (`which = "lrmsd"`) band is not yet derived and errors.
#'
#' Band construction mirrors the per-family `predict_*_ml` functions:
#' - **Nested models** -- the parameter arm is differentiated at the fit's point for all
#'   four variants (only the selected column differs), while the SPM arm is evaluated at
#'   each variant's own `(a1, a2)` (MM `(0,0)`, MS `(a1,0)`, MA `(0,a2)`, MSA `(a1,a2)`).
#'   MM has a zero parameter band but a nonzero SPM band.
#' - **Contributions** -- the SPM arm is the differenced-contribution variance across the
#'   nested models on the same scan (retaining the between-model covariance), from
#'   `var_spm_nphi()` at the fit point; the parameter arm is per contribution.
#'
#' @param fit A list from [fit_lrmsd_msa_site()] (site) or [fit_lrmsd_msa_mode()] (mode),
#'   carrying `a1`, `a2`, `cov`.
#' @param spm The `spm` object from [generate_spm_data()] (the same one used for the fit).
#' @param which Must be `"nlrmsd"` (default). `"lrmsd"` is accepted by the signature but
#'   errors -- the uncentred component band is not yet derived.
#' @return A list with two tibbles (`$site`, `$mode`). Each holds the axis key (`site`,
#'   `pdb_site` for site; `mode` for mode), the four nested-model columns
#'   (`nlrmsd_mm`...) and the three contribution columns (`nphi_*`), each with its `_se`.
#' @seealso [calculate_decomposition()] (point values, no bands);
#'   [predict_profiles()] (the banded profile these contributions sum to).
#' @family api
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' ml  <- fit_lrmsd_msa_site(spm, znb_profile$pdb_site, znb_profile$lrmsd_i_obs)
#' predict_decomposition(ml, spm, which = "nlrmsd")$site
#' }
#' @export
predict_decomposition <- function(fit, spm, which = c("lrmsd", "nlrmsd")) {
  which <- match.arg(which)
  if (which == "lrmsd") {
    stop("predict_decomposition(which = \"lrmsd\") is to be developed: the uncentred ",
         "component band (standard error) is not yet derived. Use which = \"nlrmsd\".")
  }
  validate_ml_fit(fit, "fit_lrmsd_msa_site() / fit_lrmsd_msa_mode()")
  theta_hat <- as_theta(fit$a1, fit$a2)
  variant_a <- list(mm = c(0, 0), ms = c(fit$a1, 0),
                    ma = c(0, fit$a2), msa = c(fit$a1, fit$a2))

  out <- lapply(axis_branches(spm), function(b) {
    # --- nested models: param arm at the fit point (per-variant column); SPM at variant point
    nested_mean <- nlrmsd_nested_models(b$dr2_mat, spm$energy_data, fit$a1, fit$a2)
    nested_cols <- do.call(c, lapply(c("mm", "ms", "ma", "msa"), function(v) {
      f  <- function(theta) lrmsd_nested_models(b$dr2_mat, spm$energy_data, theta[1], 2^theta[2] - 1)[[v]]
      vp <- var_param_delta(f, theta_hat, fit$cov, centred = TRUE)
      a  <- variant_a[[v]]
      w  <- weights_jm(spm$energy_data, a[1], a[2])
      vs <- var_spm_nlrmsd(b$dr2_mat, w)
      se_cols(nested_mean[[v]], vp + vs, paste0("nlrmsd_", v))
    }))

    # --- contributions: SPM arm once at the fit point (differenced, cross term kept)
    comp_mean <- nlrmsd_msa_decomposition(b$dr2_mat, spm$energy_data, fit$a1, fit$a2)
    spm_var   <- var_spm_nphi(b$dr2_mat, spm$energy_data, fit$a1, fit$a2)
    comp_cols <- do.call(c, lapply(c("mut", "stab", "act"), function(v) {
      f  <- function(theta) lrmsd_msa_decomposition(b$dr2_mat, spm$energy_data, theta[1], 2^theta[2] - 1)[[paste0("phi_", v)]]
      vp <- var_param_delta(f, theta_hat, fit$cov, centred = TRUE)
      vs <- spm_var[[v]]
      se_cols(comp_mean[[paste0("nphi_", v)]], vp + vs, paste0("nphi_", v))
    }))

    key_profile(c(nested_cols, comp_cols), spm, b$axis)
  })
  out
}
