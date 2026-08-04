# MSA prediction layer — band-calculation helpers (all @noRd, internal).
# Split out of R/predict.R so the delta-method band machinery lives in one place:
# the SPM-sampling variance primitives and the parameter-uncertainty / band assembly
# used by the predict_*_ml functions. @family prediction.

# ---- ML delta-method prediction ---------------------------------------------------
# Predict from an _ml fit: evaluate the forward maps (calculate_*) at the fit's point
# estimate and propagate the fit's asymptotic covariance through the forward map by the
# delta method to a symmetric error band. The fit carries `cov` on the
# t = (a1, b = log2(a2+1)) scale, so we differentiate g(t) = calc(a1=t1, a2=2^t2-1)
# w.r.t. t directly and sandwich with `cov` -- no covariance transform. Bands are
# symmetric on the reported scale (g(t_hat) +/- z*se), mirroring the fit's own SEs.

#' Central-difference Jacobian of a vector-valued forward map
#'
#' Numerically differentiates `f` (which returns a per-element numeric vector, e.g. a
#' per-site profile) with respect to the 2-vector `t = (a1, log2(a2+1))`. Column `j`
#' of the returned Jacobian is the partial derivative w.r.t. `t[j]`. Uses a symmetric
#' two-sided difference, so each dimension costs two `f`-evaluations (four forward-map
#' calls total for the 2-vector).
#'
#' @param f Forward map: a function of the parameter 2-vector `t` returning a numeric
#'   vector (length `nelement`).
#' @param t Length-2 numeric parameter vector `(a1, log2(a2+1))` at which to differentiate.
#' @param h Finite-difference step (default `1e-5`).
#' @return An `[nelement x 2]` numeric Jacobian matrix; column `j` is `d f / d t[j]`.
#'   Forced to matrix shape even when `nelement == 1`.
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

#' Parameter-uncertainty variance arm (delta method)
#'
#' Per-element variance of `f(t)` induced by the fit's `t`-scale covariance `cov_t`,
#' by the delta method. Uncentred (`centred = FALSE`): `var = diag(J cov_t J^T)` with
#' `J` the raw Jacobian. Centred (`centred = TRUE`, for an nlrmsd quantity
#' `nq = q - mean_S(q)`): because `mean_S(q)` is itself a function of the parameters,
#' the gradient of `nq` is the COLUMN-CENTRED Jacobian `g_i - mean_S(g)`, so the
#' variance uses `Jc = J - colMeans(J)` -- NOT the raw `J` shifted by a constant.
#'
#' @param f Forward map: a function of the parameter 2-vector returning the per-element
#'   vector.
#' @param t_hat Length-2 point estimate `(a1, log2(a2+1))` at which to evaluate.
#' @param cov_t The 2x2 asymptotic covariance of the fit on the `t` scale.
#' @param centred Logical: `FALSE` for an uncentred (lrmsd) quantity, `TRUE` for a
#'   mean-centred (nlrmsd) quantity (uses the column-centred Jacobian).
#' @return A bare per-element numeric vector of variances, to be summed with the SPM arm.
#' @noRd
var_param_delta <- function(f, t_hat, cov_t, centred) {
  J <- grad_t(f, t_hat)                   # [nelement x 2]
  if (centred) J <- sweep(J, 2, colMeans(J))
  rowSums((J %*% cov_t) * J)              # diag(J cov_t J^T)
}

#' Point estimate on the fitting (`t`) scale from an ML fit
#'
#' Maps an `_ml` fit's natural-scale point estimate to the transformed coordinate
#' `t = (a1, log2(a2 + 1))` the delta method operates on.
#'
#' @param fit An `_ml` fit carrying natural-scale `a1` and `a2`.
#' @return A length-2 numeric vector `c(a1, log2(a2 + 1))`.
#' @noRd
ml_t_hat <- function(fit) c(fit$a1, log2(fit$a2 + 1))

# ---- SPM-sampling error: per-mutant contribution primitive -----------------------

#' Per-mutant contribution matrix for SPM-sampling error
#'
#' The finite-mutation sampling error of the lrmsd profile. Each response value
#' `lrmsd_v = 1/2 log(msd_v)`, with `msd_v = sum_k w_k dr2[k, ]` a weighted mean over
#' the finite mutant ensemble (rows `k` of `dr2_mat`). This returns the per-mutant
#' CONTRIBUTION matrix `h[k, ]` whose column-wise sum of squares is the delta-method
#' sampling variance of the response.
#'
#' @details
#' The contribution matrix is built as:
#' \preformatted{
#'   G   = dr2_mat / (2 msd)                          # log-sensitivity d[1/2 log dr2]
#'   G   = G - rowMeans(G)          (centred only)    # d(nlrmsd)/d(cells): subtract the
#'                                                    #   support mean this mutant drives
#'   Gc  = G - colSums(G * w)                         # deviate from the w-weighted mean
#'   h   = Gc * w                                     # contribution of mutant k
#' }
#' Then `Var_SPM(lrmsd_v[i]) = sum_k h[k,i]^2`, and for a difference of two models on
#' the SAME ensemble (a phi contrast) `Var = sum_k (h_A - h_B)[k,i]^2` -- the
#' shared-sample cross term is automatic because the same mutant `k` appears in both
#' `h`'s. Validated to ~1% against a naive uniform-resample bootstrap across all
#' profile and phi families (see `dev/reports/spm_band_validation.Rmd`).
#'
#' @param dr2_mat Numeric `[mutant k x response]` matrix of per-mutant divergences.
#' @param weights Length-`nmutant` numeric averaging weights that MUST sum to 1 (as
#'   produced by `weights_jm()`); `msd` is a weighted mean, so an unnormalised input is
#'   an upstream error and is rejected rather than silently renormalised.
#' @param centred Logical: `FALSE` selects lrmsd, `TRUE` selects nlrmsd (mean-centred
#'   over the response support).
#' @return A numeric `[mutant k x response]` contribution matrix `h`; column sums of
#'   squares give the per-response SPM-sampling variance.
#' @noRd
spm_hmat <- function(dr2_mat, weights, centred) {
  if (nrow(dr2_mat) != length(weights)) {
    stop("spm_hmat: nrow(dr2_mat) (", nrow(dr2_mat), ") must equal length(weights) (",
         length(weights), "); the mutant axes are misaligned.")
  }
  stopifnot("spm_hmat: weights must sum to 1 (normalise upstream via weights_jm())" =
              abs(sum(weights) - 1) < 1e-6)
  msd <- colSums(dr2_mat * weights)
  G   <- sweep(dr2_mat, 2, 2 * msd, "/")        # dr2 / (2 msd)
  if (centred) G <- G - rowMeans(G)             # support-centre per mutant (nlrmsd)
  Gc  <- sweep(G, 2, colSums(G * weights))      # deviate from w-weighted per-response mean
  Gc * weights
}

#' Axis-blind per-response SPM-sampling variance of the uncentred lrmsd profile
#'
#' `Var_SPM` of the UNCENTRED `lrmsd` profile, on either axis. Returns a numeric vector in
#' the column order of `dr2_mat`. The difference (phi) variance is NOT built here --
#' callers form it from two `spm_hmat()` matrices so the shared-sample covariance is
#' retained. Named for the exact quantity it produces (lrmsd, uncentred); the mean-centred
#' quantity is `var_spm_nlrmsd()`.
#'
#' @param dr2_mat Numeric `[mutant x response]` divergence matrix (`dr2_ijm` or `dr2_njm`).
#' @param weights Length-`nmutant` numeric weights.
#' @return A numeric vector of per-response SPM-sampling variances.
#' @noRd
var_spm_lrmsd <- function(dr2_mat, weights) {
  colSums(spm_hmat(dr2_mat, weights, centred = FALSE)^2)
}

#' Axis-blind per-response SPM-sampling variance of the centred nlrmsd profile
#'
#' `Var_SPM` of the mean-centred `nlrmsd` profile, on either axis. Same primitive as
#' `var_spm_lrmsd()`, support-centred contribution.
#'
#' @param dr2_mat Numeric `[mutant x response]` divergence matrix (`dr2_ijm` or `dr2_njm`).
#' @param weights Length-`nmutant` numeric weights.
#' @return A numeric vector of per-response SPM-sampling variances.
#' @noRd
var_spm_nlrmsd <- function(dr2_mat, weights) {
  colSums(spm_hmat(dr2_mat, weights, centred = TRUE)^2)
}

#' Axis-blind SPM-sampling variance of the three nlrmsd decomposition contrasts
#'
#' Each phi contribution is a linear contrast of nested models (MM/MS/MSA) on the SAME
#' scan. Builds the centred per-mutant contribution matrix `h` for each model via
#' `spm_hmat()`, then forms each contrast's sampling variance from the DIFFERENCED
#' contributions, so the shared-sample covariance (cross term) is retained. Needs the fit
#' parameters to build the MM/MS/MSA weight sets (via `weights_jm`), and returns three
#' variances rather than one.
#'
#' @param dr2_mat Numeric `[mutant x response]` divergence matrix (`dr2_ijm` or `dr2_njm`).
#' @param energy_data The per-mutant energy tibble (for the MM/MS/MSA weights).
#' @param a1 Stability selection strength (the fit's point estimate).
#' @param a2 Activity selection strength (the fit's point estimate).
#' @return A named list `mut`/`stab`/`act`, each a per-response numeric variance vector.
#' @noRd
var_spm_nphi <- function(dr2_mat, energy_data, a1, a2) {
  hmat_v <- function(A1, A2) spm_hmat(dr2_mat, weights_jm(energy_data, A1, A2), centred = TRUE)
  h_mm  <- hmat_v(0, 0); h_ms <- hmat_v(a1, 0); h_msa <- hmat_v(a1, a2)
  list(mut  = colSums(h_mm^2),               # nphi_mut  = nlrmsd_mm
       stab = colSums((h_ms  - h_mm)^2),     # nphi_stab = nlrmsd_ms  - nlrmsd_mm
       act  = colSums((h_msa - h_ms)^2))     # nphi_act  = nlrmsd_msa - nlrmsd_ms
}
