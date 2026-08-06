# Standard errors of the predicted quantities (all @noRd, internal).
#
# The package predicts THREE objects, each in an lrmsd and an nlrmsd flavour:
#
#   object       fitted values (R/model.R)                   se (here)
#   ----------   -----------------------------------------   -----------------------
#   profile      lrmsd_msa / nlrmsd_msa                      se_profile_*
#   nested       lrmsd_nested_models / nlrmsd_nested_models  se_nested_*
#   components   lrmsd_msa_decomposition / nlrmsd_*          se_components_*
#
# Every se sums TWO independent sources of uncertainty and takes the square root:
#
#   parameter arm -- the fit's (a1, a2) uncertainty propagated through the forward map by
#                    the delta method. The fit carries `cov` on the theta = (a1,
#                    log2(a2+1)) scale, so the forward map is differentiated w.r.t. theta
#                    directly and sandwiched with `cov`. No covariance transform anywhere.
#   SPM arm       -- the finite-mutation sampling error of the scan: the profile is a
#                    weighted mean over a FINITE ensemble of mutants, so it carries
#                    sampling error even at perfectly known (a1, a2).
#
# lrmsd and nlrmsd are separate functions, never a flag, because they are different
# formulae. `spm_hmat()` below is the one exception -- see its docstring.

# ---- shared primitives ------------------------------------------------------------

#' Central-difference Jacobian of a vector-valued forward map
#'
#' Numerically differentiates `f` (which returns a per-element numeric vector, e.g. a
#' per-site profile) with respect to the 2-vector `theta = (a1, log2(a2+1))`. Column `j`
#' of the returned Jacobian is the partial derivative w.r.t. `theta[j]`. Uses a symmetric
#' two-sided difference, so each dimension costs two `f`-evaluations (four forward-map
#' calls total for the 2-vector).
#'
#' @param f Forward map: a function of the parameter 2-vector `theta` returning a numeric
#'   vector (length `nelement`).
#' @param theta Length-2 numeric parameter vector `(a1, log2(a2+1))` at which to differentiate.
#' @param h Finite-difference step (default `1e-5`).
#' @return An `[nelement x 2]` numeric Jacobian matrix; column `j` is `d f / d theta[j]`.
#'   Forced to matrix shape even when `nelement == 1`.
#' @noRd
grad_theta <- function(f, theta, h = 1e-5) {
  n <- length(f(theta))
  J <- vapply(seq_along(theta), function(j) {
    theta_p <- theta; theta_p[j] <- theta_p[j] + h
    theta_m <- theta; theta_m[j] <- theta_m[j] - h
    (f(theta_p) - f(theta_m)) / (2 * h)
  }, numeric(n))
  # vapply collapses to a bare length-length(theta) vector when n == 1; force [n x length(theta)].
  matrix(J, nrow = n)
}

#' Per-mutant contribution matrix for the SPM-sampling arm
#'
#' WHAT THIS IS: row `k` holds mutant `k`'s contribution to the sampling error of each
#' response. The column-wise sum of squares IS the per-response SPM-sampling variance:
#' `Var_SPM(v) = sum_k h[k, v]^2`.
#'
#' It is returned as a matrix rather than already reduced to that variance because a
#' DIFFERENCE of two models on the same mutant ensemble (a component contrast) must
#' difference the contributions BEFORE squaring -- `Var = sum_k (h_A - h_B)[k, ]^2` -- so
#' the shared-sample cross term is retained automatically, the same mutant `k` appearing
#' in both `h`'s. Reducing to a variance first would lose that.
#'
#' @details
#' Each response value is `lrmsd_v = 1/2 log(msd_v)` with `msd_v = sum_k w_k dr2[k, v]` a
#' weighted mean over the finite mutant ensemble. The contribution matrix is built as:
#' \preformatted{
#'   G   = dr2_mat / (2 msd)                          # log-sensitivity d[1/2 log dr2]
#'   G   = G - rowMeans(G)          (centred only)    # d(nlrmsd)/d(cells): subtract the
#'                                                    #   support mean this mutant drives
#'   Gc  = G - colSums(G * w)                         # deviate from the w-weighted mean
#'   h   = Gc * w                                     # contribution of mutant k
#' }
#' Validated to ~1% against a naive uniform-resample bootstrap across all profile and
#' component families (see `dev/reports/spm_band_validation.Rmd`, which calls this
#' function by name).
#'
#' This is the one `centred` flag kept in this file: it toggles a single line of one
#' shared construction rather than selecting between two formulae, and the validation
#' report depends on the signature.
#'
#' @param dr2_mat Numeric `[mutant k x response]` matrix of per-mutant divergences.
#' @param weights Length-`nmutant` numeric averaging weights that MUST sum to 1 (as
#'   produced by `weights_jm()`); `msd` is a weighted mean, so an unnormalised input is
#'   an upstream error and is rejected rather than silently renormalised.
#' @param centred Logical: `FALSE` selects lrmsd, `TRUE` selects nlrmsd (mean-centred
#'   over the response support).
#' @return A numeric `[mutant k x response]` contribution matrix `h`.
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

#' Combine the two variance arms into a standard error
#'
#' The arms are independent, so the total variance is their sum and the se is its root.
#'
#' FAIL LOUD: a negative total variance is asserted against, not clamped. Both arms are
#' sums of squares by construction (the parameter arm is a quadratic form in a
#' positive-definite `cov`; the SPM arm is a sum of squared contributions) and the fitter
#' already rejects a non-positive-definite Hessian (`R/fitting.R`), so a negative here
#' means a real defect upstream -- a mis-wired Jacobian, a corrupted covariance, or a sign
#' error in the component differencing. Clamping to zero would produce a clean-looking
#' zero se, and zero is a PLAUSIBLE value here (the MM nested model has a zero parameter
#' arm by construction), so the corruption would not even look anomalous.
#'
#' @param var_parameter Per-element numeric variance from the parameter arm.
#' @param var_spm Per-element numeric variance from the SPM arm.
#' @param what Character label naming the quantity, used in the error message.
#' @return A bare per-element numeric vector of standard errors.
#' @noRd
se_from_variance_arms <- function(var_parameter, var_spm, what) {
  total <- var_parameter + var_spm
  if (any(!is.finite(total)) || any(total < 0)) {
    bad <- which(!is.finite(total) | total < 0)
    stop("se_from_variance_arms: non-finite or negative total variance for ", what,
         " at element(s) ", paste(utils::head(bad, 5), collapse = ", "),
         if (length(bad) > 5) paste0(" (and ", length(bad) - 5, " more)") else "",
         ". This is an upstream defect (Jacobian, covariance, or contrast differencing), ",
         "not round-off; it is asserted rather than clamped.")
  }
  sqrt(total)
}

# ===================================================================================
#  OBJECT 1: the profile
# ===================================================================================

#' Parameter-arm variance of the uncentred lrmsd profile
#'
#' @param dr2_mat Numeric `[mutant x response]` divergence matrix (`dr2_ijm` or `dr2_njm`).
#' @param energy_data The per-mutant energy tibble.
#' @param fit An ML fit carrying `a1`, `a2`, `cov`.
#' @return A per-element numeric variance vector.
#' @noRd
var_param_profile_lrmsd <- function(dr2_mat, energy_data, fit) {
  J <- grad_theta(function(theta) lrmsd_msa(dr2_mat, energy_data, theta[1], 2^theta[2] - 1),
                  c(fit$a1, log2(fit$a2 + 1)))
  rowSums((J %*% fit$cov) * J)                            # diag(J cov J^T)
}

#' Parameter-arm variance of the mean-centred nlrmsd profile
#'
#' The forward map differentiated here is the UNCENTRED `lrmsd_msa`, deliberately. The
#' nlrmsd profile is `nq = q - mean_S(q)`, and `mean_S(q)` is ITSELF a function of the
#' parameters, so the gradient of `nq` is the COLUMN-CENTRED Jacobian `g_i - mean_S(g)`,
#' not the raw gradient of a shifted quantity. Differentiating the centred map instead
#' would double-count the centring and give a wrong se.
#'
#' @inheritParams var_param_profile_lrmsd
#' @return A per-element numeric variance vector.
#' @noRd
var_param_profile_nlrmsd <- function(dr2_mat, energy_data, fit) {
  J <- grad_theta(function(theta) lrmsd_msa(dr2_mat, energy_data, theta[1], 2^theta[2] - 1),
                  c(fit$a1, log2(fit$a2 + 1)))
  J <- sweep(J, 2, colMeans(J))                           # column-centred Jacobian
  rowSums((J %*% fit$cov) * J)                            # diag(J cov J^T)
}

#' SPM-arm variance of the uncentred lrmsd profile
#'
#' @inheritParams var_param_profile_lrmsd
#' @return A per-element numeric variance vector.
#' @noRd
var_spm_profile_lrmsd <- function(dr2_mat, energy_data, fit) {
  colSums(spm_hmat(dr2_mat, weights_jm(energy_data, fit$a1, fit$a2), centred = FALSE)^2)
}

#' SPM-arm variance of the mean-centred nlrmsd profile
#'
#' @inheritParams var_param_profile_lrmsd
#' @return A per-element numeric variance vector.
#' @noRd
var_spm_profile_nlrmsd <- function(dr2_mat, energy_data, fit) {
  colSums(spm_hmat(dr2_mat, weights_jm(energy_data, fit$a1, fit$a2), centred = TRUE)^2)
}

#' Standard error of the predicted uncentred lrmsd profile
#'
#' @inheritParams var_param_profile_lrmsd
#' @return A bare per-element numeric vector of standard errors.
#' @noRd
se_profile_lrmsd <- function(dr2_mat, energy_data, fit) {
  se_from_variance_arms(var_param_profile_lrmsd(dr2_mat, energy_data, fit),
                        var_spm_profile_lrmsd(dr2_mat, energy_data, fit),
                        "lrmsd_msa")
}

#' Standard error of the predicted mean-centred nlrmsd profile
#'
#' @inheritParams var_param_profile_lrmsd
#' @return A bare per-element numeric vector of standard errors.
#' @noRd
se_profile_nlrmsd <- function(dr2_mat, energy_data, fit) {
  se_from_variance_arms(var_param_profile_nlrmsd(dr2_mat, energy_data, fit),
                        var_spm_profile_nlrmsd(dr2_mat, energy_data, fit),
                        "nlrmsd_msa")
}

# ===================================================================================
#  OBJECT 2: the four nested models (MM, MS, MA, MSA)
# ===================================================================================
#
# The two arms are evaluated at DIFFERENT points, deliberately:
#
#   parameter arm -- at the FIT's estimate for all four variants, because that is the
#                    only point where a covariance is defined. Only the selected column
#                    of the forward map differs between variants.
#   SPM arm       -- at EACH VARIANT'S OWN (a1, a2), because each variant averages the
#                    finite scan with its own weights, and that is what its sampling
#                    error depends on.
#
# MM therefore has a ZERO parameter arm (it does not depend on the parameters at all) but
# a NONZERO SPM arm.

#' Parameter-arm variances of the four mean-centred nested models
#'
#' Column-centred Jacobian, for the reason given in `var_param_profile_nlrmsd()`.
#' All four are differentiated at the fit's estimate; only the selected variant differs.
#'
#' @inheritParams var_param_profile_lrmsd
#' @return A named list `mm`/`ms`/`ma`/`msa` of per-element variance vectors.
#' @noRd
var_param_nested_nlrmsd <- function(dr2_mat, energy_data, fit) {
  theta_hat <- c(fit$a1, log2(fit$a2 + 1))

  J_mm <- grad_theta(function(theta) lrmsd_nested_models(dr2_mat, energy_data, theta[1], 2^theta[2] - 1)$mm,  theta_hat)
  J_ms <- grad_theta(function(theta) lrmsd_nested_models(dr2_mat, energy_data, theta[1], 2^theta[2] - 1)$ms,  theta_hat)
  J_ma <- grad_theta(function(theta) lrmsd_nested_models(dr2_mat, energy_data, theta[1], 2^theta[2] - 1)$ma,  theta_hat)
  J_msa<- grad_theta(function(theta) lrmsd_nested_models(dr2_mat, energy_data, theta[1], 2^theta[2] - 1)$msa, theta_hat)

  J_mm  <- sweep(J_mm,  2, colMeans(J_mm))                # column-centred Jacobians
  J_ms  <- sweep(J_ms,  2, colMeans(J_ms))
  J_ma  <- sweep(J_ma,  2, colMeans(J_ma))
  J_msa <- sweep(J_msa, 2, colMeans(J_msa))

  list(mm  = rowSums((J_mm  %*% fit$cov) * J_mm),         # diag(J cov J^T)
       ms  = rowSums((J_ms  %*% fit$cov) * J_ms),
       ma  = rowSums((J_ma  %*% fit$cov) * J_ma),
       msa = rowSums((J_msa %*% fit$cov) * J_msa))
}

#' SPM-arm variances of the four mean-centred nested models
#'
#' Each variant's weights come from its OWN selection strengths, not the fit's -- the
#' four `(a1, a2)` pairs below ARE the MM/MS/MA/MSA definition.
#'
#' @inheritParams var_param_profile_lrmsd
#' @return A named list `mm`/`ms`/`ma`/`msa` of per-element variance vectors.
#' @noRd
var_spm_nested_nlrmsd <- function(dr2_mat, energy_data, fit) {
  list(
    mm  = colSums(spm_hmat(dr2_mat, weights_jm(energy_data, 0,      0),      centred = TRUE)^2),
    ms  = colSums(spm_hmat(dr2_mat, weights_jm(energy_data, fit$a1, 0),      centred = TRUE)^2),
    ma  = colSums(spm_hmat(dr2_mat, weights_jm(energy_data, 0,      fit$a2), centred = TRUE)^2),
    msa = colSums(spm_hmat(dr2_mat, weights_jm(energy_data, fit$a1, fit$a2), centred = TRUE)^2)
  )
}

#' Standard errors of the four mean-centred nested models
#'
#' @inheritParams var_param_profile_lrmsd
#' @return A named list `mm`/`ms`/`ma`/`msa` of per-element standard-error vectors.
#' @noRd
se_nested_nlrmsd <- function(dr2_mat, energy_data, fit) {
  var_parameter <- var_param_nested_nlrmsd(dr2_mat, energy_data, fit)
  var_spm       <- var_spm_nested_nlrmsd(dr2_mat, energy_data, fit)

  list(
    mm  = se_from_variance_arms(var_parameter$mm,  var_spm$mm,  "nlrmsd_mm"),
    ms  = se_from_variance_arms(var_parameter$ms,  var_spm$ms,  "nlrmsd_ms"),
    ma  = se_from_variance_arms(var_parameter$ma,  var_spm$ma,  "nlrmsd_ma"),
    msa = se_from_variance_arms(var_parameter$msa, var_spm$msa, "nlrmsd_msa")
  )
}

#' Standard errors of the four uncentred nested models -- TO BE DEVELOPED
#'
#' The uncentred nested-model standard error has not been derived. This function exists
#' so the gap has a named home rather than a guard at the top of the calling verb: the
#' lrmsd branch is undeveloped, not absent.
#'
#' @inheritParams var_param_profile_lrmsd
#' @return Never returns; stops.
#' @noRd
se_nested_lrmsd <- function(dr2_mat, energy_data, fit) {
  stop("se_nested_lrmsd() is to be developed: the uncentred nested-model standard ",
       "error is not yet derived. Use the nlrmsd family.")
}

# ===================================================================================
#  OBJECT 3: the three decomposition components (mut, stab, act)
# ===================================================================================
#
# The two arms have DIFFERENT SHAPES here, and the reason is mathematical:
#
#   SPM arm       -- CANNOT be done per component. Each component is a contrast of nested
#                    models on the SAME mutant scan, so its sampling variance must
#                    difference the per-mutant contributions BEFORE squaring, which means
#                    holding h_mm, h_ms and h_msa at once.
#   parameter arm -- CAN be done per component: the forward map returns the already-formed
#                    contrast, so the Jacobian is the Jacobian OF the difference and the
#                    between-model covariance falls out of the delta method. No
#                    differencing needed.

#' Parameter-arm variances of the three mean-centred components
#'
#' Column-centred Jacobian, for the reason given in `var_param_profile_nlrmsd()`. No
#' differencing here: `lrmsd_msa_decomposition()` already returns the formed contrast, so
#' each Jacobian is the Jacobian of a difference (see the block comment above).
#'
#' @inheritParams var_param_profile_lrmsd
#' @return A named list `mut`/`stab`/`act` of per-element variance vectors.
#' @noRd
var_param_components_nlrmsd <- function(dr2_mat, energy_data, fit) {
  theta_hat <- c(fit$a1, log2(fit$a2 + 1))

  J_mut  <- grad_theta(function(theta) lrmsd_msa_decomposition(dr2_mat, energy_data, theta[1], 2^theta[2] - 1)$phi_mut,  theta_hat)
  J_stab <- grad_theta(function(theta) lrmsd_msa_decomposition(dr2_mat, energy_data, theta[1], 2^theta[2] - 1)$phi_stab, theta_hat)
  J_act  <- grad_theta(function(theta) lrmsd_msa_decomposition(dr2_mat, energy_data, theta[1], 2^theta[2] - 1)$phi_act,  theta_hat)

  J_mut  <- sweep(J_mut,  2, colMeans(J_mut))             # column-centred Jacobians
  J_stab <- sweep(J_stab, 2, colMeans(J_stab))
  J_act  <- sweep(J_act,  2, colMeans(J_act))

  list(mut  = rowSums((J_mut  %*% fit$cov) * J_mut),      # diag(J cov J^T)
       stab = rowSums((J_stab %*% fit$cov) * J_stab),
       act  = rowSums((J_act  %*% fit$cov) * J_act))
}

#' SPM-arm variances of the three mean-centred components
#'
#' `mut` is a LEVEL (`nphi_mut = nlrmsd_mm`); `stab` and `act` are DIFFERENCES, and their
#' contributions are differenced before squaring so the shared-sample cross term is
#' retained. `h_ma` is never built: the sequential split uses MM/MS/MSA only, so the MA
#' variant -- reported as its own nested column -- does not enter any contrast.
#'
#' @inheritParams var_param_profile_lrmsd
#' @return A named list `mut`/`stab`/`act` of per-element variance vectors.
#' @noRd
var_spm_components_nlrmsd <- function(dr2_mat, energy_data, fit) {
  h_mm  <- spm_hmat(dr2_mat, weights_jm(energy_data, 0,      0),      centred = TRUE)
  h_ms  <- spm_hmat(dr2_mat, weights_jm(energy_data, fit$a1, 0),      centred = TRUE)
  h_msa <- spm_hmat(dr2_mat, weights_jm(energy_data, fit$a1, fit$a2), centred = TRUE)

  list(mut  = colSums(h_mm^2),               # nphi_mut  = nlrmsd_mm
       stab = colSums((h_ms  - h_mm)^2),     # nphi_stab = nlrmsd_ms  - nlrmsd_mm
       act  = colSums((h_msa - h_ms)^2))     # nphi_act  = nlrmsd_msa - nlrmsd_ms
}

#' Standard errors of the three mean-centred components
#'
#' @inheritParams var_param_profile_lrmsd
#' @return A named list `mut`/`stab`/`act` of per-element standard-error vectors.
#' @noRd
se_components_nlrmsd <- function(dr2_mat, energy_data, fit) {
  var_parameter <- var_param_components_nlrmsd(dr2_mat, energy_data, fit)
  var_spm       <- var_spm_components_nlrmsd(dr2_mat, energy_data, fit)

  list(
    mut  = se_from_variance_arms(var_parameter$mut,  var_spm$mut,  "nphi_mut"),
    stab = se_from_variance_arms(var_parameter$stab, var_spm$stab, "nphi_stab"),
    act  = se_from_variance_arms(var_parameter$act,  var_spm$act,  "nphi_act")
  )
}

#' Standard errors of the three uncentred components -- TO BE DEVELOPED
#'
#' The uncentred component standard error has not been derived. This function exists so
#' the gap has a named home rather than a guard at the top of the calling verb: the lrmsd
#' branch is undeveloped, not absent.
#'
#' @inheritParams var_param_profile_lrmsd
#' @return Never returns; stops.
#' @noRd
se_components_lrmsd <- function(dr2_mat, energy_data, fit) {
  stop("se_components_lrmsd() is to be developed: the uncentred component standard ",
       "error is not yet derived. Use the nlrmsd family.")
}
