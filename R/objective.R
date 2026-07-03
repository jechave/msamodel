# MSA fitting objectives
# Pluggable criteria a fitter optimizes (currently the log-likelihood; future
# criteria such as RMSE or robust objectives belong here). Keeping the objective
# decoupled from the fitter is what lets a future stochastic-likelihood tree swap
# in a different criterion without touching the fit layer.

#' Log-likelihood of an observed per-site divergence profile
#'
#' Scores how well the model at a given selection strength `(a1, a2)` reproduces
#' an observed per-site divergence profile. The model prediction and the
#' observations are each mean-centred and compared under a Gaussian noise model
#' whose scale is estimated from the residuals (profiled out), giving a single
#' log-likelihood value. Both fitters use this same score: [fit_lrmsd_i_msa_ml()]
#' maximises it (point estimate) and [fit_lrmsd_i_msa_agq()] integrates it over the
#' selection strengths by adaptive quadrature (posterior).
#'
#' @param spm_pp Preprocessed single-point-mutation data, the output of
#'   [preprocess_spm()]. Used to predict the divergence profile and to map the
#'   observed PDB residue numbers to internal site indices.
#' @param observed_data A tibble of the observed profile, with columns `pdb_site`
#'   (PDB residue number) and `lrmsd_i_obs` (observed log structural divergence).
#' @param a1 Stability selection strength (non-negative).
#' @param a2 Activity selection strength (non-negative).
#' @return A single numeric value: the log-likelihood of `observed_data` under the
#'   model at `(a1, a2)`.
#' @seealso [fit_lrmsd_i_msa_ml()] and [fit_lrmsd_i_msa_agq()], which optimise and
#'   integrate this score; [calculate_loglik_lrmsd_n_msa()] for the mode form.
#' @family objective
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' calculate_loglik_lrmsd_i_msa(pp, znb_profile, a1 = 1, a2 = 1)
#' }
#' @export
calculate_loglik_lrmsd_i_msa <- function(spm_pp, observed_data, a1, a2) {

  # Generate model predictions (keyed by the internal response-site index i)
  predictions <- calculate_lrmsd_i_msa(spm_pp, a1, a2)

  # Translate user-supplied pdb_site to the internal index i via the model's
  # site_map. pdb_site is the structure-anchored key; i is model-internal.
  site_map <- spm_pp$site_map
  observations <- observed_data %>%
    dplyr::select(pdb_site, lrmsd_i_obs)

  unknown <- setdiff(observations$pdb_site, site_map$pdb_site)
  if (length(unknown) > 0) {
    stop("observed_data has pdb_site value(s) not present in the model: ",
         paste(unknown, collapse = ", "))
  }

  observations <- observations %>%
    inner_join(site_map, by = "pdb_site") %>%
    dplyr::select(i, lrmsd_i_obs)

  # Match predictions with observations
  comparison <- observations %>%
    inner_join(predictions, by = "i") %>%
    mutate(
      nlrmsd_i_obs = lrmsd_i_obs - mean(lrmsd_i_obs),
      nlrmsd_i_msa = lrmsd_i_msa - mean(lrmsd_i_msa)
    )

  # Calculate residuals
  residuals <- comparison$nlrmsd_i_obs - comparison$nlrmsd_i_msa

  # Estimate sigma from residuals. sigma is profiled out at each (a1, a2): the
  # value below is the closed-form ML estimate for r_i ~ N(0, sigma^2), namely
  # sqrt(mean(r^2)) (divisor n). (Previously sd(residuals), divisor n-1, which is
  # not the profile MLE; the (a1, a2) argmax is unchanged but sigma-derived
  # quantities — logLik, SEs, posterior width — were slightly off.)
  sigma <- sqrt(mean(residuals^2))

  # Calculate log-likelihood using residuals
  log_lik <- sum(dnorm(residuals, 0, sigma, log = TRUE))

  return(log_lik)
}

#' Log-likelihood of an observed per-mode divergence profile
#'
#' Mode-indexed counterpart of [calculate_loglik_lrmsd_i_msa()]: it scores the
#' model at a given selection strength `(a1, a2)` against an observed divergence
#' profile, but over the structure's normal modes rather than its residues. As in
#' the site form, prediction and observations are mean-centred and compared under a
#' Gaussian noise model whose scale is profiled out. Because modes are not anchored
#' to residues, the mode index `n` is used directly, with no PDB-residue mapping.
#'
#' @param spm_pp_mode Preprocessed single-point-mutation data in mode form, the
#'   output of [preprocess_spm_mode()].
#' @param observed_data A tibble of the observed mode profile, with columns `n`
#'   (mode index) and `lrmsd_n_obs` (observed log structural divergence per mode).
#' @param a1 Stability selection strength (non-negative).
#' @param a2 Activity selection strength (non-negative).
#' @return A single numeric value: the log-likelihood of `observed_data` under the
#'   model at `(a1, a2)`.
#' @seealso [fit_lrmsd_n_msa_ml()], which maximises this score;
#'   [calculate_loglik_lrmsd_i_msa()] for the site form.
#' @family objective
#' @examples
#' \dontrun{
#' pp <- preprocess_spm_mode(znb_spm)
#' calculate_loglik_lrmsd_n_msa(pp, znb_profile_n, a1 = 1, a2 = 1)
#' }
#' @export
calculate_loglik_lrmsd_n_msa <- function(spm_pp_mode, observed_data, a1, a2) {

  # Generate model predictions (keyed by the response-mode index n)
  predictions <- calculate_lrmsd_n_msa(spm_pp_mode, a1, a2)

  # Modes are not structure-anchored: n is the model index directly (no site_map).
  observations <- observed_data %>%
    dplyr::select(n, lrmsd_n_obs)

  unknown <- setdiff(observations$n, predictions$n)
  if (length(unknown) > 0) {
    stop("observed_data has mode index(es) not present in the model: ",
         paste(unknown, collapse = ", "))
  }

  # Match predictions with observations
  comparison <- observations %>%
    inner_join(predictions, by = "n") %>%
    mutate(
      nlrmsd_n_obs = lrmsd_n_obs - mean(lrmsd_n_obs),
      nlrmsd_n_msa = lrmsd_n_msa - mean(lrmsd_n_msa)
    )

  # Calculate residuals
  residuals <- comparison$nlrmsd_n_obs - comparison$nlrmsd_n_msa

  # Profiled-out sigma: closed-form ML estimate sqrt(mean(r^2)) (divisor n), same
  # as the site form (see calculate_loglik_lrmsd_i_msa).
  sigma <- sqrt(mean(residuals^2))

  log_lik <- sum(dnorm(residuals, 0, sigma, log = TRUE))

  return(log_lik)
}
