#' MSA fitting objectives
#' Pluggable criteria a fitter optimizes (currently the log-likelihood; future
#' criteria such as RMSE or robust objectives belong here)

#' Log-likelihood of the observed lrmsd_i profile given model parameters (site form)
#'
#' Profiled Gaussian log-likelihood for the site-form (`lrmsd_i`) profile. Shared
#' objective of both fitters — [fit_lrmsd_i_msa_ml()] maximises it,
#' [fit_lrmsd_i_msa_mcmc()] samples it — so it carries no method token.
#'
#' @param spm_pp Preprocessed data from [preprocess_spm()] (must include `site_map`,
#'   used to translate `pdb_site` to the internal `i`).
#' @param observed_data Tibble with columns `pdb_site` (PDB residue number) and
#'   `lrmsd_i_obs` (observed log structural divergence).
#' @param a1 Stability selection parameter
#' @param a2 Activity selection parameter
#' @return Log-likelihood value
#' @family objective
#' @export
calculate_loglik_lrmsd_i_msa <- function(spm_pp, observed_data, a1, a2) {

  # Generate model predictions (keyed by the internal response-site index i)
  dr2_i_msa <- calculate_dr2_i_msa(spm_pp, a1, a2)  %>%
    rename(dr2_i_msa = dr2_i)

  predictions <- dr2_i_msa %>%
    mutate(
      lrmsd_i_msa = log(sqrt(dr2_i_msa)),
    ) %>%
    dplyr::select(i, lrmsd_i_msa)

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

#' Log-likelihood of the observed lrmsd_n profile given model parameters (mode form)
#'
#' Mode counterpart of [calculate_loglik_lrmsd_i_msa()]. Profiled Gaussian
#' log-likelihood for the mode-form (`lrmsd_n`) profile; shared objective of the mode
#' fitter [fit_lrmsd_n_msa_ml()], so it carries no method token. Modes are not
#' structure-anchored, so the response index `n` is the model index directly — there
#' is no `site_map` / `pdb_site` translation (the only difference from the site form).
#'
#' @param spm_pp_mode Preprocessed data from [preprocess_spm_mode()] (energy_data +
#'   the `dr2_njm` matrix; no `site_map`).
#' @param observed_data Tibble with columns `n` (mode index) and `lrmsd_n_obs`
#'   (observed log structural divergence per mode).
#' @param a1 Stability selection parameter
#' @param a2 Activity selection parameter
#' @return Log-likelihood value
#' @family objective
#' @export
calculate_loglik_lrmsd_n_msa <- function(spm_pp_mode, observed_data, a1, a2) {

  # Generate model predictions (keyed by the response-mode index n)
  predictions <- calculate_dr2_n_msa(spm_pp_mode, a1, a2) %>%
    mutate(lrmsd_n_msa = log(sqrt(dr2_n))) %>%
    dplyr::select(n, lrmsd_n_msa)

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
