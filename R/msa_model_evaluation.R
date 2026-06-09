#' MSA Model Evaluation
#' Functions to compute MSA model outputs (dr2_i, nlrmsd, likelihood)

#' Calculate dr2_i for a single selection parameter point using matrix operations
#'
#' @param spm_energies_and_dr2mat Output from preprocess_spm
#' @param a1 Stability selection parameter value
#' @param a2 Activity selection parameter value
#' @return Tibble with columns i (site number) and dr2_i (weighted average)
#' @family evaluation
#' @export
calculate_dr2i_msa <- function(spm_energies_and_dr2mat, a1, a2) {
  energy_data <- spm_energies_and_dr2mat$energy_data
  dr2mat <- spm_energies_and_dr2mat$dr2mat

  # Fixation probabilities
  pstab_jm <- pmin(exp(-a1 * energy_data$ddg_jm), 1)
  pact_jm <- pmin(exp(-a2 * energy_data$ddgact_jm), 1)
  pfix_jm <- pstab_jm * pact_jm
  weights_jm <- pfix_jm / sum(pfix_jm)

  # Weighted average dr2_i
  dr2_i <- colSums(dr2mat * weights_jm)

  tibble(i = as.integer(colnames(dr2mat)), dr2_i = dr2_i)
}

#' Calculate log-likelihood of observed data given model parameters (matrix approach)
#'
#' @param spm_energies_and_dr2mat Preprocessed data from preprocess_spm
#'   (must include `site_map`, used to translate `pdb_site` to the internal `i`).
#' @param observed_data Tibble with columns `pdb_site` (PDB residue number) and
#'   `lrmsd_obs` (observed log structural divergence).
#' @param a1 Stability selection parameter
#' @param a2 Activity selection parameter
#' @return Log-likelihood value
#' @family evaluation
#' @export
calculate_loglik_msa <- function(spm_energies_and_dr2mat, observed_data, a1, a2) {

  # Generate model predictions (keyed by the internal response-site index i)
  dr2i_msa <- calculate_dr2i_msa(spm_energies_and_dr2mat, a1, a2)  %>%
    rename(dr2_msa = dr2_i)

  predictions <- dr2i_msa %>%
    mutate(
      lrmsd_msa = log(sqrt(dr2_msa)),
    ) %>%
    dplyr::select(i, lrmsd_msa)

  # Translate user-supplied pdb_site to the internal index i via the model's
  # site_map. pdb_site is the structure-anchored key; i is model-internal.
  site_map <- spm_energies_and_dr2mat$site_map
  observations <- observed_data %>%
    dplyr::select(pdb_site, lrmsd_obs)

  unknown <- setdiff(observations$pdb_site, site_map$pdb_site)
  if (length(unknown) > 0) {
    stop("observed_data has pdb_site value(s) not present in the model: ",
         paste(unknown, collapse = ", "))
  }

  observations <- observations %>%
    inner_join(site_map, by = "pdb_site") %>%
    dplyr::select(i, lrmsd_obs)

  # Match predictions with observations
  comparison <- observations %>%
    inner_join(predictions, by = "i") %>%
    mutate(
      nlrmsd_obs = lrmsd_obs - mean(lrmsd_obs),
      nlrmsd_msa = lrmsd_msa - mean(lrmsd_msa)
    )

  # Calculate residuals
  residuals <- comparison$nlrmsd_obs - comparison$nlrmsd_msa

  # Estimate sigma from residuals
  sigma <- sd(residuals)

  # Calculate log-likelihood using residuals
  log_lik <- sum(dnorm(residuals, 0, sigma, log = TRUE))

  return(log_lik)
}
