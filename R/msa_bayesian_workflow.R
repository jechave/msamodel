#' Run the full Bayesian MSA analysis from raw mutation data to decomposition
#'
#' One call that takes a single-point-mutation scan and an observed divergence
#' profile through the whole site-form analysis: it preprocesses the scan, samples
#' the posterior of the selection strengths `(a1, a2)`, predicts the divergence
#' profiles, summarises both, and splits the divergence into its mutation,
#' stability, and activity contributions. It is a convenience wrapper over the
#' individual steps ([preprocess_spm()], [fit_lrmsd_i_msa_mcmc()],
#' [calculate_prediction_samples()] and the summary/decomposition functions).
#'
#' @param spm A single-point-mutation scan, as returned by [generate_spm_data()]
#'   (e.g. the bundled `znb_spm`).
#' @param observed_data A tibble of the observed profile, with columns `pdb_site`
#'   (PDB residue number) and `lrmsd_i_obs` (observed log structural divergence).
#'   Site-keyed elements of the result carry both `i` and `pdb_site`.
#' @param n_mcmc_iter Number of MCMC iterations to run.
#' @param n_burnin Number of initial iterations to discard before collecting samples.
#' @param a1_prior_range Length-2 `c(min, max)` uniform prior range for `a1`.
#' @param log2_a2_plus1_prior_range Length-2 `c(min, max)` uniform prior range for
#'   `log2(a2 + 1)`.
#' @param fix_a1 Optional value at which to hold `a1` fixed instead of sampling it.
#' @param fix_a2 Optional value at which to hold `a2` fixed instead of sampling it.
#' @return A named list with the inputs and every analysis stage: `observed_data`,
#'   `parameter_samples`, `parameter_summary`, `prediction_samples`,
#'   `prediction_summary`, `decomposition_samples`, and `decomposition_summary`.
#' @seealso the wrapped steps [fit_lrmsd_i_msa_mcmc()],
#'   [calculate_prediction_summary()], [calculate_decomposition_summary()].
#' @family fitting
#' @examples
#' \dontrun{
#' analysis <- run_msa_bayesian_analysis(znb_spm, znb_profile)
#' analysis$parameter_summary
#' }
#' @export
run_msa_bayesian_analysis <- function(spm,
                                     observed_data,
                                     n_mcmc_iter = 10000,
                                     n_burnin = 2000,
                                     a1_prior_range = c(0, 10),
                                     log2_a2_plus1_prior_range = c(0, 13),
                                     fix_a1 = NULL,
                                     fix_a2 = NULL) {
  # Preprocess data
  message("Preprocessing data...")
  spm_pp <- preprocess_spm(spm)

  # Run MCMC
  message("Running MCMC...")
  parameter_samples <- fit_lrmsd_i_msa_mcmc(
    spm_pp = spm_pp,
    observed_data = observed_data,
    n_iter = n_mcmc_iter,
    burn_in = n_burnin,
    fix_a1 = fix_a1,
    fix_a2 = fix_a2,
    a1_prior_range = a1_prior_range,
    log2_a2_plus1_prior_range = log2_a2_plus1_prior_range
  )

  # Calculate parameter summary
  message("Calculating parameter summary...")
  parameter_summary <- calculate_parameter_summary(parameter_samples)

  # Calculate predictions
  message("Calculating prediction samples...")
  prediction_samples <- calculate_prediction_samples(
    spm_pp = spm_pp,
    parameter_samples = parameter_samples
  )

  # Calculate prediction summary
  message("Calculating prediction summary...")
  prediction_summary <- calculate_prediction_summary(prediction_samples)

  # Calculate the phi decomposition
  message("Calculating decomposition samples...")
  decomposition_samples <- calculate_decomposition_samples(prediction_samples)

  # Calculate decomposition summary
  message("Calculating decomposition summary...")
  decomposition_summary <- calculate_decomposition_summary(decomposition_samples)

  # Return results
  message("Analysis complete.")
  return(list(
    observed_data = observed_data,
    parameter_samples = parameter_samples,
    parameter_summary = parameter_summary,
    prediction_samples = prediction_samples,
    prediction_summary = prediction_summary,
    decomposition_samples = decomposition_samples,
    decomposition_summary = decomposition_summary
  ))
}
