#' Run Bayesian analysis for MSA model, including the phi decomposition
#'
#' @param spm SPM data
#' @param observed_data Observed data with columns `pdb_site` (PDB residue number)
#'   and `lrmsd_obs` (observed log structural divergence). `pdb_site` is the
#'   structure-anchored site key; the package maps it to the internal response-site
#'   index. Site-keyed elements of the returned list carry both `i` and `pdb_site`.
#' @param n_mcmc_iter Number of MCMC iterations
#' @param n_burnin Number of burn-in iterations to discard
#' @param a1_prior_range Vector of length 2 giving [min, max] for a1 prior
#' @param log2_a2_plus1_prior_range Vector of length 2 giving [min, max] for log2(a2 + 1) prior
#' @param fix_a1 Optional fixed value for a1
#' @param fix_a2 Optional fixed value for a2
#' @return List with MCMC results, including the phi decomposition: observed_data,
#'   parameter_samples, parameter_summary, prediction_samples, prediction_summary,
#'   decomposition_samples, decomposition_summary
#' @family fitting
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
  spm_energies_and_dr2mat <- preprocess_spm(spm)

  # Run MCMC
  message("Running MCMC...")
  parameter_samples <- run_mcmc_msa(
    spm_energies_and_dr2mat = spm_energies_and_dr2mat,
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
    spm_energies_and_dr2mat = spm_energies_and_dr2mat,
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
