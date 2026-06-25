# Functions to fit the MSA model using Bayesian inference

#' Fit the selection strengths to an observed profile by MCMC (Bayesian)
#'
#' Estimates the two selection strengths `(a1, a2)` from an observed per-site
#' divergence profile by sampling their posterior with a Metropolis-Hastings
#' Markov chain, under uniform priors over the supplied ranges. This is the
#' Bayesian counterpart of the point-estimate fitter [fit_lrmsd_i_msa_ml()]; both
#' use the same likelihood, [calculate_loglik_lrmsd_i_msa()]. Either selection
#' strength can be held fixed to fit the other alone.
#'
#' @param spm_pp Preprocessed single-point-mutation data, the output of
#'   [preprocess_spm()].
#' @param observed_data A tibble of the observed profile, with columns `pdb_site`
#'   and `lrmsd_i_obs`, as for [calculate_loglik_lrmsd_i_msa()].
#' @param n_iter Number of MCMC iterations to run.
#' @param burn_in Number of initial iterations to discard before collecting samples.
#' @param fix_a1 Optional value at which to hold `a1` fixed instead of sampling it.
#' @param fix_a2 Optional value at which to hold `a2` fixed instead of sampling it.
#' @param a1_prior_range Length-2 `c(min, max)` uniform prior range for `a1`.
#' @param log2_a2_plus1_prior_range Length-2 `c(min, max)` uniform prior range for
#'   `log2(a2 + 1)` (the coordinate in which activity selection is sampled).
#' @return A tibble of posterior draws with one row per retained iteration and
#'   columns `sample_id`, `a1`, and `a2`.
#' @seealso [fit_lrmsd_i_msa_ml()] (the point-estimate counterpart);
#'   [calculate_loglik_lrmsd_i_msa()] (the shared likelihood);
#'   [run_msa_bayesian_analysis()] (the end-to-end workflow that wraps this).
#' @family fitting
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' post <- fit_lrmsd_i_msa_mcmc(pp, znb_profile, n_iter = 2000, burn_in = 500)
#' colMeans(post[c("a1", "a2")])
#' }
#' @export
fit_lrmsd_i_msa_mcmc <- function(spm_pp,
                         observed_data,
                         n_iter = 10000,
                         burn_in = 2000,
                         fix_a1 = NULL,
                         fix_a2 = NULL,
                         a1_prior_range = c(0, 10),
                         log2_a2_plus1_prior_range = c(0, 13)) {
  # Validate prior ranges
  if (length(a1_prior_range) != 2 || a1_prior_range[1] >= a1_prior_range[2]) {
    stop("a1_prior_range must be a vector of length 2 with min < max")
  }
  if (length(log2_a2_plus1_prior_range) != 2 || log2_a2_plus1_prior_range[1] >= log2_a2_plus1_prior_range[2]) {
    stop("log2_a2_plus1_prior_range must be a vector of length 2 with min < max")
  }

  # Create normalized version of observed data (without modifying input)
  normalized_observed_data <- observed_data %>%
    mutate(nlrmsd_i_obs = lrmsd_i_obs - mean(lrmsd_i_obs))

  # Initialize parameters
  a1_min <- a1_prior_range[1]
  a1_max <- a1_prior_range[2]
  log2_a2_plus1_min <- log2_a2_plus1_prior_range[1]
  log2_a2_plus1_max <- log2_a2_plus1_prior_range[2]

  current_a1 <- if (!is.null(fix_a1)) fix_a1 else runif(1, a1_min, a1_max)
  current_log2_a2_plus1 <- if (!is.null(fix_a2)) log2(fix_a2 + 1) else
    runif(1, log2_a2_plus1_min, log2_a2_plus1_max)
  current_a2 <- 2^current_log2_a2_plus1 - 1

  # Initial log posterior
  current_log_lik <- calculate_loglik_lrmsd_i_msa(spm_pp, normalized_observed_data,
                                          current_a1, current_a2)
  current_log_post <- current_log_lik +
    (if (!is.null(fix_a1)) 0 else dunif(current_a1, a1_min, a1_max, log = TRUE)) +
    (if (!is.null(fix_a2)) 0 else dunif(current_log2_a2_plus1, log2_a2_plus1_min, log2_a2_plus1_max, log = TRUE))

  # Storage for post-burn-in samples
  n_samples <- n_iter - burn_in
  samples <- tibble(
    sample_id = 1:n_samples,
    a1 = numeric(n_samples),
    a2 = numeric(n_samples)
  )
  sample_id <- 1
  accepts <- 0

  # MCMC loop
  for (i in 1:n_iter) {
    if (i %% 100 == 0) {
      rate <- if (i > 1) accepts / (i - 1) else 0
      message("Iteration ", i, "/", n_iter, " - Acceptance rate: ", round(rate, 3))
    }

    # Propose new parameters
    proposed_a1 <- if (!is.null(fix_a1)) fix_a1 else current_a1 + rnorm(1, 0, 0.3)
    proposed_log2_a2_plus1 <- if (!is.null(fix_a2)) log2(fix_a2 + 1) else
      current_log2_a2_plus1 + rnorm(1, 0, 0.3)
    proposed_a2 <- 2^proposed_log2_a2_plus1 - 1

    # Check bounds
    if ((!is.null(fix_a1) || (proposed_a1 >= a1_min && proposed_a1 <= a1_max)) &&
        (!is.null(fix_a2) || (proposed_log2_a2_plus1 >= log2_a2_plus1_min &&
                              proposed_log2_a2_plus1 <= log2_a2_plus1_max))) {
      proposed_log_lik <- calculate_loglik_lrmsd_i_msa(spm_pp, normalized_observed_data,
                                               proposed_a1, proposed_a2)
      proposed_log_post <- proposed_log_lik +
        (if (!is.null(fix_a1)) 0 else dunif(proposed_a1, a1_min, a1_max, log = TRUE)) +
        (if (!is.null(fix_a2)) 0 else dunif(proposed_log2_a2_plus1, log2_a2_plus1_min,
                                            log2_a2_plus1_max, log = TRUE))

      log_accept_ratio <- proposed_log_post - current_log_post
      if (!is.na(log_accept_ratio) && log(runif(1)) < log_accept_ratio) {
        current_a1 <- proposed_a1
        current_a2 <- proposed_a2
        current_log2_a2_plus1 <- proposed_log2_a2_plus1
        current_log_post <- proposed_log_post
        current_log_lik <- proposed_log_lik
        accepts <- accepts + 1
      }
    }

    # Store post-burn-in
    if (i > burn_in) {
      samples$a1[sample_id] <- current_a1
      samples$a2[sample_id] <- current_a2
      sample_id <- sample_id + 1
    }
  }

  return(samples)
}



#' Predicted divergence profiles for every posterior sample
#'
#' Turns a posterior sample of selection strengths into a posterior sample of
#' predicted profiles: for each draw of `(a1, a2)` it evaluates the four nested
#' model variants -- mutation-only, plus-stability, plus-activity, and the full
#' model -- with [calculate_lrmsd_i_nested_models()] and stacks the per-site
#' results. Feed the output to [calculate_prediction_summary()] for a summary, or
#' to [calculate_decomposition_samples()] for the contribution decomposition.
#'
#' @param spm_pp Preprocessed single-point-mutation data, the output of
#'   [preprocess_spm()].
#' @param parameter_samples A tibble of posterior draws with columns `sample_id`,
#'   `a1`, `a2`, as returned by [fit_lrmsd_i_msa_mcmc()].
#' @return A tibble with one row per (sample, site): `sample_id`, the site index
#'   `i`, the PDB residue number `pdb_site`, and the four predicted profiles
#'   `lrmsd_i_mm`, `lrmsd_i_ms`, `lrmsd_i_ma`, `lrmsd_i_msa`.
#' @seealso [calculate_lrmsd_i_nested_models()] (evaluated per draw);
#'   [calculate_prediction_summary()] and [calculate_decomposition_samples()]
#'   (consume the output).
#' @family fitting
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' post <- fit_lrmsd_i_msa_mcmc(pp, znb_profile, n_iter = 2000, burn_in = 500)
#' preds <- calculate_prediction_samples(pp, post)
#' head(preds)
#' }
#' @export
calculate_prediction_samples <- function(spm_pp, parameter_samples) {
  map_dfr(seq_len(nrow(parameter_samples)), function(j) {
    if (j %% 100 == 0) message("Processing sample ", j, "/", nrow(parameter_samples))
    sample_row <- parameter_samples[j, ]

    calculate_lrmsd_i_nested_models(spm_pp, sample_row$a1, sample_row$a2) %>%
      mutate(sample_id = sample_row$sample_id) %>%
      select(sample_id, i, pdb_site, everything())
  })
}




#' Summarise the posterior of the selection strengths
#'
#' Reduces the posterior draws of the selection strengths to one row per parameter
#' (`a1`, `a2`), reporting the posterior mean, standard deviation, median, and 95%
#' credible interval. A `sample_id` column, if present, is ignored.
#'
#' @param parameter_samples A tibble of posterior draws with columns `a1` and `a2`
#'   (and optionally `sample_id`), as returned by [fit_lrmsd_i_msa_mcmc()].
#' @return A tibble with one row per parameter and columns `parameter`, `mean`,
#'   `sd`, `median`, `lower`, and `upper`.
#' @seealso [fit_lrmsd_i_msa_mcmc()] (produces the draws);
#'   [calculate_prediction_summary()] (the per-site analogue).
#' @family fitting
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' post <- fit_lrmsd_i_msa_mcmc(pp, znb_profile, n_iter = 2000, burn_in = 500)
#' calculate_parameter_summary(post)
#' }
#' @export
calculate_parameter_summary <- function(parameter_samples) {
  parameter_samples %>%
    # Remove sample_id column if it exists
    select(-any_of("sample_id")) %>%
    pivot_longer(everything(), names_to = "parameter") %>%
    group_by(parameter) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      lower = quantile(value, 0.025, na.rm = TRUE),
      upper = quantile(value, 0.975, na.rm = TRUE)
    )
}

#' Summarise the predicted profiles across posterior samples
#'
#' Reduces the per-sample predicted profiles to a posterior summary, with one row
#' per site and model variant. For each variant it reports the mean, standard
#' deviation, median, and 95% credible interval over the samples, keeping both the
#' site index and the PDB residue number.
#'
#' @param prediction_samples A tibble of per-sample predictions in wide form (with
#'   `i`, `pdb_site`, and the four `lrmsd_i_*` columns), as returned by
#'   [calculate_prediction_samples()].
#' @return A long-format tibble with columns `i`, `pdb_site`, `variable` (the model
#'   variant), and the summary columns `mean`, `sd`, `median`, `lower`, and `upper`.
#' @seealso [calculate_prediction_samples()] (produces the input);
#'   [calculate_parameter_summary()] (the parameter analogue).
#' @family fitting
#' @examples
#' \dontrun{
#' analysis <- run_msa_bayesian_analysis(znb_spm, znb_profile)
#' head(analysis$prediction_summary)
#' }
#' @export
calculate_prediction_summary <- function(prediction_samples) {
  # Get all lrmsd columns
  lrmsd_cols <- names(prediction_samples)[grep("^lrmsd_", names(prediction_samples))]

  # Create a long format result
  result <- prediction_samples %>%
    # Step 1: Group by residue position (carry pdb_site through; 1:1 with i)
    group_by(i, pdb_site) %>%
    # Step 2: For each residue, calculate summary stats for each model
    summarise(across(all_of(lrmsd_cols),
                     list(mean = ~mean(., na.rm = TRUE),
                          sd = ~sd(., na.rm = TRUE),
                          median = ~median(., na.rm = TRUE),
                          lower = ~quantile(., 0.025, na.rm = TRUE),
                          upper = ~quantile(., 0.975, na.rm = TRUE))),
              .groups = "drop") %>%
    # Step 3: Convert to long format - first for the different variables
    pivot_longer(
      cols = -c(i, pdb_site),
      names_to = c("variable", ".value"),
      names_pattern = "(.+)_(.+)"
    )

  # Ensure proper column order
  result <- result %>%
    select(i, pdb_site, variable, mean, sd, median, lower, upper)

  return(result)
}
