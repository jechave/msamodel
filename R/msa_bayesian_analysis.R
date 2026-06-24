#' Functions to fit the MSA model using Bayesian inference

#' Maximum-likelihood fit of the lrmsd_i MSA model by MCMC (Metropolis-Hastings)
#'
#' Bayesian (posterior-sampling) fitter for the site-form (`lrmsd_i`) profile.
#' Point-estimate counterpart: [fit_lrmsd_i_msa_ml()]. Shares the same objective
#' [calculate_loglik_lrmsd_i_msa()].
#'
#' @param spm_pp Preprocessed data from preprocess_spm
#' @param observed_data Tibble with columns `pdb_site` and `lrmsd_i_obs`
#'   (`pdb_site` is mapped to the internal index by \code{calculate_loglik_lrmsd_i_msa})
#' @param n_iter Number of MCMC iterations
#' @param burn_in Number of burn-in iterations to discard
#' @param fix_a1 Optional fixed value for a1
#' @param fix_a2 Optional fixed value for a2
#' @param a1_prior_range Vector of length 2 giving [min, max] for a1 prior
#' @param log2_a2_plus1_prior_range Vector of length 2 giving [min, max] for log2(a2 + 1) prior
#' @return MCMC samples tibble with sample_id, a1, and a2 columns
#' @family fitting
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



#' Calculate lrmsd predictions for the four nested models, per posterior sample
#'
#' For each MCMC sample, evaluates the four nested-model lrmsd profiles via the
#' single-source-of-truth [calculate_lrmsd_i_nested_models()] at that sample's
#' (a1, a2), and stacks the results.
#'
#' @param spm_pp Preprocessed data from preprocess_spm
#' @param parameter_samples MCMC samples tibble with sample_id, a1, a2 columns
#' @return Tibble with per-sample lrmsd predictions in wide format
#'   (`sample_id`, `i`, `pdb_site`, `lrmsd_i_mm`, `lrmsd_i_ms`, `lrmsd_i_ma`,
#'   `lrmsd_i_msa`), with both the internal response-site index `i` and the
#'   structure-anchored `pdb_site`.
#' @family fitting
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




#' Compute summary statistics for MCMC parameter samples
#'
#' @param parameter_samples MCMC samples tibble with columns a1, a2 (and possibly sample_id)
#' @return Tibble with mean, sd, median, and 95% CI for each parameter
#' @family fitting
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

#' Summarize lrmsd predictions across samples in long format
#'
#' @param prediction_samples Tibble with per-sample predictions in wide format
#'   (with both `i` and `pdb_site`, as returned by calculate_prediction_samples)
#' @return A tibble with predictions in long format (i, pdb_site, variable, mean, sd, median, lower, upper)
#' @family fitting
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
