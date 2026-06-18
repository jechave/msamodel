#' Functions to fit the MSA model using Bayesian inference

#' Run MCMC using Metropolis-Hastings algorithm (matrix-based version)
#'
#' @param spm_pp Preprocessed data from preprocess_spm
#' @param observed_data Tibble with columns `pdb_site` and `lrmsd_obs`
#'   (`pdb_site` is mapped to the internal index by \code{calculate_loglik_msa})
#' @param n_iter Number of MCMC iterations
#' @param burn_in Number of burn-in iterations to discard
#' @param fix_a1 Optional fixed value for a1
#' @param fix_a2 Optional fixed value for a2
#' @param a1_prior_range Vector of length 2 giving [min, max] for a1 prior
#' @param log2_a2_plus1_prior_range Vector of length 2 giving [min, max] for log2(a2 + 1) prior
#' @return MCMC samples tibble with sample_id, a1, and a2 columns
#' @family fitting
#' @export
run_mcmc_msa <- function(spm_pp,
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
    mutate(nlrmsd_obs = lrmsd_obs - mean(lrmsd_obs))

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
  current_log_lik <- calculate_loglik_msa(spm_pp, normalized_observed_data,
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
      proposed_log_lik <- calculate_loglik_msa(spm_pp, normalized_observed_data,
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



#' Calculate lrmsd predictions for MSA, MS, MA, and MM models
#'
#' @param spm_pp Preprocessed data from preprocess_spm
#' @param parameter_samples MCMC samples tibble with sample_id, a1, a2 columns
#' @return Tibble with per-sample lrmsd predictions in wide format, with both the
#'   internal response-site index `i` and the structure-anchored `pdb_site`.
#' @family fitting
#' @export
calculate_prediction_samples <- function(spm_pp, parameter_samples) {
  # Generate the long format predictions
  long_predictions <- map_dfr(seq_len(nrow(parameter_samples)), function(j) {
    if (j %% 100 == 0) message("Processing sample ", j, "/", nrow(parameter_samples))
    sample_row <- parameter_samples[j, ]

    # Extract parameters
    sample_id <- sample_row$sample_id
    a1_param <- sample_row$a1
    a2_param <- sample_row$a2

    # MSA model (full model)
    msa <- calculate_dr2_i_msa(spm_pp, a1_param, a2_param) %>%
      mutate(
        lrmsd = log(sqrt(dr2_i)),
        sample_id = sample_id,
        model = "msa"
      ) %>%
      select(i, lrmsd, sample_id, model)

    # MS model (a2 = 0)
    ms <- calculate_dr2_i_msa(spm_pp, a1_param, 0) %>%
      mutate(
        lrmsd = log(sqrt(dr2_i)),
        sample_id = sample_id,
        model = "ms"
      ) %>%
      select(i, lrmsd, sample_id, model)

    # MA model (a1 = 0)
    ma <- calculate_dr2_i_msa(spm_pp, 0, a2_param) %>%
      mutate(
        lrmsd = log(sqrt(dr2_i)),
        sample_id = sample_id,
        model = "ma"
      ) %>%
      select(i, lrmsd, sample_id, model)

    # MM model (a1 = 0, a2 = 0)
    mm <- calculate_dr2_i_msa(spm_pp, 0, 0) %>%
      mutate(
        lrmsd = log(sqrt(dr2_i)),
        sample_id = sample_id,
        model = "mm"
      ) %>%
      select(i, lrmsd, sample_id, model)

    bind_rows(msa, ms, ma, mm)
  })

  # Pivot to wide format
  wide_predictions <- long_predictions %>%
    pivot_wider(
      id_cols = c(sample_id, i),
      names_from = model,
      values_from = lrmsd,
      names_glue = "lrmsd_{model}"
    )

  # Attach the structure-anchored pdb_site alongside the internal index i.
  wide_predictions <- wide_predictions %>%
    left_join(spm_pp$site_map, by = "i") %>%
    select(sample_id, i, pdb_site, everything())

  return(wide_predictions)
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
