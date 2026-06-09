test_that("run_msa_bayesian_analysis returns the 7 v0.1 result elements", {
  # Integration test: this is where the real preprocess -> mcmc -> predict ->
  # decompose pipeline runs end-to-end. Structural assertion only (numeric output
  # is stochastic and expensive).
  res <- suppressMessages(
    run_msa_bayesian_analysis(znb_spm, znb_profile,
                              n_mcmc_iter = 200, n_burnin = 50)
  )
  expect_type(res, "list")
  expect_named(res, c(
    "observed_data",
    "parameter_samples", "parameter_summary",
    "prediction_samples", "prediction_summary",
    "decomposition_samples", "decomposition_summary"
  ))
  # model_comparison_* must NOT reappear in v0.1
  expect_false(any(grepl("model_comparison", names(res))))
})
