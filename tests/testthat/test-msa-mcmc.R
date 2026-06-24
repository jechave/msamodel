test_that("fit_lrmsd_i_msa_mcmc validates prior ranges", {
  pp <- preprocess_spm(znb_spm)
  expect_error(
    fit_lrmsd_i_msa_mcmc(pp, znb_profile, n_iter = 5, burn_in = 1,
                 a1_prior_range = c(0)),
    "a1_prior_range"
  )
  expect_error(
    fit_lrmsd_i_msa_mcmc(pp, znb_profile, n_iter = 5, burn_in = 1,
                 a1_prior_range = c(5, 1)),
    "a1_prior_range"
  )
  expect_error(
    fit_lrmsd_i_msa_mcmc(pp, znb_profile, n_iter = 5, burn_in = 1,
                 log2_a2_plus1_prior_range = c(3, 3)),
    "log2_a2_plus1_prior_range"
  )
})

test_that("fixing a1 and a2 yields constant parameter columns", {
  # Deterministic: when both parameters are fixed, every sampled value must equal
  # the fixed value. Catches a break in the fix-parameter branch; no seed needed.
  pp <- preprocess_spm(znb_spm)
  samples <- fit_lrmsd_i_msa_mcmc(pp, znb_profile, n_iter = 30, burn_in = 10,
                          fix_a1 = 3, fix_a2 = 7)
  expect_true(all(samples$a1 == 3))
  expect_true(all(samples$a2 == 7))
})

test_that("fit_lrmsd_i_msa_mcmc returns a tibble of the expected shape", {
  pp <- preprocess_spm(znb_spm)
  samples <- fit_lrmsd_i_msa_mcmc(pp, znb_profile, n_iter = 50, burn_in = 10,
                          fix_a1 = 2, fix_a2 = 3)
  expect_named(samples, c("sample_id", "a1", "a2"))
  expect_equal(nrow(samples), 40L)
})
