# The observed_data site-key contract: users supply {pdb_site, lrmsd_obs}; the
# package maps pdb_site -> internal index i.

test_that("unknown pdb_site in observed_data is an error, not a silent drop", {
  pp <- preprocess_spm(znb_spm)
  bad <- znb_profile
  bad$pdb_site[1] <- 999999L
  expect_error(
    calculate_loglik_msa(pp, bad, a1 = 2, a2 = 5),
    "not present in the model"
  )
})

test_that("partial site coverage gives a single finite loglik", {
  # znb_profile covers 225 of the model's 228 sites; the inner join keeps the
  # intersection and must still produce one finite value (no NA).
  pp <- preprocess_spm(znb_spm)
  expect_lt(nrow(znb_profile), nrow(pp$site_map))   # genuinely a subset
  ll <- calculate_loglik_msa(pp, znb_profile, a1 = 2, a2 = 5)
  expect_length(ll, 1L)
  expect_true(is.finite(ll))
})

test_that("pdb_site contract gives the same loglik as a manual i-keyed join", {
  # Independent-route check (NOT circular): map pdb_site -> i by hand and emulate
  # the original i-keyed likelihood, then compare to the contract path.
  pp <- preprocess_spm(znb_spm)
  a1 <- 2; a2 <- 5

  ll_contract <- calculate_loglik_msa(pp, znb_profile, a1, a2)

  obs_i <- dplyr::inner_join(znb_profile, pp$site_map, by = "pdb_site")
  pred <- calculate_dr2_i_msa(pp, a1, a2)
  pred$lrmsd_msa <- log(sqrt(pred$dr2_i))
  cmp <- dplyr::inner_join(obs_i[, c("i", "lrmsd_obs")],
                           pred[, c("i", "lrmsd_msa")], by = "i")
  res <- (cmp$lrmsd_obs - mean(cmp$lrmsd_obs)) -
         (cmp$lrmsd_msa - mean(cmp$lrmsd_msa))
  ll_manual <- sum(stats::dnorm(res, 0, stats::sd(res), log = TRUE))

  expect_equal(ll_contract, ll_manual)
})

test_that("site-keyed workflow outputs carry pdb_site", {
  res <- suppressMessages(
    run_msa_bayesian_analysis(znb_spm, znb_profile,
                              n_mcmc_iter = 60, n_burnin = 20)
  )
  expect_contains(names(res$prediction_summary), "pdb_site")
  expect_contains(names(res$decomposition_summary), "pdb_site")
})
