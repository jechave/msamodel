# The observed_data site-key contract: users supply {pdb_site, lrmsd_i_obs}; the
# package maps pdb_site -> internal index i.

test_that("unknown pdb_site in observed_data is an error, not a silent drop", {
  pp <- znb_spm
  bad <- znb_profile
  bad$pdb_site[1] <- 999999L
  expect_error(
    msamodel:::calculate_loglik_lrmsd_i_msa(pp, bad, a1 = 2, a2 = 5),
    "not present in the model"
  )
})

test_that("partial site coverage gives a single finite loglik", {
  # znb_profile covers 225 of the model's 228 sites; the inner join keeps the
  # intersection and must still produce one finite value (no NA).
  pp <- znb_spm
  expect_lt(nrow(znb_profile), nrow(pp$site_map))   # genuinely a subset
  ll <- msamodel:::calculate_loglik_lrmsd_i_msa(pp, znb_profile, a1 = 2, a2 = 5)
  expect_length(ll, 1L)
  expect_true(is.finite(ll))
})

test_that("pdb_site contract gives the same loglik as a manual i-keyed join", {
  # Independent-route check (NOT circular): map pdb_site -> i by hand and emulate
  # the original i-keyed likelihood, then compare to the contract path.
  pp <- znb_spm
  a1 <- 2; a2 <- 5

  ll_contract <- msamodel:::calculate_loglik_lrmsd_i_msa(pp, znb_profile, a1, a2)

  obs_i <- dplyr::inner_join(znb_profile, pp$site_map, by = "pdb_site")
  pred <- calculate_dr2_i_msa(pp, a1, a2)
  pred$lrmsd_i_msa <- log(sqrt(pred$dr2_i))
  cmp <- dplyr::inner_join(obs_i[, c("i", "lrmsd_i_obs")],
                           pred[, c("i", "lrmsd_i_msa")], by = "i")
  res <- (cmp$lrmsd_i_obs - mean(cmp$lrmsd_i_obs)) -
         (cmp$lrmsd_i_msa - mean(cmp$lrmsd_i_msa))
  ll_manual <- sum(stats::dnorm(res, 0, sqrt(mean(res^2)), log = TRUE))

  expect_equal(ll_contract, ll_manual)
})
