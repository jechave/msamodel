# The observation site-key contract: users supply a (pdb_site, lrmsd_obs) VECTOR PAIR;
# the package maps pdb_site -> the internal site index via site_map.

test_that("unknown pdb_site in observed data is an error, not a silent drop", {
  pp <- znb_spm
  bad_site <- znb_profile$pdb_site
  bad_site[1] <- 999999L
  expect_error(
    msamodel:::resolve_site_obs(pp, bad_site, znb_profile$lrmsd_i_obs),
    "not present in the model"
  )
})

test_that("partial site coverage gives a single finite loglik", {
  # znb_profile covers 225 of the model's 228 sites; the inner join keeps the
  # intersection and must still produce one finite value (no NA).
  pp <- znb_spm
  expect_lt(nrow(znb_profile), nrow(pp$site_map))   # genuinely a subset
  ll <- msamodel:::loglik_lrmsd_msa(pp$dr2_ijm, pp$energy_data, msamodel:::resolve_site_obs(pp, znb_profile$pdb_site, znb_profile$lrmsd_i_obs), a1 = 2, a2 = 5)
  expect_length(ll, 1L)
  expect_true(is.finite(ll))
})

test_that("pdb_site contract gives the same loglik as a manual index-keyed join", {
  # Independent-route check (NOT circular): map pdb_site -> the internal site index by
  # hand and emulate the index-keyed likelihood, then compare to the contract path.
  pp <- znb_spm
  a1 <- 2; a2 <- 5

  ll_contract <- msamodel:::loglik_lrmsd_msa(pp$dr2_ijm, pp$energy_data, msamodel:::resolve_site_obs(pp, znb_profile$pdb_site, znb_profile$lrmsd_i_obs), a1, a2)

  obs_i <- dplyr::inner_join(znb_profile, pp$site_map, by = "pdb_site")
  dr2_i <- dr2_msa(pp$dr2_ijm, pp$energy_data, a1, a2)
  pred <- tibble::tibble(site = seq_along(dr2_i), lrmsd_msa = log(sqrt(dr2_i)))
  cmp <- dplyr::inner_join(obs_i[, c("site", "lrmsd_i_obs")],
                           pred[, c("site", "lrmsd_msa")], by = "site")
  res <- (cmp$lrmsd_i_obs - mean(cmp$lrmsd_i_obs)) -
         (cmp$lrmsd_msa - mean(cmp$lrmsd_msa))
  ll_manual <- sum(stats::dnorm(res, 0, sqrt(mean(res^2)), log = TRUE))

  expect_equal(ll_contract, ll_manual)
})
