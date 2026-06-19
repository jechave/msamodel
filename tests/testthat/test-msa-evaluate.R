test_that("calculate_dr2_i_msa returns one finite, positive dr2_i per site", {
  pp <- preprocess_spm(znb_spm)
  d <- calculate_dr2_i_msa(pp, a1 = 2, a2 = 5)
  expect_named(d, c("i", "dr2_i"))
  expect_equal(nrow(d), 228L)
  expect_true(all(is.finite(d$dr2_i)))
  expect_true(all(d$dr2_i > 0))
})

test_that("calculate_lrmsd_i_nested_models builds the four variants at the right (a1,a2)", {
  pp <- preprocess_spm(znb_spm)
  a1 <- 2; a2 <- 5
  nested <- calculate_lrmsd_i_nested_models(pp, a1, a2)

  expect_named(nested, c("i", "pdb_site", "lrmsd_i_mm", "lrmsd_i_ms",
                         "lrmsd_i_ma", "lrmsd_i_msa"))
  expect_equal(nrow(nested), 228L)

  # Each variant = log(sqrt(dr2_i)) of calculate_dr2_i_msa at its (a1,a2) point.
  # Independent route (recompute the forward map directly), not circular.
  expect_equal(nested$lrmsd_i_mm,  log(sqrt(calculate_dr2_i_msa(pp, 0,  0 )$dr2_i)))
  expect_equal(nested$lrmsd_i_ms,  log(sqrt(calculate_dr2_i_msa(pp, a1, 0 )$dr2_i)))
  expect_equal(nested$lrmsd_i_ma,  log(sqrt(calculate_dr2_i_msa(pp, 0,  a2)$dr2_i)))
  expect_equal(nested$lrmsd_i_msa, log(sqrt(calculate_dr2_i_msa(pp, a1, a2)$dr2_i)))
})

test_that("calculate_loglik_msa matches a frozen reference value", {
  # Frozen literal: captured once from the current implementation. NOT recomputed
  # here -- this catches a real change in the likelihood math or the pdb_site->i
  # join, not a re-derivation.
  pp <- preprocess_spm(znb_spm)
  ll <- calculate_loglik_msa(pp, znb_profile, a1 = 2, a2 = 5)
  expect_equal(ll, -184.3241923285, tolerance = 1e-6)
})
