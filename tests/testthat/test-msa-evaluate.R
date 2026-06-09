test_that("calculate_dr2i_msa returns one finite, positive dr2_i per site", {
  pp <- preprocess_spm(znb_spm)
  d <- calculate_dr2i_msa(pp, a1 = 2, a2 = 5)
  expect_named(d, c("i", "dr2_i"))
  expect_equal(nrow(d), 228L)
  expect_true(all(is.finite(d$dr2_i)))
  expect_true(all(d$dr2_i > 0))
})

test_that("calculate_loglik_msa matches a frozen reference value", {
  # Frozen literal: captured once from the current implementation. NOT recomputed
  # here -- this catches a real change in the likelihood math or the pdb_site->i
  # join, not a re-derivation.
  pp <- preprocess_spm(znb_spm)
  ll <- calculate_loglik_msa(pp, znb_profile, a1 = 2, a2 = 5)
  expect_equal(ll, -184.3241923285, tolerance = 1e-6)
})
