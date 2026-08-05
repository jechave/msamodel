# Site-axis forward-map / nested-model / decomposition PROPERTIES, tested on the
# @noRd primitives the leaf verbs are built from (the old exported calculate_*_i_*
# grid wrappers were deleted; the math they wrapped lives in these primitives). Each
# check is an independent-route property with a negative control, not a re-derivation
# by the code under test. The loglik frozen-reference (a KEEP @noRd) stays here too.

test_that("dr2_msa returns one finite, positive dr2 per site", {
  pp <- znb_spm
  d <- dr2_msa(pp$dr2_ijm, pp$energy_data, a1 = 2, a2 = 5)
  expect_length(d, 228L)
  expect_true(all(is.finite(d)))
  expect_true(all(d > 0))
})

test_that("lrmsd_nested_models builds the four variants at the right (a1,a2)", {
  pp <- znb_spm
  a1 <- 2; a2 <- 5
  nested <- lrmsd_nested_models(pp$dr2_ijm, pp$energy_data, a1, a2)

  expect_named(nested, c("mm", "ms", "ma", "msa"))
  # Each variant = log(sqrt(dr2)) of dr2_msa at its own (a1,a2) point.
  # Independent route (recompute the forward map directly), not circular.
  # Negative control: the wrong (a1,a2) per variant would break these equalities.
  expect_equal(nested$mm,  log(sqrt(dr2_msa(pp$dr2_ijm, pp$energy_data, 0,  0 ))))
  expect_equal(nested$ms,  log(sqrt(dr2_msa(pp$dr2_ijm, pp$energy_data, a1, 0 ))))
  expect_equal(nested$ma,  log(sqrt(dr2_msa(pp$dr2_ijm, pp$energy_data, 0,  a2))))
  expect_equal(nested$msa, log(sqrt(dr2_msa(pp$dr2_ijm, pp$energy_data, a1, a2))))
})

test_that("nlrmsd_msa centres the uncentred profile over the full support", {
  pp <- znb_spm
  a1 <- 2; a2 <- 5
  nc <- nlrmsd_msa(pp$dr2_ijm, pp$energy_data, a1, a2)

  # The centred value is the uncentred lrmsd minus its own mean over ALL model
  # residues. Negative control: an uncentred return (nlrmsd == lrmsd) fails, since
  # the profile mean is far from zero.
  lr <- lrmsd_msa(pp$dr2_ijm, pp$energy_data, a1, a2)
  expect_equal(nc, lr - mean(lr))
  expect_equal(mean(nc), 0)          # centred: mean is (machine) zero
})

test_that("nlrmsd_msa agrees with predict_profiles's point profile", {
  # The forward-map primitive and the public predictor must return the SAME centred
  # profile. Negative control: centring over the wrong support (or not at all) would
  # make the predictor -- which centres over the full model support -- disagree.
  pp <- znb_spm
  ml <- fit_lrmsd_i_msa_ml(pp, znb_profile$pdb_site, znb_profile$lrmsd_i_obs)
  fwd  <- nlrmsd_msa(pp$dr2_ijm, pp$energy_data, ml$a1, ml$a2)
  pred <- predict_profiles(ml, pp, which = "nlrmsd")$site
  expect_equal(fwd, pred$nlrmsd_msa)
})

test_that("nlrmsd_msa_decomposition contributions sum to the centred profile", {
  pp <- znb_spm
  a1 <- 2; a2 <- 5
  d <- nlrmsd_msa_decomposition(pp$dr2_ijm, pp$energy_data, a1, a2)

  expect_named(d, c("nphi_mut", "nphi_stab", "nphi_act"))
  # The three centred contributions reconstruct the centred full-model profile.
  # Independent route: nlrmsd_msa from its own primitive, not from this decomposition.
  # Negative control: dropping the centring on any nphi term breaks the equality.
  prof <- nlrmsd_msa(pp$dr2_ijm, pp$energy_data, a1, a2)
  expect_equal(d$nphi_mut + d$nphi_stab + d$nphi_act, prof)
})

test_that("nlrmsd_nested_models centres each variant independently", {
  pp <- znb_spm
  a1 <- 2; a2 <- 5
  nc <- nlrmsd_nested_models(pp$dr2_ijm, pp$energy_data, a1, a2)

  expect_named(nc, c("mm", "ms", "ma", "msa"))
  # Each variant centred by ITS OWN mean (independent route via the uncentred twin).
  # Negative control: centring every column by a single shared mean would fail.
  un <- lrmsd_nested_models(pp$dr2_ijm, pp$energy_data, a1, a2)
  expect_equal(nc$mm,  un$mm  - mean(un$mm))
  expect_equal(nc$msa, un$msa - mean(un$msa))
})

test_that("calculate_loglik_lrmsd_i_msa matches a frozen reference value", {
  # Frozen literal: captured once from the current implementation. NOT recomputed
  # here -- this catches a real change in the likelihood math or the pdb_site->i
  # join, not a re-derivation. (Updated 2026-06-24 when sigma was corrected from
  # sd(residuals) to the profile MLE sqrt(mean(residuals^2)); the value shifted by
  # a constant ~+0.0011 and the (a1,a2) argmax is unchanged.)
  pp <- znb_spm
  ll <- msamodel:::calculate_loglik_lrmsd_i_msa(pp, znb_profile$pdb_site, znb_profile$lrmsd_i_obs, a1 = 2, a2 = 5)
  expect_equal(ll, -184.3230779142, tolerance = 1e-6)
})
