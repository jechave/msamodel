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

test_that("calculate_nlrmsd_i_msa centres the uncentred profile over the full support", {
  pp <- preprocess_spm(znb_spm)
  a1 <- 2; a2 <- 5
  nc <- calculate_nlrmsd_i_msa(pp, a1, a2)

  expect_named(nc, c("i", "nlrmsd_i_msa"))
  expect_equal(nrow(nc), 228L)

  # Independent route: the centred value is the uncentred lrmsd minus its own mean
  # over ALL model residues. Negative control: an uncentred return (nlrmsd == lrmsd)
  # would fail this, since the profile mean is far from zero.
  lr <- calculate_lrmsd_i_msa(pp, a1, a2)$lrmsd_i_msa
  expect_equal(nc$nlrmsd_i_msa, lr - mean(lr))
  expect_equal(mean(nc$nlrmsd_i_msa), 0)          # centred: mean is (machine) zero
})

test_that("calculate_nlrmsd_i_msa agrees with predict_nlrmsd_i_msa_ml's point profile", {
  # The forward map and the predictor must return the SAME centred profile (the whole
  # point of the twin). uncertainty='none' strips the band so only the mean remains.
  # Negative control: if the twin centred over the wrong support (or not at all), the
  # predictor -- which centres over the full model support -- would disagree.
  # Compare VALUES: the forward map carries the dr2 column names on the value vector
  # (existing convention), whereas the predictor's band assembler unname()s them --
  # an intended difference in the names attribute, not the numbers.
  pp <- preprocess_spm(znb_spm)
  ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)
  fwd  <- calculate_nlrmsd_i_msa(pp, ml$a1, ml$a2)
  pred <- predict_nlrmsd_i_msa_ml(ml, pp, uncertainty = "none")
  expect_equal(unname(fwd$nlrmsd_i_msa), pred$nlrmsd_i_msa_mean)
})

test_that("calculate_nlrmsd_i_msa_decomposition contributions sum to the centred profile", {
  pp <- preprocess_spm(znb_spm)
  a1 <- 2; a2 <- 5
  d <- calculate_nlrmsd_i_msa_decomposition(pp, a1, a2)

  expect_named(d, c("i", "pdb_site", "nphi_mut", "nphi_stab", "nphi_act"))
  # The three centred contributions reconstruct the centred full-model profile.
  # Independent route: nlrmsd_i_msa from its own twin, not from this decomposition.
  # Negative control: dropping the centring on any nphi term breaks the equality.
  prof <- calculate_nlrmsd_i_msa(pp, a1, a2)$nlrmsd_i_msa
  expect_equal(d$nphi_mut + d$nphi_stab + d$nphi_act, prof)
})

test_that("calculate_nlrmsd_i_nested_models centres each variant independently", {
  pp <- preprocess_spm(znb_spm)
  a1 <- 2; a2 <- 5
  nc <- calculate_nlrmsd_i_nested_models(pp, a1, a2)

  expect_named(nc, c("i", "pdb_site", "nlrmsd_i_mm", "nlrmsd_i_ms",
                     "nlrmsd_i_ma", "nlrmsd_i_msa"))
  # Each variant centred by ITS OWN mean (independent route via the uncentred twin).
  # Negative control: centring every column by a single shared mean would fail.
  un <- calculate_lrmsd_i_nested_models(pp, a1, a2)
  expect_equal(nc$nlrmsd_i_mm,  un$lrmsd_i_mm  - mean(un$lrmsd_i_mm))
  expect_equal(nc$nlrmsd_i_msa, un$lrmsd_i_msa - mean(un$lrmsd_i_msa))
})

test_that("calculate_loglik_lrmsd_i_msa matches a frozen reference value", {
  # Frozen literal: captured once from the current implementation. NOT recomputed
  # here -- this catches a real change in the likelihood math or the pdb_site->i
  # join, not a re-derivation. (Updated 2026-06-24 when sigma was corrected from
  # sd(residuals) to the profile MLE sqrt(mean(residuals^2)); the value shifted by
  # a constant ~+0.0011 and the (a1,a2) argmax is unchanged.)
  pp <- preprocess_spm(znb_spm)
  ll <- msamodel:::calculate_loglik_lrmsd_i_msa(pp, znb_profile, a1 = 2, a2 = 5)
  expect_equal(ll, -184.3230779142, tolerance = 1e-6)
})
