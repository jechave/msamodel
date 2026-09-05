# Structural invariants.
#
# These hold for ANY valid input, at any (a1, a2), and follow from how the quantities are
# defined rather than from what this particular protein produced. That makes them the one
# part of the suite that survives a re-freeze: if a snapshot is deliberately accepted
# after an intended change, these still constrain the new values.
#
# A drift pin says "this number is what it was". An invariant says "these numbers are
# consistent with each other". The second catches a class the first cannot: a change that
# moves every value coherently -- a re-freeze accepted too readily, say -- leaves the
# pins green by construction but breaks an identity if the change was wrong.

test_that("the decomposition sums exactly to the profile it splits", {
  # phi_mut + phi_stab + phi_act == the full-model profile, by construction of the
  # sequential split. Exact, not approximate: it is a telescoping sum.
  # Negative control: swapping a minus for a plus in any phi term breaks it.
  for (metric in c("lrmsd", "nlrmsd")) {
    d <- calculate_decomposition(spm, 1.3, 0.7, metric)
    pre <- if (metric == "lrmsd") "" else "n"
    for (axis in c("site", "mode")) {
      tb <- d[[axis]]
      total <- tb[[paste0(pre, "phi_mut")]] + tb[[paste0(pre, "phi_stab")]] +
               tb[[paste0(pre, "phi_act")]]
      # Measured max relative deviation 3.8e-14 (round-off from the telescoping sum);
      # 1e-12 sits above that and far below anything a sign or term error would produce.
      expect_equal(total, tb[[paste0(pre, "lrmsd_msa")]],
                   tolerance = 1e-12, info = paste(metric, axis))
    }
  }
})

test_that("the mutation contribution IS the mutation-only nested model", {
  # phi_mut is defined as the MM variant, so the two columns must be identical -- value
  # AND standard error. A mis-wired arm in predict_se.R would separate them.
  ml <- fit_lrmsd_msa_site(spm, obs_site$pdb_site, obs_site$lrmsd_obs)
  pd <- predict_decomposition(ml, spm, "nlrmsd")
  expect_identical(pd$site$nphi_mut,    pd$site$nlrmsd_mm)
  expect_identical(pd$site$nphi_mut_se, pd$site$nlrmsd_mm_se)
  expect_identical(pd$mode$nphi_mut,    pd$mode$nlrmsd_mm)
  expect_identical(pd$mode$nphi_mut_se, pd$mode$nlrmsd_mm_se)
})

test_that("a centred profile has mean zero and differs from its uncentred twin by that mean", {
  # nlrmsd = lrmsd - mean(lrmsd) over the full model support.
  # Negative control: returning the uncentred profile unchanged fails both assertions,
  # since the profile mean is far from zero.
  p_l <- calculate_profiles(spm, 1.3, 0.7, "lrmsd")$site$lrmsd_msa
  p_n <- calculate_profiles(spm, 1.3, 0.7, "nlrmsd")$site$nlrmsd_msa
  # Measured: p_n - (p_l - mean(p_l)) is exactly zero (same computation, two routes),
  # so this is expect_identical -- no slack at all. mean(p_n) is not exactly 0 (it is
  # ~1e-16, the round-off of summing the centred values), so that one keeps a tolerance.
  expect_equal(mean(p_n), 0, tolerance = 1e-12)
  expect_identical(p_n, p_l - mean(p_l))
})

test_that("switching a selection pressure off reproduces the corresponding nested model", {
  # Setting a1 = 0 must give exactly the activity-only variant, and a2 = 0 the
  # stability-only one -- the nested models ARE the full model at those points.
  # This ties calculate_profiles and calculate_decomposition together through the
  # meaning of the parameters, not through a shared implementation detail.
  d <- calculate_decomposition(spm, 1.3, 0.7, "lrmsd")$site
  # Measured: all three agree exactly (bit-identical), so no tolerance is warranted.
  expect_identical(calculate_profiles(spm, 0, 0.7, "lrmsd")$site$lrmsd_msa, d$lrmsd_ma)
  expect_identical(calculate_profiles(spm, 1.3, 0, "lrmsd")$site$lrmsd_msa, d$lrmsd_ms)
  expect_identical(calculate_profiles(spm, 0, 0,   "lrmsd")$site$lrmsd_msa, d$lrmsd_mm)
})

test_that("the fit object is internally consistent", {
  ml <- fit_lrmsd_msa_site(spm, obs_site$pdb_site, obs_site$lrmsd_obs)

  expect_equal(ml$convergence, 0L)              # L-BFGS-B converged
  expect_equal(ml$cov, t(ml$cov))               # a covariance matrix is symmetric
  expect_true(ml$se_a1 > 0 && ml$se_a2 > 0)
  expect_equal(ml$gof$logLik, ml$logLik)        # the one field carried in both places

  # D2 is deviance-explained; deviance ties the profiled scale to the residual sum of
  # squares. Both are recomputed here from the primitives stored beside them, so a
  # clamped D2 or a sigma_hat that is sd() rather than the MLE goes red.
  #
  # D2 is expect_identical: the right-hand side is literally the expression the fitter
  # evaluates, so the two are the same arithmetic on the same doubles.
  #
  # deviance is NOT. The fitter computes `sum(resid^2)` and `sqrt(mean(resid^2))`
  # separately, so recovering one from the other runs a division, a sqrt, a square and
  # a multiply that `sum()` never did -- mathematically an identity, in floating point
  # a round trip. It was expect_identical until 2026-09-05, when CI failed on macOS and
  # Windows at a relative difference of 4e-15 while passing on Linux; over 1000 random
  # samples the bit-exact version fails 56% of the time, so the old comment's "each
  # side is the same arithmetic" was simply wrong and passing here had been luck.
  # expect_equal's default tolerance (~1.5e-8) still catches what this line is for:
  # a sigma_hat built from sd() is off by 5e-3, five orders of magnitude above it.
  expect_identical(ml$gof$D2, 1 - ml$gof$deviance / ml$gof$null_deviance)
  expect_equal(ml$gof$deviance, ml$gof$nobs * ml$gof$sigma_hat^2)
  expect_lt(ml$gof$D2, 1)
})

test_that("the two axes describe the same displacement", {
  # The site and mode forms are one displacement written in two bases, so for every
  # mutant the totals agree. This is the one check that the scan's two output matrices
  # are consistent with each other -- an axis-assembly bug that filled dr2mat_mode from
  # the wrong mutants would break it while leaving each matrix individually plausible.
  #
  # This is an EXACT identity -- one displacement expressed in two orthonormal bases --
  # so it should hold to floating-point round-off and nothing looser. The default
  # tolerance is far too slack for it: measured, it lets a single cell of dr2mat_mode be
  # scaled by 1.001 without failing.
  #
  # Measured agreement over the 1070 mutants: max absolute difference 3.2e-15, max
  # relative difference 3.0e-14 -- about 134x machine epsilon, i.e. ordinary accumulated
  # round-off over the 315 summed terms. 1e-12 sits safely above that while staying ~10
  # orders of magnitude tighter than the default, so a real assembly or alignment error
  # cannot hide under it.
  expect_equal(rowSums(spm$dr2mat_site), rowSums(spm$dr2mat_mode), tolerance = 1e-12)
})
