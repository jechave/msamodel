# Goodness of fit is attached by the fitter as `fit$gof` -- there is no exported
# accessor. Scope is deliberately narrow: the D2/AIC/BIC arithmetic is a pure
# three-line function whose correctness was proved once in the scratchpad (the
# D2 == 1 - Var(resid)/Var(obs) identity and the AIC/BIC formulas) -- a one-time check,
# NOT a permanent regression guard, so it is not re-encoded here.
#
# The old "fail loud on a fit missing the primitives" test was DELETED with
# validate_gof_fit(): the fitter now always attaches $gof, so a fit without it is
# unrepresentable. A test for an unreachable state cannot fail for a real reason.

test_that("both fitters attach a well-formed $gof (both axes)", {
  cols <- c("D2", "AIC", "BIC", "logLik", "deviance", "null_deviance",
            "nobs", "k", "sigma_hat")

  gi <- fit_lrmsd_msa_site(znb_spm, znb_profile$pdb_site, znb_profile$lrmsd_i_obs)$gof
  expect_named(gi, cols)
  expect_equal(nrow(gi), 1L)
  expect_true(all(vapply(gi, is.finite, logical(1))))

  gn <- fit_lrmsd_msa_mode(znb_spm, znb_profile_n$n, znb_profile_n$lrmsd_n_obs)$gof
  expect_named(gn, cols)
  expect_equal(nrow(gn), 1L)
  expect_true(all(vapply(gn, is.finite, logical(1))))

  # nobs is the MATCHED count, not the model's full support: znb_profile covers a
  # subset of the 228 sites, whereas every one of the 678 modes is observed.
  expect_equal(gi$nobs, nrow(znb_profile))
  expect_equal(gn$nobs, nrow(znb_profile_n))
  expect_equal(gi$k, 3L)                    # a1, a2, profiled sigma
})

test_that("D2 relates deviance to the null the way the formula claims", {
  # Independent route: D2 is recomputed from the two deviances stored beside it, so a
  # mis-wired numerator/denominator (or a clamped D2) goes red. Negative control:
  # returning deviance/null_deviance instead of 1 - that ratio fails this.
  g <- fit_lrmsd_msa_site(znb_spm, znb_profile$pdb_site, znb_profile$lrmsd_i_obs)$gof
  expect_equal(g$D2, 1 - g$deviance / g$null_deviance)
  expect_lt(g$D2, 1)                        # at most 1, never exactly (noise remains)

  # deviance == nobs * sigma_hat^2 ties the profiled scale to the residual sum of
  # squares; an inconsistent sigma_hat (e.g. sd() instead of the MLE) breaks this.
  expect_equal(g$deviance, g$nobs * g$sigma_hat^2)
})
