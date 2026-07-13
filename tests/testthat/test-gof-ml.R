# Goodness-of-fit accessors for the ML arm (gof_lrmsd_i_msa_ml / gof_lrmsd_n_msa_ml).
# Scope is deliberately narrow: the D2/AIC/BIC arithmetic is a pure three-line function
# whose correctness was proved once in the scratchpad (the D2 == 1 - Var(resid)/Var(obs)
# identity and the AIC/BIC formulas) -- that is a one-time check, NOT a permanent
# regression guard, so it is not re-encoded here. What is kept:
#   - the fail-loud validation branch (a real stop() that must not silently return NA),
#   - a smoke test that the fitter -> accessor pipeline is wired (primitives stored,
#     row shape flows through) on both axes.

test_that("gof accessors fail loud on a fit missing the primitives", {
  # A raw pre-GoF fit (no deviance/null_deviance/nobs/k) must error, not return NA.
  raw <- list(a1 = 1, a2 = 1, logLik = -100)
  expect_error(gof_lrmsd_i_msa_ml(raw), "goodness-of-fit primitives")
  expect_error(gof_lrmsd_i_msa_ml(raw), "fit_lrmsd_i_msa_ml")   # names its producer
  expect_error(gof_lrmsd_i_msa_ml(42), "goodness-of-fit primitives")  # not-a-list
  expect_error(gof_lrmsd_n_msa_ml(raw), "fit_lrmsd_n_msa_ml")   # mode names its own
})

test_that("fitter -> gof accessor is wired on both axes (smoke)", {
  cols <- c("D2", "AIC", "BIC", "logLik", "deviance", "null_deviance", "nobs", "k")

  gi <- gof_lrmsd_i_msa_ml(fit_lrmsd_i_msa_ml(preprocess_spm(znb_spm), znb_profile))
  expect_named(gi, cols)
  expect_true(all(vapply(gi, is.finite, logical(1))))

  gn <- gof_lrmsd_n_msa_ml(fit_lrmsd_n_msa_ml(preprocess_spm_mode(znb_spm), znb_profile_n))
  expect_named(gn, cols)
  expect_true(all(vapply(gn, is.finite, logical(1))))
})
