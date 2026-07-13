# ML delta-method predictors (predict_*_ml, both axes). Scope is deliberately narrow:
# the delta arithmetic (means == forward map, band == g +/- z*se, se == fit$se_a2) was
# proved once in the scratchpad -- a one-time correctness check, NOT a permanent
# regression guard, so it is not re-encoded here. Kept:
#   - the fail-loud validation branch (a real stop() on a fit missing a1/a2/cov),
#   - a wiring smoke test (fitter -> predictor returns the right-named row) for all six
#     predictors, which can fail if a forward map or the fit contract drifts.

test_that("predict_*_ml fail loud on a fit missing the delta-method primitives", {
  pp  <- preprocess_spm(znb_spm)
  raw <- list(a1 = 1, a2 = 1)                 # no cov
  expect_error(predict_lrmsd_i_msa_ml(raw, pp), "2x2 cov")
  expect_error(predict_lrmsd_i_msa_ml(raw, pp), "fit_lrmsd_i_msa_ml")

  bad_cov <- list(a1 = 1, a2 = 1, cov = matrix(0, 3, 3))   # wrong-shape cov
  expect_error(predict_lrmsd_i_nested_models_ml(bad_cov, pp), "2x2 cov")

  expect_error(predict_decomposition_n_msa_ml(42, preprocess_spm_mode(znb_spm)),
               "2x2 cov")
  # bad level on an otherwise-valid fit
  ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)
  expect_error(predict_lrmsd_i_msa_ml(ml, pp, level = 1.5), "level must be")
})

test_that("fitter -> predict_*_ml is wired on both axes (smoke)", {
  pp  <- preprocess_spm(znb_spm)
  ml  <- fit_lrmsd_i_msa_ml(pp, znb_profile)
  ppm <- preprocess_spm_mode(znb_spm)
  mln <- fit_lrmsd_n_msa_ml(ppm, znb_profile_n)

  bandcols <- function(prefix) as.vector(t(outer(prefix, c("_mean", "_lower", "_upper"), paste0)))

  # site
  expect_named(predict_lrmsd_i_msa_ml(ml, pp),
               c("i", "pdb_site", bandcols("lrmsd_i_msa"), bandcols("nlrmsd_i_msa")))
  expect_named(predict_lrmsd_i_nested_models_ml(ml, pp),
               c("i", "pdb_site", bandcols(c("lrmsd_i_mm", "lrmsd_i_ms",
                                             "lrmsd_i_ma", "lrmsd_i_msa"))))
  expect_named(predict_decomposition_i_msa_ml(ml, pp),
               c("i", "pdb_site", bandcols(c("phi_mut", "phi_stab", "phi_act"))))

  # mode (no pdb_site)
  expect_named(predict_lrmsd_n_msa_ml(mln, ppm),
               c("n", bandcols("lrmsd_n_msa"), bandcols("nlrmsd_n_msa")))
  expect_named(predict_lrmsd_n_nested_models_ml(mln, ppm),
               c("n", bandcols(c("lrmsd_n_mm", "lrmsd_n_ms",
                                 "lrmsd_n_ma", "lrmsd_n_msa"))))
  expect_named(predict_decomposition_n_msa_ml(mln, ppm),
               c("n", bandcols(c("phi_mut", "phi_stab", "phi_act"))))

  # every band value finite, and lower <= mean <= upper on the full profile
  p <- predict_lrmsd_i_msa_ml(ml, pp)
  expect_true(all(is.finite(p$lrmsd_i_msa_mean)))
  expect_true(all(p$lrmsd_i_msa_lower <= p$lrmsd_i_msa_mean &
                  p$lrmsd_i_msa_mean  <= p$lrmsd_i_msa_upper))
})
