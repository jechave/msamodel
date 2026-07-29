# ML delta-method predictors (predict_*_ml, both axes). Scope is deliberately narrow.
# The uncentred delta arithmetic (means == forward map, band == g +/- z*se) was proved
# once in the scratchpad -- a one-time check, not re-encoded here. Kept / added:
#   - the fail-loud validation branch (a real stop() on a fit missing a1/a2/cov),
#   - a wiring smoke test (fitter -> predictor returns the right-named row) for all ten
#     predictors, which can fail if a forward map or the fit contract drifts,
#   - the lrmsd-vs-nlrmsd band split: the centred (nlrmsd) band must use the
#     column-centred gradient, so it is NOT the uncentred (lrmsd) band. This is the
#     regression the split fixes; guarded here with a negative control below.

test_that("predict_*_ml fail loud on a fit missing the delta-method primitives", {
  pp  <- znb_spm
  raw <- list(a1 = 1, a2 = 1)                 # no cov
  expect_error(predict_lrmsd_i_msa_ml(raw, pp), "2x2 cov")
  expect_error(predict_lrmsd_i_msa_ml(raw, pp), "fit_lrmsd_i_msa_ml")

  bad_cov <- list(a1 = 1, a2 = 1, cov = matrix(0, 3, 3))   # wrong-shape cov
  expect_error(predict_nlrmsd_i_nested_models_ml(bad_cov, pp), "2x2 cov")

  expect_error(predict_nlrmsd_n_msa_decomposition_ml(42, znb_spm),
               "2x2 cov")
  # bad level on an otherwise-valid fit
  ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)
  expect_error(predict_nlrmsd_i_msa_ml(ml, pp, level = 1.5), "level must be")
})

test_that("fitter -> predict_*_ml is wired on both axes (smoke)", {
  pp  <- znb_spm
  ml  <- fit_lrmsd_i_msa_ml(pp, znb_profile)
  ppm <- znb_spm
  mln <- fit_lrmsd_n_msa_ml(ppm, znb_profile_n)

  bandcols <- function(prefix) as.vector(t(outer(prefix, c("_mean", "_lower", "_upper"), paste0)))

  # site: uncentred / centred profile
  expect_named(predict_lrmsd_i_msa_ml(ml, pp),
               c("i", "pdb_site", bandcols("lrmsd_i_msa")))
  expect_named(predict_nlrmsd_i_msa_ml(ml, pp),
               c("i", "pdb_site", bandcols("nlrmsd_i_msa")))
  # site: nested models
  expect_named(predict_lrmsd_i_nested_models_ml(ml, pp),
               c("i", "pdb_site", bandcols(c("lrmsd_i_mm", "lrmsd_i_ms",
                                             "lrmsd_i_ma", "lrmsd_i_msa"))))
  expect_named(predict_nlrmsd_i_nested_models_ml(ml, pp),
               c("i", "pdb_site", bandcols(c("nlrmsd_i_mm", "nlrmsd_i_ms",
                                             "nlrmsd_i_ma", "nlrmsd_i_msa"))))
  # site: decomposition (centred-only)
  expect_named(predict_nlrmsd_i_msa_decomposition_ml(ml, pp),
               c("i", "pdb_site", bandcols(c("nphi_mut", "nphi_stab", "nphi_act"))))

  # mode (no pdb_site)
  expect_named(predict_lrmsd_n_msa_ml(mln, ppm),
               c("n", bandcols("lrmsd_n_msa")))
  expect_named(predict_nlrmsd_n_msa_ml(mln, ppm),
               c("n", bandcols("nlrmsd_n_msa")))
  expect_named(predict_lrmsd_n_nested_models_ml(mln, ppm),
               c("n", bandcols(c("lrmsd_n_mm", "lrmsd_n_ms",
                                 "lrmsd_n_ma", "lrmsd_n_msa"))))
  expect_named(predict_nlrmsd_n_nested_models_ml(mln, ppm),
               c("n", bandcols(c("nlrmsd_n_mm", "nlrmsd_n_ms",
                                 "nlrmsd_n_ma", "nlrmsd_n_msa"))))
  expect_named(predict_nlrmsd_n_msa_decomposition_ml(mln, ppm),
               c("n", bandcols(c("nphi_mut", "nphi_stab", "nphi_act"))))

  # every band value finite, and lower <= mean <= upper on the full profile
  p <- predict_lrmsd_i_msa_ml(ml, pp)
  expect_true(all(is.finite(p$lrmsd_i_msa_mean)))
  expect_true(all(p$lrmsd_i_msa_lower <= p$lrmsd_i_msa_mean &
                  p$lrmsd_i_msa_mean  <= p$lrmsd_i_msa_upper))
})

test_that("nlrmsd band uses the CENTRED gradient, not the lrmsd (raw) band", {
  # The regression the split fixes: the centred (nlrmsd) band must be computed from the
  # column-centred gradient g - mean_S(g), so it is strictly different from -- here,
  # narrower than -- the uncentred (lrmsd) band. A vertical shift would leave the width
  # equal; that is the bug this guards against. Pinned to the PARAMETER arm: the
  # centred-vs-raw-gradient property is a statement about the parameter band (the SPM
  # arm has its own separate centring, tested in test-predict-uncertainty.R).
  pp <- znb_spm
  ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)
  z  <- stats::qnorm(0.975)

  pl <- predict_lrmsd_i_msa_ml(ml, pp, uncertainty = "parameter")
  pn <- predict_nlrmsd_i_msa_ml(ml, pp, uncertainty = "parameter")
  w_lrmsd  <- pl$lrmsd_i_msa_upper  - pl$lrmsd_i_msa_mean
  w_nlrmsd <- pn$nlrmsd_i_msa_upper - pn$nlrmsd_i_msa_mean

  # they must NOT be equal (negative control: if the code shifted the lrmsd band by a
  # constant instead of centring the gradient, these widths would be identical)
  expect_gt(max(abs(w_lrmsd - w_nlrmsd)), 1e-6)
  # and the centred band is narrower on average
  expect_lt(mean(w_nlrmsd), mean(w_lrmsd))

  # centred mean == uncentred mean minus its own mean over the full model support
  expect_equal(pn$nlrmsd_i_msa_mean,
               pl$lrmsd_i_msa_mean - mean(pl$lrmsd_i_msa_mean))

  # cross-check the nlrmsd width against a direct centred-gradient sandwich
  t_hat <- c(ml$a1, log2(ml$a2 + 1))
  f <- function(t) calculate_lrmsd_i_msa(pp, t[1], 2^t[2] - 1)$lrmsd_i_msa
  h <- 1e-5
  J <- sapply(1:2, function(j) {
    tp <- t_hat; tp[j] <- tp[j] + h; tm <- t_hat; tm[j] <- tm[j] - h
    (f(tp) - f(tm)) / (2 * h)
  })
  Jc <- sweep(J, 2, colMeans(J))
  se_ref <- unname(sqrt(rowSums((Jc %*% ml$cov) * Jc)))  # J carries dr2 col names; band doesn't
  expect_equal(w_nlrmsd, z * se_ref)
})

test_that("MM's centred PARAMETER band is exactly zero-width (constant gradient)", {
  # MM = (a1,a2)=(0,0) is fixed, so its parameter gradient is identically zero and
  # centring leaves it zero: the nlrmsd_i_mm PARAMETER band has no width. Pinned to
  # uncertainty="parameter" -- under the default "both" MM's band is nonzero (its SPM
  # arm), which is tested in test-predict-uncertainty.R. Negative control: nlrmsd_i_msa
  # (the full model) has a nonzero parameter band, so a "0 == 0" tautology cannot pass.
  pp <- znb_spm
  ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)
  nn <- predict_nlrmsd_i_nested_models_ml(ml, pp, uncertainty = "parameter")

  expect_equal(nn$nlrmsd_i_mm_upper, nn$nlrmsd_i_mm_lower)   # zero-width
  expect_gt(max(abs(nn$nlrmsd_i_msa_upper - nn$nlrmsd_i_msa_lower)), 1e-6)  # not vacuous
})
