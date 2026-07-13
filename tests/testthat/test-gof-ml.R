# Goodness-of-fit accessors for the ML arm (gof_lrmsd_i_msa_ml / gof_lrmsd_n_msa_ml).
# The GoF numbers are pure arithmetic on the fit's stored primitives, so the formula
# tests use a HAND-BUILT fit-shaped list with known primitives -- the expected D2/AIC/
# BIC are computed from those literals independently of the accessor, not re-derived
# by the code under test.

test_that("gof_lrmsd_i_msa_ml returns the documented glance-style row", {
  pp <- preprocess_spm(znb_spm)
  ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)
  g  <- gof_lrmsd_i_msa_ml(ml)

  expect_s3_class(g, "tbl_df")
  expect_equal(nrow(g), 1L)
  expect_named(g, c("D2", "AIC", "BIC", "logLik",
                    "deviance", "null_deviance", "nobs", "k"))
  expect_true(all(vapply(g, is.finite, logical(1))))
})

test_that("D2 / AIC / BIC match their definitions from stored primitives", {
  # Hand-built fit with KNOWN primitives; expected values computed here by hand, not
  # by the accessor. logLik = -100, deviance = 30, null_deviance = 120, nobs = 50, k = 3.
  fit <- list(logLik = -100, deviance = 30, null_deviance = 120, nobs = 50L, k = 3L)
  g   <- gof_lrmsd_i_msa_ml(fit)

  expect_equal(g$D2,  1 - 30 / 120)          # = 0.75
  expect_equal(g$AIC, -2 * (-100) + 2 * 3)   # = 206
  expect_equal(g$BIC, -2 * (-100) + 3 * log(50))
  # primitives pass through unchanged
  expect_equal(g$logLik, -100)
  expect_equal(g$deviance, 30)
  expect_equal(g$null_deviance, 120)
  expect_equal(g$nobs, 50L)
  expect_equal(g$k, 3L)
})

test_that("D2 is capped at 1 but unbounded below (returned unclamped)", {
  # deviance == null_deviance -> exactly the flat null -> D2 == 0.
  expect_equal(
    gof_lrmsd_i_msa_ml(list(logLik = -1, deviance = 10, null_deviance = 10,
                            nobs = 20L, k = 3L))$D2,
    0
  )
  # deviance > null_deviance (prediction worse than flat) -> D2 < 0, NOT floored.
  expect_lt(
    gof_lrmsd_i_msa_ml(list(logLik = -1, deviance = 25, null_deviance = 10,
                            nobs = 20L, k = 3L))$D2,
    0
  )
  # deviance == 0 (perfect) -> D2 == 1 (the upper cap).
  expect_equal(
    gof_lrmsd_i_msa_ml(list(logLik = -1, deviance = 0, null_deviance = 10,
                            nobs = 20L, k = 3L))$D2,
    1
  )
})

test_that("gof accessors fail loud on a fit missing the primitives", {
  # A raw pre-GoF fit (no deviance/null_deviance/nobs/k) must error, not return NA.
  raw <- list(a1 = 1, a2 = 1, logLik = -100)
  expect_error(gof_lrmsd_i_msa_ml(raw),
               "goodness-of-fit primitives")
  expect_error(gof_lrmsd_i_msa_ml(raw),
               "fit_lrmsd_i_msa_ml")
  # not-a-list
  expect_error(gof_lrmsd_i_msa_ml(42), "goodness-of-fit primitives")
  # mode accessor names its own producer
  expect_error(gof_lrmsd_n_msa_ml(raw), "fit_lrmsd_n_msa_ml")
})

test_that("gof_lrmsd_n_msa_ml (mode) returns the same-shaped row", {
  ppm <- preprocess_spm_mode(znb_spm)
  ml  <- fit_lrmsd_n_msa_ml(ppm, znb_profile_n)
  g   <- gof_lrmsd_n_msa_ml(ml)

  expect_named(g, c("D2", "AIC", "BIC", "logLik",
                    "deviance", "null_deviance", "nobs", "k"))
  expect_true(all(vapply(g, is.finite, logical(1))))
})
