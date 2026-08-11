# Drift-guards for the mode ML fit arm (fit_lrmsd_msa_mode) and its objective
# (loglik_lrmsd_msa). Regression / contract tests: once trusted, later
# edits must not change results or failure behaviour. They are NOT claims that the
# model is scientifically correct. The fit target is the SYNTHETIC znb_profile_n.

test_that("fit_lrmsd_msa_mode returns the documented list shape", {
  pp <- znb_spm
  ml <- fit_lrmsd_msa_mode(pp, znb_profile_n$mode, znb_profile_n$lrmsd_obs)

  expect_named(ml, c("a1", "a2", "logLik", "cov", "se_a1", "se_a2",
                     "convergence", "gof", "call"))
  expect_identical(as.character(ml$call[[1]]), "fit_lrmsd_msa_mode")
  expect_length(ml$a1, 1L)
  expect_length(ml$a2, 1L)
  expect_equal(dim(ml$cov), c(2L, 2L))
  expect_equal(ml$cov, t(ml$cov))                 # symmetric
  expect_true(all(is.finite(ml$cov)))
  expect_true(ml$se_a1 > 0 && ml$se_a2 > 0)
  expect_true(is.finite(ml$logLik))
  expect_equal(ml$convergence, 0L)                # L-BFGS-B converged
})

test_that("mode ML sits at a local max of its own objective (consistency)", {
  # Replaces the retired 81x81 grid search (~65 s). Mode analogue of the site-arm
  # check in test-fit-ml.R: evaluate the SAME objective the mode fit optimises on a
  # tiny grid centred on the fit's own (a1, a2); optim's logLik must be no worse than
  # any nearby point. Self-contained -- depends only on fit_lrmsd_msa_mode +
  # loglik_lrmsd_msa. Honest scope: a local-max CONSISTENCY check, NOT
  # independent correctness (the grid recomputes the same objective). WHERE the
  # optimum is, is pinned by the frozen-reference test below; THAT it is a max there,
  # is pinned here.
  pp <- znb_spm
  ml <- fit_lrmsd_msa_mode(pp, znb_profile_n$mode, znb_profile_n$lrmsd_obs)

  a1g <- ml$a1 + c(-0.1, -0.05, 0, 0.05, 0.1)
  bg  <- log2(ml$a2 + 1) + c(-0.3, -0.15, 0, 0.15, 0.3)
  G   <- expand.grid(a1 = a1g, b = bg)
  ll  <- apply(G, 1L, function(r)
    msamodel:::loglik_lrmsd_msa(pp$dr2_njm, pp$energy_data, msamodel:::resolve_mode_obs(seq_len(ncol(pp$dr2_njm)), znb_profile_n$mode, znb_profile_n$lrmsd_obs), a1 = r[["a1"]], a2 = 2^r[["b"]] - 1))

  expect_gte(ml$logLik, max(ll) - 1e-8)
})

test_that("fit_lrmsd_msa_mode matches frozen reference values", {
  # Frozen literals captured once from the current implementation (devtools state,
  # 2026-06-24) on the synthetic znb_profile_n. Catches drift in the mode fit math /
  # optimiser path, not a re-derivation.
  pp <- znb_spm
  ml <- fit_lrmsd_msa_mode(pp, znb_profile_n$mode, znb_profile_n$lrmsd_obs)

  expect_equal(ml$a1,        0.449221, tolerance = 1e-4)
  expect_equal(ml$a2,       40.819573, tolerance = 1e-3)
  expect_equal(ml$logLik, -155.553290, tolerance = 1e-4)
  expect_equal(ml$gof$sigma_hat,  0.304371, tolerance = 1e-4)
  expect_equal(ml$se_a1,      0.057756, tolerance = 1e-3)
  expect_equal(ml$se_a2,      6.313771, tolerance = 1e-2)
})

test_that("fit_lrmsd_msa_mode validates box bounds and init (fail loud)", {
  pp <- znb_spm
  expect_error(fit_lrmsd_msa_mode(pp, znb_profile_n$mode, znb_profile_n$lrmsd_obs, a1_range = c(5)),
               "a1_range must be a vector of length 2")
  expect_error(fit_lrmsd_msa_mode(pp, znb_profile_n$mode, znb_profile_n$lrmsd_obs, a1_range = c(10, 0)),
               "a1_range must be a vector of length 2 with min < max")
  expect_error(fit_lrmsd_msa_mode(pp, znb_profile_n$mode, znb_profile_n$lrmsd_obs, log2_a2_plus1_range = c(13, 0)),
               "log2_a2_plus1_range must be a vector of length 2 with min < max")
  expect_error(fit_lrmsd_msa_mode(pp, znb_profile_n$mode, znb_profile_n$lrmsd_obs, init = c(1)),
               "init must be a length-2 numeric")
  expect_error(fit_lrmsd_msa_mode(pp, znb_profile_n$mode, znb_profile_n$lrmsd_obs, init = c(99, 1)),
               "init must lie within the box")
  expect_error(fit_lrmsd_msa_mode(pp, znb_profile_n$mode, znb_profile_n$lrmsd_obs, init = c(a1 = 1, b = 1)),
               "init must be an unnamed positional")
})

test_that("unknown mode index in observed_data is an error, not a silent drop", {
  # The n-coverage guard, mode analogue of the site pdb_site guard.
  pp <- znb_spm
  bad_mode <- znb_profile_n$mode
  bad_mode[1] <- 999999L
  expect_error(
    msamodel:::resolve_mode_obs(seq_len(ncol(pp$dr2_njm)), bad_mode,
                                znb_profile_n$lrmsd_obs),
    "mode index\\(es\\) not present in the model"
  )
  expect_error(fit_lrmsd_msa_mode(pp, bad_mode, znb_profile_n$lrmsd_obs),
               "not present in the model")
})

test_that("loglik_lrmsd_msa is invariant to a constant shift in lrmsd_obs", {
  # Both profiles are mean-centered, so adding a constant to the observed column
  # leaves the log-likelihood unchanged. Independent property, not a re-derivation.
  pp <- znb_spm
  ll1 <- msamodel:::loglik_lrmsd_msa(pp$dr2_njm, pp$energy_data, msamodel:::resolve_mode_obs(seq_len(ncol(pp$dr2_njm)), znb_profile_n$mode, znb_profile_n$lrmsd_obs), a1 = 2, a2 = 5)
  shifted <- znb_profile_n$lrmsd_obs + 3.7
  ll2 <- msamodel:::loglik_lrmsd_msa(
    pp$dr2_njm, pp$energy_data,
    msamodel:::resolve_mode_obs(seq_len(ncol(pp$dr2_njm)), znb_profile_n$mode, shifted),
    a1 = 2, a2 = 5)
  expect_equal(ll1, ll2)
})

test_that("znb_profile_n is reproducible from its seeded recipe", {
  # Determinism guard (hard project rule): re-derive the synthetic mode profile from
  # its data-raw recipe and compare to the shipped file.
  #
  # tolerance = 1e-12, not exact equality: the shipped profile is a CSV, and a CSV
  # round-trip is NOT bit-exact for doubles (write_csv emits ~15 significant digits,
  # so the last bit is lost -- measured at 8.9e-16 max here). That is a property of
  # the file format, not drift in the recipe. Real drift moves numbers far above
  # 1e-12; a failure here is a recipe change, not round-off.
  ppm <- znb_spm
  site_fit <- fit_lrmsd_msa_site(znb_spm, znb_profile$pdb_site, znb_profile$lrmsd_obs)
  dr2_n <- dr2_msa(ppm$dr2_njm, ppm$energy_data, site_fit$a1, site_fit$a2)
  set.seed(2025)
  expected_obs <- log(sqrt(dr2_n)) + rnorm(length(dr2_n), 0, 0.30)
  expect_equal(znb_profile_n$lrmsd_obs, expected_obs, tolerance = 1e-12)
  expect_equal(znb_profile_n$mode, seq_along(dr2_n))
})
