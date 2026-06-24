# Drift-guards for the mode ML fit arm (fit_lrmsd_n_msa_ml) and its objective
# (calculate_loglik_lrmsd_n_msa). Regression / contract tests: once trusted, later
# edits must not change results or failure behaviour. They are NOT claims that the
# model is scientifically correct. The fit target is the SYNTHETIC znb_profile_n.

test_that("fit_lrmsd_n_msa_ml returns the documented list shape", {
  pp <- preprocess_spm_mode(znb_spm)
  ml <- fit_lrmsd_n_msa_ml(pp, znb_profile_n)

  expect_named(ml, c("a1", "a2", "logLik", "sigma_hat", "cov",
                     "se_a1", "se_a2", "convergence", "par_fit"))
  expect_length(ml$a1, 1L)
  expect_length(ml$a2, 1L)
  expect_equal(dim(ml$cov), c(2L, 2L))
  expect_equal(ml$cov, t(ml$cov))                 # symmetric
  expect_true(all(is.finite(ml$cov)))
  expect_true(ml$se_a1 > 0 && ml$se_a2 > 0)
  expect_true(is.finite(ml$logLik))
  expect_equal(ml$convergence, 0L)                # L-BFGS-B converged
  expect_equal(ml$par_fit[["b"]], log2(ml$a2 + 1))
})

test_that("mode ML estimate maximises the likelihood (independent grid route)", {
  # Independent route: a fine grid-search max of the SAME mode likelihood, in the
  # same (a1, b) coordinates. The optimiser's logLik must be at least as high as the
  # grid's, and its (a1, a2) must land in the same basin.
  pp <- preprocess_spm_mode(znb_spm)
  ml <- fit_lrmsd_n_msa_ml(pp, znb_profile_n)

  a1g <- seq(0, 10, length.out = 81)
  bg  <- seq(0, 13, length.out = 81)
  G   <- expand.grid(a1 = a1g, b = bg)
  ll  <- apply(G, 1L, function(r)
    calculate_loglik_lrmsd_n_msa(pp, znb_profile_n, a1 = r[["a1"]], a2 = 2^r[["b"]] - 1))
  best <- G[which.max(ll), ]

  expect_gte(ml$logLik, max(ll) - 1e-8)
  expect_lt(abs(ml$a1 - best$a1), a1g[2] - a1g[1])
  expect_lt(abs(ml$par_fit[["b"]] - best$b), bg[2] - bg[1])
})

test_that("fit_lrmsd_n_msa_ml matches frozen reference values", {
  # Frozen literals captured once from the current implementation (devtools state,
  # 2026-06-24) on the synthetic znb_profile_n. Catches drift in the mode fit math /
  # optimiser path, not a re-derivation.
  pp <- preprocess_spm_mode(znb_spm)
  ml <- fit_lrmsd_n_msa_ml(pp, znb_profile_n)

  expect_equal(ml$a1,        0.449221, tolerance = 1e-4)
  expect_equal(ml$a2,       40.819573, tolerance = 1e-3)
  expect_equal(ml$logLik, -155.553290, tolerance = 1e-4)
  expect_equal(ml$sigma_hat,  0.304371, tolerance = 1e-4)
  expect_equal(ml$se_a1,      0.057756, tolerance = 1e-3)
  expect_equal(ml$se_a2,      6.313771, tolerance = 1e-2)
})

test_that("fit_lrmsd_n_msa_ml validates box bounds and init (fail loud)", {
  pp <- preprocess_spm_mode(znb_spm)
  expect_error(fit_lrmsd_n_msa_ml(pp, znb_profile_n, a1_range = c(5)),
               "a1_range must be a vector of length 2")
  expect_error(fit_lrmsd_n_msa_ml(pp, znb_profile_n, a1_range = c(10, 0)),
               "a1_range must be a vector of length 2 with min < max")
  expect_error(fit_lrmsd_n_msa_ml(pp, znb_profile_n, log2_a2_plus1_range = c(13, 0)),
               "log2_a2_plus1_range must be a vector of length 2 with min < max")
  expect_error(fit_lrmsd_n_msa_ml(pp, znb_profile_n, init = c(1)),
               "init must be a length-2 numeric")
  expect_error(fit_lrmsd_n_msa_ml(pp, znb_profile_n, init = c(99, 1)),
               "init must lie within the box")
})

test_that("unknown mode index in observed_data is an error, not a silent drop", {
  # The n-coverage guard, mode analogue of the site pdb_site guard.
  pp <- preprocess_spm_mode(znb_spm)
  bad <- znb_profile_n
  bad$n[1] <- 999999L
  expect_error(
    calculate_loglik_lrmsd_n_msa(pp, bad, a1 = 2, a2 = 5),
    "mode index\\(es\\) not present in the model"
  )
  expect_error(fit_lrmsd_n_msa_ml(pp, bad), "not present in the model")
})

test_that("calculate_loglik_lrmsd_n_msa is invariant to a constant shift in lrmsd_n_obs", {
  # Both profiles are mean-centered, so adding a constant to the observed column
  # leaves the log-likelihood unchanged. Independent property, not a re-derivation.
  pp <- preprocess_spm_mode(znb_spm)
  ll1 <- calculate_loglik_lrmsd_n_msa(pp, znb_profile_n, a1 = 2, a2 = 5)
  shifted <- znb_profile_n
  shifted$lrmsd_n_obs <- shifted$lrmsd_n_obs + 3.7
  ll2 <- calculate_loglik_lrmsd_n_msa(pp, shifted, a1 = 2, a2 = 5)
  expect_equal(ll1, ll2)
})

test_that("znb_profile_n is reproducible from its seeded recipe", {
  # Determinism guard (hard project rule): re-derive the synthetic fixture from its
  # data-raw recipe and compare to the embedded data, to machine precision.
  ppm <- preprocess_spm_mode(znb_spm)
  site_fit <- fit_lrmsd_i_msa_ml(preprocess_spm(znb_spm), znb_profile)
  truth <- calculate_dr2_n_msa(ppm, site_fit$a1, site_fit$a2)
  set.seed(2025)
  expected_obs <- log(sqrt(truth$dr2_n)) + rnorm(nrow(truth), 0, 0.30)
  expect_equal(znb_profile_n$lrmsd_n_obs, expected_obs, tolerance = 1e-12)
  expect_equal(znb_profile_n$n, truth$n)
})
