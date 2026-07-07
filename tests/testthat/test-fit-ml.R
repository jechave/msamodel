# Drift-guards for the ML fit arm (fit_lrmsd_i_msa_ml). These are regression / contract
# tests: once the estimator is trusted, later edits must not change its results
# or its failure behaviour. They are NOT claims that the model is scientifically
# correct.

test_that("fit_lrmsd_i_msa_ml returns the documented list shape", {
  pp <- preprocess_spm(znb_spm)
  ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)

  expect_named(ml, c("a1", "a2", "logLik", "sigma_hat", "cov",
                     "se_a1", "se_a2", "convergence"))
  expect_length(ml$a1, 1L)
  expect_length(ml$a2, 1L)
  expect_equal(dim(ml$cov), c(2L, 2L))
  expect_equal(ml$cov, t(ml$cov))                 # symmetric
  expect_true(all(is.finite(ml$cov)))
  expect_true(ml$se_a1 > 0 && ml$se_a2 > 0)
  expect_true(is.finite(ml$logLik))
  expect_equal(ml$convergence, 0L)                # L-BFGS-B converged
})

test_that("ML sits at a local max of its own objective (consistency)", {
  # Replaces the retired 81x81 grid search (~58 s). Evaluate the SAME objective the
  # ML fit optimises on a tiny grid centred on the fit's own (a1, a2): optim's logLik
  # must be no worse than any nearby point, i.e. it sits at a local max. Self-
  # contained -- depends only on fit_lrmsd_i_msa_ml + calculate_loglik_lrmsd_i_msa,
  # NOT on AGQ or any other function. Honest scope: this is a local-max CONSISTENCY
  # check on the ML path, NOT an independent-correctness proof (the grid recomputes
  # the same objective). WHERE the optimum is, is pinned by the frozen-reference test
  # below; THAT it is a max there, is pinned here.
  pp <- preprocess_spm(znb_spm)
  ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)

  a1g <- ml$a1 + c(-0.1, -0.05, 0, 0.05, 0.1)
  bg  <- log2(ml$a2 + 1) + c(-0.3, -0.15, 0, 0.15, 0.3)
  G   <- expand.grid(a1 = a1g, b = bg)
  ll  <- apply(G, 1L, function(r)
    msamodel:::calculate_loglik_lrmsd_i_msa(pp, znb_profile, a1 = r[["a1"]], a2 = 2^r[["b"]] - 1))

  expect_gte(ml$logLik, max(ll) - 1e-8)
})

test_that("fit_lrmsd_i_msa_ml matches frozen reference values", {
  # Frozen literals captured once from the current implementation (devtools state,
  # 2026-06-24). Catches drift in the fit math / optimiser path, not a re-derivation.
  pp <- preprocess_spm(znb_spm)
  ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)

  expect_equal(ml$a1,     0.458002, tolerance = 1e-4)
  expect_equal(ml$a2,    42.302191, tolerance = 1e-3)
  expect_equal(ml$logLik, -138.545817, tolerance = 1e-4)
  expect_equal(ml$sigma_hat, 0.447903, tolerance = 1e-4)
  expect_equal(ml$se_a1,  0.122275, tolerance = 1e-3)
  expect_equal(ml$se_a2, 10.186796, tolerance = 1e-2)
})

test_that("fit_lrmsd_i_msa_ml validates box bounds (fail loud)", {
  pp <- preprocess_spm(znb_spm)
  expect_error(fit_lrmsd_i_msa_ml(pp, znb_profile, a1_range = c(5)),
               "a1_range must be a vector of length 2")
  expect_error(fit_lrmsd_i_msa_ml(pp, znb_profile, a1_range = c(10, 0)),
               "a1_range must be a vector of length 2 with min < max")
  expect_error(fit_lrmsd_i_msa_ml(pp, znb_profile, log2_a2_plus1_range = c(13, 0)),
               "log2_a2_plus1_range must be a vector of length 2 with min < max")
  expect_error(fit_lrmsd_i_msa_ml(pp, znb_profile, init = c(1)),
               "init must be a length-2 numeric")
  expect_error(fit_lrmsd_i_msa_ml(pp, znb_profile, init = c(99, 1)),
               "init must lie within the box")
})

test_that("fit_lrmsd_i_msa_ml inherits the pdb_site contract from the likelihood", {
  pp <- preprocess_spm(znb_spm)

  # Unknown pdb_site -> error (not a silent drop), inherited from calculate_loglik_lrmsd_i_msa.
  bad <- znb_profile
  bad$pdb_site[1] <- 999999L
  expect_error(fit_lrmsd_i_msa_ml(pp, bad), "not present in the model")

  # Partial coverage still yields a finite fit.
  expect_lt(nrow(znb_profile), nrow(pp$site_map))   # genuinely a subset
  ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)
  expect_true(is.finite(ml$a1) && is.finite(ml$a2) && is.finite(ml$logLik))
})
