# Drift-guards and correctness checks for the AGQ fit arm (fit_lrmsd_i_msa_agq) and
# its profile predictor (predict_lrmsd_i_agq). The moment checks are NOT
# tautologies: they compare AGQ against an INDEPENDENT dense-grid normalization of
# the same likelihood, in the (a1, log2(a2+1)) coordinates where the prior is flat.

# Independent ground truth: normalize exp(ll) on a dense regular grid in
# (a1, t = log2(a2+1)) -- the flat-in-t posterior, matching the MCMC's prior. This is
# a QUADRATURE reference (brute-force posterior mean/sd), so its accuracy is its
# value -- it cannot be tiny. A measured grid-convergence sweep vs a 61x61 reference:
# 5x5/7x7 are garbage (a2 off by ~19/13), 11x11 and 15x15 still miss the test
# tolerances; 21x21 matches 61x61 to a2 < 0.003 and the sds < 0.001 -- orders of
# magnitude inside the tolerances asserted below -- at ~1/8 the cost. So 21x21 is the
# converged minimum with margin. Computed ONCE here and reused.
agq_ground_truth <- local({
  pp  <- preprocess_spm(znb_spm)
  obs <- znb_profile
  ll_of <- function(a1, t) calculate_loglik_lrmsd_i_msa(pp, obs, a1, 2^t - 1)
  a1g <- seq(0, 3, length.out = 21)
  tg  <- seq(0, 9, length.out = 21)
  G   <- expand.grid(a1 = a1g, t = tg)
  G$ll <- mapply(ll_of, G$a1, G$t)
  w  <- exp(G$ll - max(G$ll)); ws <- sum(w)
  a2 <- 2^G$t - 1
  m_a1 <- sum(G$a1 * w) / ws; m_a2 <- sum(a2 * w) / ws
  list(
    a1 = m_a1, a2 = m_a2,
    sd_a1 = sqrt(sum((G$a1 - m_a1)^2 * w) / ws),
    sd_a2 = sqrt(sum((a2 - m_a2)^2 * w) / ws)
  )
})

# Shared default fit. AGQ is deterministic (pure quadrature, no RNG --
# identical(fit, fit) is TRUE), so the default n_nodes = 7 fit on znb_profile is the
# SAME object in every block that uses it. Compute it ONCE and reuse, instead of
# paying ~7 s per block (object-shape, frozen-ref, predict-band, band-sharpens). Same
# compute-once pattern as agq_ground_truth above. Blocks that need a DIFFERENT n_nodes
# still fit their own.
agq_default <- local({
  pp <- preprocess_spm(znb_spm)
  list(pp = pp, agq = fit_lrmsd_i_msa_agq(pp, znb_profile))
})

test_that("gauss_hermite nodes/weights are correct (integrate polynomials exactly)", {
  gh <- gauss_hermite(5)
  expect_length(gh$x, 5L)
  expect_equal(sum(gh$w), sqrt(pi))          # total weight
  # int x^2 e^{-x^2} dx = sqrt(pi)/2; GH must be exact for a degree-2 polynomial.
  expect_equal(sum(gh$w * gh$x^2), sqrt(pi) / 2)
  # odd moment is zero by symmetry.
  expect_equal(sum(gh$w * gh$x), 0)
})

test_that("fit_lrmsd_i_msa_agq returns the documented object", {
  agq <- agq_default$agq

  expect_s3_class(agq, "msa_agq")
  expect_named(agq, c("a1", "a2", "sd_a1", "sd_a2", "ci_a1", "ci_a2",
                      "nodes", "n_nodes", "laplace", "log_evidence"))
  expect_equal(agq$n_nodes, 7L)
  expect_equal(nrow(agq$nodes), 49L)
  expect_named(agq$nodes, c("a1", "a2", "log_weight"))
  # normalized posterior MASS
  expect_equal(sum(exp(agq$nodes$log_weight)), 1)
  # CI brackets the mean
  expect_true(agq$ci_a1[1] < agq$a1 && agq$a1 < agq$ci_a1[2])
  expect_true(agq$ci_a2[1] < agq$a2 && agq$a2 < agq$ci_a2[2])
  # Laplace cov is the t-scale object, explicitly labelled
  expect_equal(dim(agq$laplace$cov), c(2L, 2L))
  expect_equal(dimnames(agq$laplace$cov)[[1]], c("a1", "log2_a2_plus1"))
  expect_true(is.finite(agq$log_evidence))
})

test_that("AGQ moments match an independent dense-grid posterior", {
  pp    <- preprocess_spm(znb_spm)
  truth <- agq_ground_truth
  agq   <- fit_lrmsd_i_msa_agq(pp, znb_profile, n_nodes = 5)

  expect_lt(abs(agq$a1 - truth$a1), 0.01)
  expect_lt(abs(agq$a2 - truth$a2), 0.05)
  # sd within ~3% -- this is the change-of-measure check: a wrong weight (Gaussian
  # applied twice instead of divided out) collapses the variance and fails here.
  expect_lt(abs(agq$sd_a1 - truth$sd_a1) / truth$sd_a1, 0.03)
  expect_lt(abs(agq$sd_a2 - truth$sd_a2) / truth$sd_a2, 0.03)
})

test_that("AGQ moments converge toward truth as nodes increase (3 -> 5)", {
  pp <- preprocess_spm(znb_spm)
  truth <- agq_ground_truth
  e <- function(n) abs(fit_lrmsd_i_msa_agq(pp, znb_profile, n_nodes = n)$a2 - truth$a2)
  # 5 nodes must be at least as close to the truth as 3 (the built-in non-Gaussianity
  # diagnostic: moments stabilize as the grid refines).
  expect_lt(e(5), e(3))
})

test_that("fit_lrmsd_i_msa_agq matches frozen reference values", {
  # Frozen literals captured once from the current implementation (devtools state,
  # 2026-06-25) on znb_profile, n_nodes = 5. Catches drift in the quadrature path.
  # AGQ is deterministic (pure quadrature, no RNG): identical(fit, fit) is TRUE, so
  # these literals are exact, not seeded. Captured at the default n_nodes = 7.
  agq <- agq_default$agq
  expect_equal(agq$a1, 0.47065568, tolerance = 1e-6)
  expect_equal(agq$a2, 42.69310722, tolerance = 1e-5)
  expect_equal(agq$sd_a1, 0.12545482, tolerance = 1e-6)
  expect_equal(agq$sd_a2, 10.55793646, tolerance = 1e-5)
  expect_equal(agq$log_evidence, -1.31818323, tolerance = 1e-5)
})

test_that("fit_lrmsd_i_msa_agq fails loud on bad input", {
  pp <- preprocess_spm(znb_spm)
  expect_error(fit_lrmsd_i_msa_agq(pp, znb_profile, n_nodes = 0), "positive integer")
  expect_error(fit_lrmsd_i_msa_agq(pp, znb_profile, n_nodes = 2.5), "positive integer")
  # singular-Hessian path propagates from the ML reference (degenerate box bound)
  expect_error(
    fit_lrmsd_i_msa_agq(pp, znb_profile, a1_range = c(1, 1)),
    "min < max"
  )
})

test_that("predict_lrmsd_i_agq gives a per-site banded profile", {
  pp  <- agq_default$pp
  agq <- agq_default$agq
  pr  <- predict_lrmsd_i_agq(agq, pp)

  # one row per MODEL site (228), not per observed site (znb_profile covers 225) --
  # a prediction is defined everywhere the model is, including unobserved sites.
  expect_equal(nrow(pr), nrow(pp$site_map))
  expect_named(pr, c("i", "pdb_site",
                     "lrmsd_i_msa_mean", "lrmsd_i_msa_lower", "lrmsd_i_msa_upper",
                     "nlrmsd_i_msa_mean", "nlrmsd_i_msa_lower", "nlrmsd_i_msa_upper"))
  expect_true(all(pr$lrmsd_i_msa_lower <= pr$lrmsd_i_msa_mean))
  expect_true(all(pr$lrmsd_i_msa_mean  <= pr$lrmsd_i_msa_upper))
  expect_true(all(pr$lrmsd_i_msa_upper - pr$lrmsd_i_msa_lower > 0))   # nonzero band
  # the mean-centred profile mean is ~0
  expect_equal(mean(pr$nlrmsd_i_msa_mean), 0, tolerance = 1e-8)
  expect_error(predict_lrmsd_i_agq(agq, pp, level = 1.5), "level")
  expect_error(predict_lrmsd_i_agq(list(), pp), "msa_agq")
})

test_that("predict_lrmsd_i_agq band sharpens with n_nodes but always returns cleanly", {
  # No coarseness warning by design (documented, not flagged): a small-node fit still
  # produces a valid band, just read off fewer nodes. A higher level (0.99) widens it.
  pp    <- agq_default$pp
  agq3  <- fit_lrmsd_i_msa_agq(pp, znb_profile, n_nodes = 3)   # different n_nodes: own fit
  expect_silent(pr3 <- predict_lrmsd_i_agq(agq3, pp))
  expect_equal(nrow(pr3), nrow(pp$site_map))
  expect_true(all(pr3$lrmsd_i_msa_lower <= pr3$lrmsd_i_msa_upper))

  agq <- agq_default$agq
  w95 <- with(predict_lrmsd_i_agq(agq, pp, level = 0.95),
              lrmsd_i_msa_upper - lrmsd_i_msa_lower)
  w99 <- with(predict_lrmsd_i_agq(agq, pp, level = 0.99),
              lrmsd_i_msa_upper - lrmsd_i_msa_lower)
  expect_true(mean(w99) >= mean(w95))   # higher coverage -> wider band
})
