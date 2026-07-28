# Mode-form structural divergence (v0.3a): preprocess_spm_mode + calculate_dr2_n_msa.
# Predict-only path (no fit), parallel to the site path in test-msa-evaluate.R.

test_that("preprocess_spm_mode returns energy_data and a [mutant x mode] dr2_njm", {
  pp <- preprocess_spm_mode(znb_spm)
  expect_named(pp, c("energy_data", "dr2_njm"))
  expect_named(pp$energy_data, c("j", "m", "ddg_jm", "ddgact_jm"))
  # 228 sites x 10 mutations (m > 0 rows) = 2280 rows; 678 modes (= 3*228 - 6).
  expect_equal(dim(pp$dr2_njm), c(2280L, 678L))
  # Mode index is the column position, not stored as colnames (which would leak onto
  # colSums-derived value vectors); calculate_dr2_n_msa recovers n via seq_len(ncol).
  expect_null(colnames(pp$dr2_njm))
})

test_that("calculate_dr2_n_msa returns one finite, positive dr2_n per mode", {
  pp <- preprocess_spm_mode(znb_spm)
  d <- calculate_dr2_n_msa(pp, a1 = 2, a2 = 5)
  expect_named(d, c("n", "dr2_n"))
  expect_equal(nrow(d), 678L)
  expect_equal(d$n, 1:678)
  expect_true(all(is.finite(d$dr2_n)))
  expect_true(all(d$dr2_n > 0))
})

test_that("per-mutant site and mode divergence agree (basis invariance)", {
  # Total structural divergence is basis-invariant: dr2n = (U^T dr)^2 with U
  # orthonormal, so sum over modes == sum over sites for every mutant. This ties
  # the two SPM columns together via a real physical invariant -- it is NOT a
  # recomputation of preprocess/evaluate (the code under test), and it would fail
  # if dr2/dr2n were ever computed from different displacements or mis-stored.
  rows <- which(znb_spm$m > 0)[c(1, 500, 2000)]
  for (r in rows) {
    expect_equal(sum(znb_spm$dr2_ijm[[r]]), sum(znb_spm$dr2_njm[[r]]), tolerance = 1e-8)
  }
})

test_that("calculate_lrmsd_n_nested_models builds the four variants at the right (a1,a2)", {
  pp <- preprocess_spm_mode(znb_spm)
  a1 <- 2; a2 <- 5
  nested <- calculate_lrmsd_n_nested_models(pp, a1, a2)

  # Mode form: keyed by n, NO pdb_site (modes are not anchored to residues).
  expect_named(nested, c("n", "lrmsd_n_mm", "lrmsd_n_ms", "lrmsd_n_ma", "lrmsd_n_msa"))
  expect_equal(nrow(nested), 678L)

  # Each variant = log(sqrt(dr2_n)) of calculate_dr2_n_msa at its (a1,a2) point.
  # Independent route (recompute the forward map directly), not circular.
  expect_equal(nested$lrmsd_n_mm,  log(sqrt(calculate_dr2_n_msa(pp, 0,  0 )$dr2_n)))
  expect_equal(nested$lrmsd_n_ms,  log(sqrt(calculate_dr2_n_msa(pp, a1, 0 )$dr2_n)))
  expect_equal(nested$lrmsd_n_ma,  log(sqrt(calculate_dr2_n_msa(pp, 0,  a2)$dr2_n)))
  expect_equal(nested$lrmsd_n_msa, log(sqrt(calculate_dr2_n_msa(pp, a1, a2)$dr2_n)))
})

test_that("calculate_nlrmsd_n_msa centres the uncentred profile and agrees with its predictor", {
  pp <- preprocess_spm_mode(znb_spm)
  a1 <- 2; a2 <- 5
  nc <- calculate_nlrmsd_n_msa(pp, a1, a2)

  expect_named(nc, c("n", "nlrmsd_n_msa"))
  expect_equal(nrow(nc), 678L)

  # Independent route: centred == uncentred lrmsd minus its full-support mean.
  # Negative control: an uncentred return would fail (profile mean far from zero).
  lr <- calculate_lrmsd_n_msa(pp, a1, a2)$lrmsd_n_msa
  expect_equal(nc$nlrmsd_n_msa, lr - mean(lr))

  # Forward map == predictor point profile (uncertainty='none' strips the band).
  mln  <- fit_lrmsd_n_msa_ml(pp, znb_profile_n)
  fwd  <- calculate_nlrmsd_n_msa(pp, mln$a1, mln$a2)
  pred <- predict_nlrmsd_n_msa_ml(mln, pp, uncertainty = "none")
  expect_equal(fwd$nlrmsd_n_msa, pred$nlrmsd_n_msa_mean)
})

test_that("calculate_nlrmsd_n_msa_decomposition contributions sum to the centred profile", {
  pp <- preprocess_spm_mode(znb_spm)
  a1 <- 2; a2 <- 5
  d <- calculate_nlrmsd_n_msa_decomposition(pp, a1, a2)

  expect_named(d, c("n", "nphi_mut", "nphi_stab", "nphi_act"))
  # Independent route: nlrmsd_n_msa from its own twin. Negative control: dropping the
  # centring on any nphi term breaks the equality.
  prof <- calculate_nlrmsd_n_msa(pp, a1, a2)$nlrmsd_n_msa
  expect_equal(d$nphi_mut + d$nphi_stab + d$nphi_act, prof)
})

test_that("dr2_n reweighting collapses the mutant axis with the same weights as the site form", {
  # Independent-route check: build the same fixation weights by hand and apply
  # colSums(dr2_njm * w). Confirms calculate_dr2_n_msa uses the axis-agnostic
  # mutant-axis weights, not a mode-specific scheme. Energies come straight from
  # the SPM, NOT from preprocess output, so this is not circular with Step 2.
  a1 <- 2; a2 <- 5
  pp <- preprocess_spm_mode(znb_spm)

  filt <- znb_spm[znb_spm$m > 0, ]
  ddg_jm    <- filt$ddg_dv_jm + filt$ddg_tds_jm
  ddgact_jm <- filt$ddgact_dv_jm + filt$ddgact_tds_jm
  pfix <- pmin(exp(-a1 * ddg_jm), 1) * pmin(exp(-a2 * ddgact_jm), 1)
  w <- pfix / sum(pfix)
  expected <- colSums(pp$dr2_njm * w)

  d <- calculate_dr2_n_msa(pp, a1, a2)
  expect_equal(unname(d$dr2_n), unname(expected))
})
