# Mode-form structural divergence: the mode matrix dr2mat_mode on the assembled spm
# object + the @noRd forward-map primitives (dr2_msa / lrmsd_nested_models /
# nlrmsd_msa / nlrmsd_msa_decomposition / nlrmsd_nested_models) that the leaf verbs
# are built from. Mode mirror of test-msa-evaluate.R. The old exported calculate_*_n_*
# grid wrappers were deleted; the math they wrapped is these primitives.

test_that("znb_spm carries a [mutant x mode] dr2mat_mode of the expected shape", {
  # 228 sites x 10 mutations (m > 0) = 2280 mutant rows; 678 modes (= 3*228 - 6).
  expect_equal(dim(znb_spm$dr2mat_mode), c(2280L, 678L))
  # Mode index is the column position, not stored as colnames (which would leak onto
  # colSums-derived value vectors); the forward map recovers n via seq_len(ncol).
  expect_null(colnames(znb_spm$dr2mat_mode))
})

test_that("dr2_msa returns one finite, positive dr2 per mode", {
  d <- dr2_msa(znb_spm$dr2mat_mode, znb_spm$energy_data, a1 = 2, a2 = 5)
  expect_length(d, 678L)
  expect_true(all(is.finite(d)))
  expect_true(all(d > 0))
})

test_that("per-mutant site and mode divergence agree (basis invariance)", {
  # Total structural divergence is basis-invariant: dr2n = (U^T dr)^2 with U
  # orthonormal, so sum over modes == sum over sites for every mutant. This ties
  # the two SPM matrices together via a real physical invariant -- it is NOT a
  # recomputation of the reshape/evaluate code, and it would fail if dr2mat_site/dr2mat_mode
  # were ever computed from different displacements or mis-stored (e.g. transposed).
  rows <- c(1L, 500L, 2000L)
  for (r in rows) {
    expect_equal(sum(znb_spm$dr2mat_site[r, ]), sum(znb_spm$dr2mat_mode[r, ]), tolerance = 1e-8)
  }
})

test_that("lrmsd_nested_models builds the four variants at the right (a1,a2) (mode)", {
  a1 <- 2; a2 <- 5
  nested <- lrmsd_nested_models(znb_spm$dr2mat_mode, znb_spm$energy_data, a1, a2)

  expect_named(nested, c("mm", "ms", "ma", "msa"))
  # Each variant = log(sqrt(dr2)) of dr2_msa at its (a1,a2) point.
  # Independent route (recompute the forward map directly), not circular.
  expect_equal(nested$mm,  log(sqrt(dr2_msa(znb_spm$dr2mat_mode, znb_spm$energy_data, 0,  0 ))))
  expect_equal(nested$ms,  log(sqrt(dr2_msa(znb_spm$dr2mat_mode, znb_spm$energy_data, a1, 0 ))))
  expect_equal(nested$ma,  log(sqrt(dr2_msa(znb_spm$dr2mat_mode, znb_spm$energy_data, 0,  a2))))
  expect_equal(nested$msa, log(sqrt(dr2_msa(znb_spm$dr2mat_mode, znb_spm$energy_data, a1, a2))))
})

test_that("nlrmsd_msa centres the uncentred profile and agrees with its predictor (mode)", {
  a1 <- 2; a2 <- 5
  nc <- nlrmsd_msa(znb_spm$dr2mat_mode, znb_spm$energy_data, a1, a2)

  # Independent route: centred == uncentred lrmsd minus its full-support mean.
  # Negative control: an uncentred return would fail (profile mean far from zero).
  lr <- lrmsd_msa(znb_spm$dr2mat_mode, znb_spm$energy_data, a1, a2)
  expect_equal(nc, lr - mean(lr))

  # Forward-map primitive == predictor point profile.
  mln  <- fit_lrmsd_msa_mode(znb_spm, znb_profile_n$mode, znb_profile_n$lrmsd_obs)
  fwd  <- nlrmsd_msa(znb_spm$dr2mat_mode, znb_spm$energy_data, mln$a1, mln$a2)
  pred <- predict_profiles(mln, znb_spm, metric = "nlrmsd")$mode
  expect_equal(fwd, pred$nlrmsd_msa)
})

test_that("nlrmsd_msa_decomposition contributions sum to the centred profile (mode)", {
  a1 <- 2; a2 <- 5
  d <- nlrmsd_msa_decomposition(znb_spm$dr2mat_mode, znb_spm$energy_data, a1, a2)

  expect_named(d, c("nphi_mut", "nphi_stab", "nphi_act"))
  # Independent route: nlrmsd_msa from its own primitive. Negative control: dropping the
  # centring on any nphi term breaks the equality.
  prof <- nlrmsd_msa(znb_spm$dr2mat_mode, znb_spm$energy_data, a1, a2)
  expect_equal(d$nphi_mut + d$nphi_stab + d$nphi_act, prof)
})

test_that("dr2_n reweighting collapses the mutant axis with the same weights as the site form", {
  # Independent-route check: build the same fixation weights by hand and apply
  # colSums(dr2mat_mode * w). Confirms dr2_msa uses the axis-agnostic mutant-axis weights,
  # not a mode-specific scheme. Energies come straight from the raw scan
  # (generate_spm_core), NOT from the assembled object under test, so this is not
  # circular with the forward map.
  a1 <- 2; a2 <- 5

  # The regenerated scan supplies the ENERGIES, which are then applied to the
  # FIXTURE's dr2mat_mode -- so the two must be the same realization. SPM_ENSEMBLE is
  # sourced from the recipe rather than copied, so it cannot drift out of step.
  source(testthat::test_path("fixtures", "make-znb-fixtures.R"))
  scan <- generate_spm_core(znb_wt, n_mutations = 10, model = "lfenm", sigma = 0.3,
                            min_sd = 2,
                            pdb_site_active = PDB_SITE_ACTIVE,
                            ensemble = SPM_ENSEMBLE)
  filt <- scan[scan$m > 0, ]
  ddg    <- filt$ddg_dv_jm + filt$ddg_tds_jm
  ddgact <- filt$ddgact_dv_jm + filt$ddgact_tds_jm
  pfix <- pmin(exp(-a1 * ddg), 1) * pmin(exp(-a2 * ddgact), 1)
  w <- pfix / sum(pfix)
  expected <- colSums(znb_spm$dr2mat_mode * w)

  d <- dr2_msa(znb_spm$dr2mat_mode, znb_spm$energy_data, a1, a2)
  expect_equal(d, expected)
})
