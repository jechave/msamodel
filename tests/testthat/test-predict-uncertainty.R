# The `uncertainty` toggle on the ten predict_*_ml band functions, and the SPM-sampling
# arm it adds. These are the CHEAP, deterministic wiring/regression checks that run every
# test(). Scope (each can fail for a real reason; see the negative controls):
#   - match.arg validation of the mode argument (a real stop()),
#   - "none" is exactly zero-width, "parameter" reproduces the PRE-REFACTOR band exactly
#     (a frozen regression lock -- the refactor must not have moved the parameter arm),
#   - MM's centred band flips from zero (parameter) to nonzero (spm/both) -- the
#     qualitative fix, with the parameter-mode zero-width as its own control,
#   - nphi_mut band == nlrmsd_mm band exactly (they ARE the same quantity), guarded by
#     phi_stab band != nlrmsd_ms band (a difference is not a single model).
#
# The delta FORMULA's correctness (delta vs an independent naive bootstrap, every family,
# both axes) is validated in dev/spm_band_validation.Rmd -> dev/preview/*.pdf. That is a
# stochastic, many-replicate computation -- deliberately NOT a suite test: it would slow
# every test()/check() with no regression it could catch that these cheap checks don't
# (the formula is fixed; wiring is what regresses, and wiring is covered here).
#
# NOT re-encoded (would be true-by-construction given the code does one `sqrt(vp + vs)`):
# "combined width^2 == spm width^2 + parameter width^2" -- a change-detector, not a
# correctness check, dropped per the test-review gate.

test_that("uncertainty must be one of the four modes (match.arg)", {
  pp <- preprocess_spm(znb_spm)
  ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)
  expect_error(predict_lrmsd_i_msa_ml(ml, pp, uncertainty = "spmm"))   # typo
  expect_error(predict_lrmsd_i_msa_ml(ml, pp, uncertainty = "all"))
  for (u in c("both", "spm", "parameter", "none")) {
    expect_no_error(predict_lrmsd_i_msa_ml(ml, pp, uncertainty = u))
  }
})

test_that('uncertainty = "none" gives an exactly zero-width band, stable schema', {
  pp <- preprocess_spm(znb_spm)
  ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)
  pn   <- predict_lrmsd_i_msa_ml(ml, pp, uncertainty = "none")
  pboth <- predict_lrmsd_i_msa_ml(ml, pp, uncertainty = "both")

  expect_equal(pn$lrmsd_i_msa_lower, pn$lrmsd_i_msa_mean)
  expect_equal(pn$lrmsd_i_msa_upper, pn$lrmsd_i_msa_mean)
  # schema identical across modes (not fewer columns for "none")
  expect_identical(names(pn), names(pboth))
  # negative control: "both" is NOT zero-width on the same column, so the equalities
  # above are a real property of "none", not a vacuous 0==0.
  expect_gt(max(pboth$lrmsd_i_msa_upper - pboth$lrmsd_i_msa_mean), 1e-6)
})

test_that('uncertainty = "parameter" reproduces the pre-refactor band exactly', {
  # Frozen regression lock. Before this work item the predictors had NO uncertainty arg
  # and their band was parameter-only. The values below were captured from that code
  # (git HEAD at the start of the item) and are hard-coded literals -- NOT recomputed by
  # the code under test -- so a change in the parameter arm makes this test go red.
  # (First six residues of the znb site fit; full vectors match to 0 as verified in the
  # scratchpad.)
  pp <- preprocess_spm(znb_spm)
  ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)
  p  <- predict_lrmsd_i_msa_ml(ml, pp, uncertainty = "parameter")
  pn <- predict_nlrmsd_i_msa_ml(ml, pp, uncertainty = "parameter")

  # frozen expected (pre-refactor lrmsd parameter band, first 6 residues). Literals are
  # 10-sig-fig rounded, so the default tolerance (~1.5e-8) is what the lock checks.
  expect_equal(head(p$lrmsd_i_msa_lower, 6),  PREREF$lrmsd_lower)
  expect_equal(head(p$lrmsd_i_msa_upper, 6),  PREREF$lrmsd_upper)
  expect_equal(head(pn$nlrmsd_i_msa_lower, 6), PREREF$nlrmsd_lower)
  expect_equal(head(pn$nlrmsd_i_msa_upper, 6), PREREF$nlrmsd_upper)
})

test_that("MM centred band: zero under parameter, nonzero under spm/both", {
  # The qualitative fix. MM has zero parameter gradient -> zero PARAMETER band, but a
  # nonzero SPM-sampling band. Both facts asserted on the SAME column, so neither a
  # "MM always zero" nor a "MM always nonzero" bug can pass.
  pp <- preprocess_spm(znb_spm)
  ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)

  par <- predict_nlrmsd_i_nested_models_ml(ml, pp, uncertainty = "parameter")
  spm <- predict_nlrmsd_i_nested_models_ml(ml, pp, uncertainty = "spm")
  bth <- predict_nlrmsd_i_nested_models_ml(ml, pp, uncertainty = "both")

  expect_equal(par$nlrmsd_i_mm_upper, par$nlrmsd_i_mm_lower)               # parameter: zero
  expect_gt(max(spm$nlrmsd_i_mm_upper - spm$nlrmsd_i_mm_lower), 1e-6)      # spm: nonzero
  expect_gt(max(bth$nlrmsd_i_mm_upper - bth$nlrmsd_i_mm_lower), 1e-6)      # both: nonzero
})

test_that("nphi_mut band == nlrmsd_mm band (same quantity); phi_stab differs", {
  # phi_mut = mm by construction (calculate_msa_decomposition sets phi_mut = mm), so its
  # SPM band must equal the nlrmsd_mm band exactly. Negative control: phi_stab is a
  # DIFFERENCE (ms - mm), so its band must NOT equal the single-model nlrmsd_ms band --
  # if the decomposition wrongly reused a single-model formula, that control goes red.
  pp <- preprocess_spm(znb_spm)
  ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)
  dec <- predict_nlrmsd_i_msa_decomposition_ml(ml, pp, uncertainty = "spm")
  nn  <- predict_nlrmsd_i_nested_models_ml(ml, pp, uncertainty = "spm")

  w_mut  <- dec$nphi_mut_upper  - dec$nphi_mut_mean
  w_stab <- dec$nphi_stab_upper - dec$nphi_stab_mean
  w_mm   <- nn$nlrmsd_i_mm_upper - nn$nlrmsd_i_mm_mean
  w_ms   <- nn$nlrmsd_i_ms_upper - nn$nlrmsd_i_ms_mean

  expect_equal(w_mut, w_mm)                       # same quantity, identical band
  expect_gt(max(abs(w_stab - w_ms)), 1e-6)        # difference != single model
})
