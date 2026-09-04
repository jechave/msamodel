# Fail-loud paths.
#
# Package policy is that an unexpected input or state ERRORS with a clear message rather
# than being papered over -- no silent NA, no quiet fallback, no clamping. A silent NA is
# the runtime version of an unverified reassurance: it turns "something is wrong here"
# into "looks fine".
#
# These tests encode that policy. They are cheap, they exercise no science, and each one
# fails if a future edit replaces a stop() with a tolerant default.

test_that("a predictor rejects an object that is not an ML fit", {
  # Guards against a caller passing, say, the spm or a bare list of parameters and
  # getting a plausible-looking band computed from a missing covariance.
  expect_error(predict_profiles(list(a1 = 1, a2 = 1), spm),
               "must be an ML fit")
  expect_error(predict_decomposition(list(a1 = 1, a2 = 1), spm),
               "must be an ML fit")
})

test_that("an inconsistent spm is rejected rather than silently recycled", {
  # site_map is bound on POSITIONALLY to the profile body. If the two disagree, binding
  # anyway would mislabel every row -- dplyr would either recycle or error obscurely.
  bad <- spm
  bad$site_map <- bad$site_map[1:50, ]
  expect_error(calculate_profiles(bad, 1, 1), "site_map has 50 rows")

  bad_mode <- spm
  bad_mode$mode_map <- bad_mode$mode_map[1:50, , drop = FALSE]
  expect_error(calculate_profiles(bad_mode, 1, 1), "mode_map has 50 rows")
})

test_that("unimplemented lrmsd standard errors error instead of returning a wrong band", {
  # The lrmsd-metric SEs are genuinely not implemented. The risk this guards is a future
  # edit quietly returning the nlrmsd band for an lrmsd request -- a wrong error bar is
  # worse than a refusal, because it looks usable.
  ml <- fit_lrmsd_msa_site(spm, obs_site$pdb_site, obs_site$lrmsd_obs)
  expect_error(predict_decomposition(ml, spm, metric = "lrmsd"),
               "still to be developed")
})

test_that("the fitter validates its search ranges", {
  expect_error(
    fit_lrmsd_msa_site(spm, obs_site$pdb_site, obs_site$lrmsd_obs, a1_range = c(5, 1)),
    "a1_range must be"
  )
  expect_error(
    fit_lrmsd_msa_site(spm, obs_site$pdb_site, obs_site$lrmsd_obs,
                       log2_a2_plus1_range = c(9, 2)),
    "log2_a2_plus1_range must be"
  )
})

test_that("spm_hmat rejects misaligned axes and unnormalised weights", {
  # The SPM-sampling arm is a weighted mean over mutants. Unnormalised weights would
  # scale every standard error by an arbitrary factor while still producing finite,
  # plausible numbers -- the exact class of silent corruption worth an assertion.
  dr2 <- spm$dr2mat_site
  w   <- msamodel:::weights_jm(spm$energy_data, 1, 1)

  expect_error(msamodel:::spm_hmat(dr2, w[1:10], centred = FALSE),
               "mutant axes are misaligned")
  expect_error(msamodel:::spm_hmat(dr2, w * 2, centred = FALSE),
               "must sum to 1")
})
