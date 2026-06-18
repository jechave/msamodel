# Profile-invariance guard for the v0.3b dr2* rename.
#
# These snapshots are recorded from the CURRENT (pre-rename) working code and pin
# the numeric profiles the rename must NOT change. After the rename, only the
# function/column SYMBOL references below get updated; the recorded snapshot VALUES
# in _snaps/profile-invariance/ must still reproduce. A failure post-rename is a bug
# in the rename, not drift -- find the cause, do NOT snapshot_accept.
#
# style = "serialize" is used deliberately: it is exactly lossless for doubles
# (json2 is NOT -- ~3e-15 error).

test_that("dr2 site profile is unchanged by the rename", {
  pp <- preprocess_spm(znb_spm)
  prof <- calculate_dr2i_msa(pp, a1 = 2, a2 = 5)$dr2_i
  expect_snapshot_value(prof, style = "serialize")
})

test_that("dr2 mode profile is unchanged by the rename", {
  pp_mode <- preprocess_spm_mode(znb_spm)
  prof_n <- calculate_dr2n_msa(pp_mode, a1 = 2, a2 = 5)$dr2_n
  expect_snapshot_value(prof_n, style = "serialize")
})

test_that("loglik is unchanged by the rename", {
  pp <- preprocess_spm(znb_spm)
  expect_snapshot_value(
    calculate_loglik_msa(pp, znb_profile, a1 = 2, a2 = 5),
    style = "serialize"
  )
})
