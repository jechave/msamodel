# Per-element drift guard for the structural-divergence profile vectors.
#
# These snapshots pin the FULL per-element profiles -- the dr2_i (site) and dr2_n
# (mode) vectors at a fixed (a1, a2), plus the loglik scalar -- against any future
# change to the dr2_msa forward-map primitive / calculate_loglik_lrmsd_i_msa.
# This is the ONLY value-level guard on those whole vectors: the frozen loglik
# literal elsewhere is a scalar REDUCTION (a bug that permutes sites or flips signs
# symmetric about the mean leaves it unchanged while the profile drifts), and the
# nested-model tests recompute via the function under test (circular for drift). The
# snapshots catch what those cannot.
#
# (Originally written to guard the v0.3b dr2* rename; that shipped, but the standing
# job -- pinning the per-element profiles -- remains, so the snapshots stay.)
#
# A failure here is a real change in the profile values: find the cause, do NOT
# snapshot_accept unless the change is intended and understood. The recorded values
# live in _snaps/profile-invariance.md.
#
# style = "serialize" is used deliberately: it is exactly lossless for doubles
# (json2 is NOT -- ~3e-15 error).

test_that("dr2 site profile is unchanged by the rename", {
  pp <- znb_spm
  prof <- dr2_msa(pp$dr2_ijm, pp$energy_data, a1 = 2, a2 = 5)
  expect_snapshot_value(prof, style = "serialize")
})

test_that("dr2 mode profile is unchanged by the rename", {
  pp_mode <- znb_spm
  prof_n <- dr2_msa(pp_mode$dr2_njm, pp_mode$energy_data, a1 = 2, a2 = 5)
  expect_snapshot_value(prof_n, style = "serialize")
})

test_that("loglik is unchanged by the rename", {
  pp <- znb_spm
  expect_snapshot_value(
    msamodel:::calculate_loglik_lrmsd_i_msa(pp, znb_profile$pdb_site, znb_profile$lrmsd_i_obs, a1 = 2, a2 = 5),
    style = "serialize"
  )
})
