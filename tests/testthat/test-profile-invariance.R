# Per-element drift guard for the structural-divergence profile vectors.
#
# These snapshots pin the FULL per-element profiles -- the dr2_i (site) and dr2_n
# (mode) vectors at a fixed (a1, a2), plus the loglik scalar -- against any future
# change to the dr2_msa forward-map primitive / loglik_lrmsd_msa.
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
    msamodel:::loglik_lrmsd_msa(pp$dr2_ijm, pp$energy_data, msamodel:::resolve_site_obs(pp, znb_profile$pdb_site, znb_profile$lrmsd_obs), a1 = 2, a2 = 5),
    style = "serialize"
  )
})

# ---- Public API verbs: whole-output pins ------------------------------------
# The four leaf verbs were pinned only by FIRST-ELEMENT literals in test-api.R
# ([1:3] for profiles, [1] for decomposition and every _se). Those stay -- they are
# a lock a person must retype -- but they cannot see a change at site 150, and the
# argument at the top of this file applies to the verbs as much as to the
# primitives: a permutation, or a sign flip symmetric about the mean, leaves a
# spot check green while the profile drifts.
#
# predict_* is the reason this matters most: its delta-method standard errors are
# real math the primitives never perform, and nothing else re-derives them.
#
# The whole returned object is snapshotted (keys, value columns, _se columns), so a
# column silently added, dropped, reordered, or renamed also goes red here.
#
# (a1 = 1.3, a2 = 0.7) and the ML fits below match test-api.R's frozen-reference
# block, so a failure there and here points at the same change.
#
# DO NOT EDIT THESE test_that DESCRIPTIONS. A snapshot is keyed by its description:
# renaming one orphans the recorded value, and testthat then writes a NEW snapshot
# with a warning instead of failing -- the comparison silently stops happening.

test_that("calculate_profiles output is pinned (site and mode, both metrics)", {
  spm <- znb_spm
  expect_snapshot_value(calculate_profiles(spm, 1.3, 0.7, "lrmsd"),  style = "serialize")
  expect_snapshot_value(calculate_profiles(spm, 1.3, 0.7, "nlrmsd"), style = "serialize")
})

test_that("calculate_decomposition output is pinned (site and mode, both metrics)", {
  spm <- znb_spm
  expect_snapshot_value(calculate_decomposition(spm, 1.3, 0.7, "lrmsd"),  style = "serialize")
  expect_snapshot_value(calculate_decomposition(spm, 1.3, 0.7, "nlrmsd"), style = "serialize")
})

test_that("predict_profiles output including standard errors is pinned", {
  spm  <- znb_spm
  ml_i <- fit_lrmsd_msa_site(spm, znb_profile$pdb_site, znb_profile$lrmsd_obs)
  ml_n <- fit_lrmsd_msa_mode(spm, znb_profile_n$mode, znb_profile_n$lrmsd_obs)
  expect_snapshot_value(predict_profiles(ml_i, spm, "nlrmsd"), style = "serialize")
  expect_snapshot_value(predict_profiles(ml_n, spm, "nlrmsd"), style = "serialize")
})

test_that("predict_decomposition output including standard errors is pinned", {
  spm  <- znb_spm
  ml_i <- fit_lrmsd_msa_site(spm, znb_profile$pdb_site, znb_profile$lrmsd_obs)
  ml_n <- fit_lrmsd_msa_mode(spm, znb_profile_n$mode, znb_profile_n$lrmsd_obs)
  expect_snapshot_value(predict_decomposition(ml_i, spm, "nlrmsd"), style = "serialize")
  expect_snapshot_value(predict_decomposition(ml_n, spm, "nlrmsd"), style = "serialize")
})
