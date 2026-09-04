# Drift pins on the public verbs.
#
# WHAT THIS FILE CLAIMS, AND WHAT IT DOES NOT. These snapshots record that the output of
# each public verb has not moved since it was last recorded. That is the entire claim.
# A snapshot CANNOT tell you a value was ever correct -- it freezes whatever the code
# produced, and a wrong number freezes as happily as a right one. Correctness of the MSA
# model is established elsewhere (the paper, and the delta-method validation study that
# was done once and is not re-run here); these pins protect it from silent erosion.
#
# WHY THIS IS THE CORE OF THE SUITE. Every verb below reads spm, which helper-setup.R
# BUILDS rather than loads. So the chain pdb -> ENM -> scan -> forward map -> fit ->
# prediction is exercised end to end on every run, and a regression anywhere in it
# necessarily moves one of these snapshots. That is why the suite needs no separate
# guard on generate_spm(): its output has nowhere to hide.
#
# The whole returned object is snapshotted, not selected values, so a column silently
# added, dropped, reordered or renamed goes red here too.
#
# style = "serialize" is deliberate: exactly lossless for doubles (json2 carries ~3e-15
# error).
#
# ON A FAILURE: do NOT run snapshot_accept() to make it green. A red snapshot means a
# value moved -- find out which and why first. Accept only once the change is understood
# and intended, and say so in the commit.
#
# DO NOT EDIT THESE test_that DESCRIPTIONS. A snapshot is keyed by its description:
# renaming one orphans the recorded value, and testthat then writes a NEW snapshot with
# only a warning instead of failing -- the comparison silently stops happening.

# Fixed evaluation point for the calculate_* verbs. Arbitrary but non-degenerate: both
# pressures are on, so a bug that drops either arm shows up.
A1 <- 1.3
A2 <- 0.7

test_that("generate_spm output is pinned", {
  # The scan itself, on the small 2-mutation fixture. The verb pins below would also
  # move if the generator changed, but they would report it as a calculate_*/predict_*
  # failure; this one names the generator.
  #
  # The whole object is pinned -- both dr2 matrices in full, the energy table and both
  # key maps -- for the same reason as every other pin here: a partial pin leaves the
  # unpinned part free to drift.
  expect_snapshot_value(spm_small, style = "serialize")
})

test_that("calculate_profiles output is pinned (both axes, both metrics)", {
  expect_snapshot_value(calculate_profiles(spm, A1, A2, "lrmsd"),  style = "serialize")
  expect_snapshot_value(calculate_profiles(spm, A1, A2, "nlrmsd"), style = "serialize")
})

test_that("calculate_decomposition output is pinned (both axes, both metrics)", {
  expect_snapshot_value(calculate_decomposition(spm, A1, A2, "lrmsd"),  style = "serialize")
  expect_snapshot_value(calculate_decomposition(spm, A1, A2, "nlrmsd"), style = "serialize")
})

test_that("the site fit is pinned", {
  ml <- fit_lrmsd_msa_site(spm, obs_site$pdb_site, obs_site$lrmsd_obs)
  # $call records match.call() and would churn on any refactor of the wrapper's
  # signature without indicating a numeric change, so it is excluded deliberately.
  expect_snapshot_value(ml[c("a1", "a2", "logLik", "cov", "se_a1", "se_a2", "gof")],
                        style = "serialize")
})

test_that("the mode fit is pinned", {
  # MACHINERY ONLY. obs_mode is SYNTHETIC (see helper-setup.R): it was built by
  # evaluating this very model and adding seeded noise, so these numbers say nothing
  # about real per-mode divergence. What they pin is that the mode-axis code path runs
  # and its arithmetic has not drifted.
  ml <- fit_lrmsd_msa_mode(spm, obs_mode$mode, obs_mode$lrmsd_obs)
  expect_snapshot_value(ml[c("a1", "a2", "logLik", "cov", "se_a1", "se_a2", "gof")],
                        style = "serialize")
})

test_that("predict_profiles output including standard errors is pinned", {
  # The delta-method standard errors are the reason this matters most: they are real
  # math the forward map never performs, and nothing else in the suite re-derives them.
  ml_site <- fit_lrmsd_msa_site(spm, obs_site$pdb_site, obs_site$lrmsd_obs)
  ml_mode <- fit_lrmsd_msa_mode(spm, obs_mode$mode, obs_mode$lrmsd_obs)
  expect_snapshot_value(predict_profiles(ml_site, spm, "nlrmsd"), style = "serialize")
  expect_snapshot_value(predict_profiles(ml_mode, spm, "nlrmsd"), style = "serialize")
})

test_that("predict_decomposition output including standard errors is pinned", {
  ml_site <- fit_lrmsd_msa_site(spm, obs_site$pdb_site, obs_site$lrmsd_obs)
  ml_mode <- fit_lrmsd_msa_mode(spm, obs_mode$mode, obs_mode$lrmsd_obs)
  expect_snapshot_value(predict_decomposition(ml_site, spm, "nlrmsd"), style = "serialize")
  expect_snapshot_value(predict_decomposition(ml_mode, spm, "nlrmsd"), style = "serialize")
})
