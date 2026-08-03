# Leaf-API verbs (calculate_profiles / predict_profiles / calculate_decomposition /
# predict_decomposition). Scope is the ASSEMBLY these verbs add: the list(site, mode)
# shape, the axis keys, the column names per `which`, the `_se` presence on predict,
# and the branching (which validation, the to-be-developed error, the fail-loud fit
# check). The numbers were proved equal to the old exported functions once in the
# scratchpad oracle (disposable), so they are NOT re-encoded here.

test_that("calculate_profiles returns both axes with the which-selected column", {
  spm <- znb_spm
  for (which in c("lrmsd", "nlrmsd")) {
    out <- calculate_profiles(spm, a1 = 1, a2 = 1, which = which)
    expect_named(out, c("site", "mode"))
    expect_named(out$site, c("i", "pdb_site", paste0(which, "_i_msa")))
    expect_named(out$mode, c("n", paste0(which, "_n_msa")))
    # point values only: no _se on calculate_*
    expect_false(any(grepl("_se$", names(out$site))))
  }
})

test_that("predict_profiles adds an _se sibling to every value column, both axes", {
  spm <- znb_spm
  ml  <- fit_lrmsd_i_msa_ml(spm, znb_profile)
  for (which in c("lrmsd", "nlrmsd")) {
    out <- predict_profiles(ml, spm, which = which)
    expect_named(out$site, c("i", "pdb_site", paste0(which, "_i_msa"), paste0(which, "_i_msa_se")))
    expect_named(out$mode, c("n", paste0(which, "_n_msa"), paste0(which, "_n_msa_se")))
    # SEs are real non-negative widths, not NA / all-zero placeholder
    se <- out$site[[paste0(which, "_i_msa_se")]]
    expect_true(all(is.finite(se)) && all(se >= 0) && any(se > 0))
  }
})

test_that("calculate_decomposition switches nested + component family in lockstep", {
  spm <- znb_spm
  lr <- calculate_decomposition(spm, 1, 1, which = "lrmsd")
  nl <- calculate_decomposition(spm, 1, 1, which = "nlrmsd")
  expect_named(lr$site, c("i", "pdb_site", "lrmsd_i_mm", "lrmsd_i_ms", "lrmsd_i_ma",
                          "lrmsd_i_msa", "phi_mut", "phi_stab", "phi_act"))
  expect_named(nl$site, c("i", "pdb_site", "nlrmsd_i_mm", "nlrmsd_i_ms", "nlrmsd_i_ma",
                          "nlrmsd_i_msa", "nphi_mut", "nphi_stab", "nphi_act"))
  # mode branch drops pdb_site (modes are not residues)
  expect_named(nl$mode, c("n", "nlrmsd_n_mm", "nlrmsd_n_ms", "nlrmsd_n_ma",
                          "nlrmsd_n_msa", "nphi_mut", "nphi_stab", "nphi_act"))
})

test_that("predict_decomposition(nlrmsd) bands every nested + component column, both axes", {
  spm <- znb_spm
  ml  <- fit_lrmsd_i_msa_ml(spm, znb_profile)
  mln <- fit_lrmsd_n_msa_ml(spm, znb_profile_n)
  out  <- predict_decomposition(ml, spm, which = "nlrmsd")
  # every value column has an adjacent _se; 7 value cols -> 7 _se cols on each axis
  val_site <- c("nlrmsd_i_mm","nlrmsd_i_ms","nlrmsd_i_ma","nlrmsd_i_msa",
                "nphi_mut","nphi_stab","nphi_act")
  expect_named(out$site, c("i", "pdb_site",
                           as.vector(rbind(val_site, paste0(val_site, "_se")))))
  val_mode <- c("nlrmsd_n_mm","nlrmsd_n_ms","nlrmsd_n_ma","nlrmsd_n_msa",
                "nphi_mut","nphi_stab","nphi_act")
  expect_named(out$mode, c("n",
                           as.vector(rbind(val_mode, paste0(val_mode, "_se")))))
})

test_that("predict_decomposition(which='lrmsd') errors as to-be-developed", {
  spm <- znb_spm
  ml  <- fit_lrmsd_i_msa_ml(spm, znb_profile)
  expect_error(predict_decomposition(ml, spm, which = "lrmsd"), "to be developed")
  # and it must NOT be reachable via a partial-match or a silent default
  expect_error(predict_decomposition(ml, spm, which = "lr"))     # bad which -> match.arg stops
})

test_that("the verbs validate their inputs (which + fit contract)", {
  spm <- znb_spm
  expect_error(calculate_profiles(spm, 1, 1, which = "bogus"))            # match.arg
  expect_error(calculate_decomposition(spm, 1, 1, which = "bogus"))       # match.arg
  # predict_* fail loud on an object that is not an ML fit (missing a1/a2/cov)
  expect_error(predict_profiles(list(a1 = 1, a2 = 1), spm), "2x2 cov")
  expect_error(predict_decomposition(list(a1 = 1, a2 = 1), spm, "nlrmsd"), "2x2 cov")
})
