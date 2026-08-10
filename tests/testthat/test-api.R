# Leaf-API verbs (calculate_profiles / predict_profiles / calculate_decomposition /
# predict_decomposition). Scope is the ASSEMBLY these verbs add: the list(site, mode)
# shape, the axis keys, the column names per `metric`, the `_se` presence on predict,
# and the branching (metric validation, the to-be-developed error, the fail-loud fit
# check). The numbers were proved equal to the old exported functions once in the
# scratchpad oracle (disposable), so they are NOT re-encoded here.

test_that("calculate_profiles returns both axes with the metric-selected column", {
  spm <- znb_spm
  for (metric in c("lrmsd", "nlrmsd")) {
    out <- calculate_profiles(spm, a1 = 1, a2 = 1, metric = metric)
    expect_named(out, c("site", "mode"))
    expect_named(out$site, c("site", "pdb_site", paste0(metric, "_msa")))
    expect_named(out$mode, c("mode", paste0(metric, "_msa")))
    # point values only: no _se on calculate_*
    expect_false(any(grepl("_se$", names(out$site))))
  }
})

test_that("predict_profiles adds an _se sibling to every value column, both axes", {
  spm <- znb_spm
  ml  <- fit_lrmsd_msa_site(spm, znb_profile$pdb_site, znb_profile$lrmsd_obs)
  for (metric in c("lrmsd", "nlrmsd")) {
    out <- predict_profiles(ml, spm, metric = metric)
    # Value columns first, then the matching _se columns (one of each here).
    expect_named(out$site, c("site", "pdb_site", paste0(metric, "_msa"), paste0(metric, "_msa_se")))
    expect_named(out$mode, c("mode", paste0(metric, "_msa"), paste0(metric, "_msa_se")))
    # SEs are real non-negative widths, not NA / all-zero placeholder
    se <- out$site[[paste0(metric, "_msa_se")]]
    expect_true(all(is.finite(se)) && all(se >= 0) && any(se > 0))
  }
})

test_that("calculate_decomposition switches nested + component columns together", {
  spm <- znb_spm
  lr <- calculate_decomposition(spm, 1, 1, metric = "lrmsd")
  nl <- calculate_decomposition(spm, 1, 1, metric = "nlrmsd")
  expect_named(lr$site, c("site", "pdb_site", "lrmsd_mm", "lrmsd_ms", "lrmsd_ma",
                          "lrmsd_msa", "phi_mut", "phi_stab", "phi_act"))
  expect_named(nl$site, c("site", "pdb_site", "nlrmsd_mm", "nlrmsd_ms", "nlrmsd_ma",
                          "nlrmsd_msa", "nphi_mut", "nphi_stab", "nphi_act"))
  # mode branch drops pdb_site (modes are not residues)
  expect_named(nl$mode, c("mode", "nlrmsd_mm", "nlrmsd_ms", "nlrmsd_ma",
                          "nlrmsd_msa", "nphi_mut", "nphi_stab", "nphi_act"))
  # the two branches speak ONE value-column vocabulary: they differ only in the key
  expect_identical(setdiff(names(nl$site), c("site", "pdb_site")),
                   setdiff(names(nl$mode), "mode"))
})

test_that("predict_decomposition(nlrmsd) gives every nested + component column an _se, both axes", {
  spm <- znb_spm
  ml  <- fit_lrmsd_msa_site(spm, znb_profile$pdb_site, znb_profile$lrmsd_obs)
  mln <- fit_lrmsd_msa_mode(spm, znb_profile_n$mode, znb_profile_n$lrmsd_obs)
  out  <- predict_decomposition(ml, spm, metric = "nlrmsd")
  # The 7 value columns come FIRST, then the 7 _se columns in the same order -- they are
  # grouped, not interleaved. `expect_named` is order-sensitive, so this pins the layout.
  # Both branches share one value vocabulary, so one `vals` drives both assertions.
  vals <- c("nlrmsd_mm","nlrmsd_ms","nlrmsd_ma","nlrmsd_msa",
            "nphi_mut","nphi_stab","nphi_act")
  with_se <- c(vals, paste0(vals, "_se"))
  expect_named(out$site, c("site", "pdb_site", with_se))
  expect_named(out$mode, c("mode", with_se))
  # Negative control on the layout itself: the OLD interleaved order must now fail.
  expect_false(identical(names(out$site),
                         c("site", "pdb_site", as.vector(rbind(vals, paste0(vals, "_se"))))))
})

test_that("predict_decomposition(metric='lrmsd') errors as to-be-developed", {
  spm <- znb_spm
  ml  <- fit_lrmsd_msa_site(spm, znb_profile$pdb_site, znb_profile$lrmsd_obs)
  # The error now comes from the se_*_lrmsd stubs in R/predict_se.R, not from a guard at
  # the top of the verb: the uncentred branch is undeveloped, and the stub is where it
  # will be developed.
  expect_error(predict_decomposition(ml, spm, metric = "lrmsd"), "to be developed")
  expect_error(predict_decomposition(ml, spm, metric = "lr"), "to be developed")  # partial match
  expect_error(predict_decomposition(ml, spm, metric = "bogus"))                  # match.arg stops
  # The DEFAULT is nlrmsd, so a bare call reaches the working branch instead of erroring
  # on its own default (which is what it used to do).
  expect_identical(predict_decomposition(ml, spm),
                   predict_decomposition(ml, spm, metric = "nlrmsd"))
})

test_that("the verbs validate their inputs (metric + fit contract)", {
  spm <- znb_spm
  expect_error(calculate_profiles(spm, 1, 1, metric = "bogus"))            # match.arg
  expect_error(calculate_decomposition(spm, 1, 1, metric = "bogus"))       # match.arg
  # predict_* fail loud on an object that is not an ML fit (missing a1/a2/cov)
  expect_error(predict_profiles(list(a1 = 1, a2 = 1), spm), "2x2 cov")
  expect_error(predict_decomposition(list(a1 = 1, a2 = 1), spm, "nlrmsd"), "2x2 cov")
})

test_that("a metric with no branch fails loud instead of falling through to lrmsd", {
  # The verbs dispatch on `metric` with an explicit branch per metric and a closing
  # stop(), NOT an `else` that silently means lrmsd. match.arg cannot catch the failure
  # this guards: a metric ADDED to a verb's formals (drmsf, dnh, ...) whose branch was
  # never written. match.arg would accept it and the old `else` would have returned
  # lrmsd numbers under the new metric's name.
  #
  # Reaching the guard therefore requires simulating that future state -- rewriting
  # formals() to admit an unimplemented metric -- because no public call can get there.
  spm <- znb_spm
  ml  <- fit_lrmsd_msa_site(spm, znb_profile$pdb_site, znb_profile$lrmsd_obs)

  for (fname in c("calculate_profiles", "calculate_decomposition",
                  "predict_profiles", "predict_decomposition")) {
    f <- get(fname)
    formals(f)$metric <- c("drmsf", "lrmsd", "nlrmsd")   # accepted by match.arg, unimplemented
    args <- if (startsWith(fname, "predict")) list(ml, spm) else list(spm, 1, 1)
    expect_error(do.call(f, c(args, list(metric = "drmsf"))),
                 'metric "drmsf" is not available',
                 info = fname)
  }

  # Control: the same call shape on an IMPLEMENTED metric must not error, so the guard
  # is specific to the unimplemented one rather than blanket-rejecting the rewrite.
  f <- calculate_profiles
  formals(f)$metric <- c("drmsf", "lrmsd", "nlrmsd")
  expect_no_error(do.call(f, list(spm, 1, 1, metric = "nlrmsd")))
})

test_that("a site_map inconsistent with the profile length fails loud, not silently NA", {
  # prepend_site_key binds site_map on POSITIONALLY (it is already the (site, pdb_site)
  # key table in dr2_ijm column order). The previous left_join against a manufactured
  # index would have silently filled pdb_site with NA on a mismatch; this must error.
  bad <- znb_spm
  bad$site_map <- bad$site_map[1:100, ]
  expect_error(calculate_profiles(bad, 1, 1), "site_map has 100 rows")
  # Control: the untampered object still works, so the guard is not blanket.
  expect_no_error(calculate_profiles(znb_spm, 1, 1))
})

test_that("a mode_map inconsistent with the profile length fails loud, not silently NA", {
  # The mode axis used to MANUFACTURE its key (seq_len(n)), so nothing could be
  # inconsistent and nothing was checked. It now binds the STORED spm$mode_map
  # positionally, exactly as the site axis binds site_map -- so it needs the same guard.
  spm <- znb_spm
  ml  <- fit_lrmsd_msa_site(spm, znb_profile$pdb_site, znb_profile$lrmsd_obs)
  bad <- spm
  bad$mode_map <- bad$mode_map[1:100, , drop = FALSE]
  expect_error(calculate_profiles(bad, 1, 1), "mode_map has 100 rows")
  expect_error(predict_profiles(ml, bad), "mode_map has 100 rows")
  # Control: the untampered object still works, so the guard is not blanket.
  expect_no_error(predict_profiles(ml, spm))
})

# ---- Frozen numeric reference ----------------------------------------------------
# The standalone value guard for the four verbs. Until the old exported grid functions
# are deleted, their equality with these verbs is checked by a disposable scratchpad
# oracle (refactor-invariance, one-time). THIS block is the permanent reference that
# survives that deletion: it pins the verbs' output to literals frozen from znb_spm at a
# fixed (a1, a2) and a fixed deterministic ML fit, so a future refactor that silently
# corrupts a value goes red here. Expected values are hand-frozen literals, NOT recomputed
# by the code under test. Regenerate ONLY on an intended, reviewed output change.

test_that("calculate_* reproduce frozen reference values (znb_spm at a1=1.3, a2=0.7)", {
  spm <- znb_spm
  tol <- 1e-8

  cp_l <- calculate_profiles(spm, 1.3, 0.7, "lrmsd")
  cp_n <- calculate_profiles(spm, 1.3, 0.7, "nlrmsd")
  expect_equal(cp_l$site$lrmsd_msa[1:3],
               c(-2.38860098918393, -3.41729557242227, -3.40211980354922), tolerance = tol)
  expect_equal(cp_l$mode$lrmsd_msa[1:3],
               c(-1.5687305271568, -1.82327102574257, -2.21838411202063), tolerance = tol)
  expect_equal(cp_n$site$nlrmsd_msa[1:3],
               c(1.81659196313855, 0.787897379900211, 0.803073148773266), tolerance = tol)
  expect_equal(cp_n$mode$nlrmsd_msa[1:3],
               c(4.02007678194882, 3.76553628336305, 3.37042319708499), tolerance = tol)
  expect_equal(cp_l$site$pdb_site[1:3], c(20L, 21L, 22L))
  expect_equal(nrow(cp_l$site), 228L)
  expect_equal(nrow(cp_l$mode), 678L)

  cd_l <- calculate_decomposition(spm, 1.3, 0.7, "lrmsd")
  cd_n <- calculate_decomposition(spm, 1.3, 0.7, "nlrmsd")
  expect_equal(cd_l$site$lrmsd_mm[1], -2.42154799236766, tolerance = tol)
  expect_equal(cd_l$site$phi_stab[1],    0.0287338535954338, tolerance = tol)
  expect_equal(cd_l$site$phi_act[1],     0.00421314958829955, tolerance = tol)
  expect_equal(cd_n$site$nphi_stab[1],   0.379280436055831, tolerance = tol)
  expect_equal(cd_n$site$nphi_act[1],    0.00933343562163675, tolerance = tol)
  expect_equal(cd_n$mode$nphi_stab[1],   0.534774823663534, tolerance = tol)
  expect_equal(cd_n$mode$nphi_act[1],   -0.00628475618290865, tolerance = tol)
})

test_that("predict_* reproduce frozen reference values (fixed ML fit; _se pinned)", {
  spm  <- znb_spm
  ml_i <- fit_lrmsd_msa_site(spm, znb_profile$pdb_site, znb_profile$lrmsd_obs)     # deterministic; a1>0 and a2>0
  ml_n <- fit_lrmsd_msa_mode(spm, znb_profile_n$mode, znb_profile_n$lrmsd_obs)
  tol  <- 1e-6                                      # se's carry Hessian round-off; looser
  expect_true(ml_i$a1 > 0 && ml_i$a2 > 0)          # both arms live (guards the arm-point mix)

  pp_i <- predict_profiles(ml_i, spm, "nlrmsd")
  pp_n <- predict_profiles(ml_n, spm, "nlrmsd")
  expect_equal(pp_i$site$nlrmsd_msa[1],    1.81683984987468, tolerance = tol)
  expect_equal(pp_i$site$nlrmsd_msa_se[1], 0.110443916713325, tolerance = tol)
  expect_equal(pp_n$mode$nlrmsd_msa[1],    3.2741642946694,  tolerance = tol)
  expect_equal(pp_n$mode$nlrmsd_msa_se[1], 0.0826541993570476, tolerance = tol)

  pd_i <- predict_decomposition(ml_i, spm, "nlrmsd")
  pd_n <- predict_decomposition(ml_n, spm, "nlrmsd")
  # nested-model _se (the load-bearing columns: a mis-wired arm/point would hide here)
  expect_equal(pd_i$site$nlrmsd_mm_se[1],  0.102662779822919, tolerance = tol)
  expect_equal(pd_i$site$nlrmsd_ms_se[1],  0.108697825909678, tolerance = tol)
  expect_equal(pd_i$site$nlrmsd_msa_se[1], 0.110443916713325, tolerance = tol)
  # component _se
  expect_equal(pd_i$site$nphi_stab_se[1],    0.0851289549563722, tolerance = tol)
  expect_equal(pd_i$site$nphi_act_se[1],     0.035098062683972,  tolerance = tol)
  expect_equal(pd_n$mode$nlrmsd_mm_se[1],  0.0738708293607008, tolerance = tol)
  expect_equal(pd_n$mode$nphi_stab_se[1],    0.041478920753972,  tolerance = tol)

  # Exact structural invariants (hold regardless of the frozen numbers): the mutation
  # contribution IS the MM nested model, value and se.
  expect_identical(pd_i$site$nphi_mut,    pd_i$site$nlrmsd_mm)
  expect_identical(pd_i$site$nphi_mut_se, pd_i$site$nlrmsd_mm_se)
})
