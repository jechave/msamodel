test_that("calculate_msa_decomposition is a pure vector function (4 in, 3 out)", {
  # Plain numeric vectors -- no tibble, no column names, no key. This is the point
  # of the refactor: the math knows nothing about sites/modes.
  mm  <- c(0.1, 0.2, 0.3)
  ms  <- c(0.2, 0.3, 0.4)
  ma  <- c(0.3, 0.4, 0.5)   # unused by sequential; reserved for Shapley
  msa <- c(0.5, 0.6, 0.7)
  out <- calculate_msa_decomposition(mm, ms, ma, msa)

  expect_type(out, "list")
  expect_named(out, c("phi_mut", "phi_stab", "phi_act"))

  # Sequential formula: phi_mut = mm, phi_stab = ms - mm, phi_act = msa - ms.
  # Differs from symmetric Shapley (phi_stab = 0.5*(ms-mm + msa-ma) = .15 here);
  # these assertions FAIL under Shapley, so they pin the right form.
  expect_equal(out$phi_mut,  mm)
  expect_equal(out$phi_stab, ms - mm)
  expect_equal(out$phi_act,  msa - ms)
  expect_equal(out$phi_stab[1], 0.1)   # Shapley would give 0.15
  expect_equal(out$phi_act[1], 0.3)    # Shapley would give 0.25

  # Telescoping: terms sum to the full predicted divergence; ma does not affect it.
  expect_equal(out$phi_mut + out$phi_stab + out$phi_act, msa)
})

test_that("calculate_decomposition_i_msa packages nested-models + the kernel", {
  pp <- znb_spm
  got <- calculate_decomposition_i_msa(pp, a1 = 1, a2 = 1)

  # Shape / keys: one row per site, site keys carried, phi columns present.
  expect_named(got, c("i", "pdb_site", "phi_mut", "phi_stab", "phi_act"))
  expect_equal(got$i, sort(got$i))   # ordered by site

  # It equals the manual two-step it packages (its contract).
  nm  <- calculate_lrmsd_i_nested_models(pp, a1 = 1, a2 = 1)
  ref <- calculate_msa_decomposition(nm$lrmsd_i_mm, nm$lrmsd_i_ms,
                                     nm$lrmsd_i_ma, nm$lrmsd_i_msa)
  expect_equal(got$phi_mut,  ref$phi_mut)
  expect_equal(got$phi_stab, ref$phi_stab)
  expect_equal(got$phi_act,  ref$phi_act)

  # Independent check: the three contributions sum to the full-model profile.
  # (a1=1, a2=1 keeps all three non-trivial -- not a degenerate zero point.)
  expect_equal(got$phi_mut + got$phi_stab + got$phi_act, nm$lrmsd_i_msa)
  expect_true(all(got$phi_stab != 0) && all(got$phi_act != 0))
})

test_that("calculate_decomposition_n_msa packages nested-models + the kernel (mode)", {
  pp <- znb_spm
  got <- calculate_decomposition_n_msa(pp, a1 = 1, a2 = 1)

  # Shape / keys: one row per mode, keyed on n, no pdb_site (modes aren't residues).
  expect_named(got, c("n", "phi_mut", "phi_stab", "phi_act"))
  expect_false("pdb_site" %in% names(got))
  expect_equal(got$n, sort(got$n))

  # It equals the manual two-step it packages (its contract).
  nm  <- calculate_lrmsd_n_nested_models(pp, a1 = 1, a2 = 1)
  ref <- calculate_msa_decomposition(nm$lrmsd_n_mm, nm$lrmsd_n_ms,
                                     nm$lrmsd_n_ma, nm$lrmsd_n_msa)
  expect_equal(got$phi_mut,  ref$phi_mut)
  expect_equal(got$phi_stab, ref$phi_stab)
  expect_equal(got$phi_act,  ref$phi_act)

  # Independent check: the three contributions sum to the full-model profile.
  expect_equal(got$phi_mut + got$phi_stab + got$phi_act, nm$lrmsd_n_msa)
})
