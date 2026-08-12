test_that("decompose_nested is a pure vector function (4 in, 3 out)", {
  # Plain numeric vectors -- no tibble, no column names, no key. This is the point
  # of the refactor: the math knows nothing about sites/modes.
  mm  <- c(0.1, 0.2, 0.3)
  ms  <- c(0.2, 0.3, 0.4)
  ma  <- c(0.3, 0.4, 0.5)   # unused by sequential; reserved for Shapley
  msa <- c(0.5, 0.6, 0.7)
  out <- decompose_nested(mm, ms, ma, msa)

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

test_that("lrmsd_msa_decomposition packages nested-models + the kernel (site)", {
  pp <- znb_spm
  got <- lrmsd_msa_decomposition(pp$dr2mat_site, pp$energy_data, a1 = 1, a2 = 1)

  expect_named(got, c("phi_mut", "phi_stab", "phi_act"))

  # It equals the manual two-step it packages (its contract): nested models, then kernel.
  nm  <- lrmsd_nested_models(pp$dr2mat_site, pp$energy_data, a1 = 1, a2 = 1)
  ref <- decompose_nested(nm$mm, nm$ms, nm$ma, nm$msa)
  expect_equal(got$phi_mut,  ref$phi_mut)
  expect_equal(got$phi_stab, ref$phi_stab)
  expect_equal(got$phi_act,  ref$phi_act)

  # Independent check: the three contributions sum to the full-model profile.
  # (a1=1, a2=1 keeps all three non-trivial -- not a degenerate zero point.)
  expect_equal(got$phi_mut + got$phi_stab + got$phi_act, nm$msa)
  expect_true(all(got$phi_stab != 0) && all(got$phi_act != 0))
})

test_that("lrmsd_msa_decomposition packages nested-models + the kernel (mode)", {
  pp <- znb_spm
  got <- lrmsd_msa_decomposition(pp$dr2mat_mode, pp$energy_data, a1 = 1, a2 = 1)

  expect_named(got, c("phi_mut", "phi_stab", "phi_act"))

  # It equals the manual two-step it packages (its contract).
  nm  <- lrmsd_nested_models(pp$dr2mat_mode, pp$energy_data, a1 = 1, a2 = 1)
  ref <- decompose_nested(nm$mm, nm$ms, nm$ma, nm$msa)
  expect_equal(got$phi_mut,  ref$phi_mut)
  expect_equal(got$phi_stab, ref$phi_stab)
  expect_equal(got$phi_act,  ref$phi_act)

  # Independent check: the three contributions sum to the full-model profile.
  expect_equal(got$phi_mut + got$phi_stab + got$phi_act, nm$msa)
})
