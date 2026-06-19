# Helper: a small synthetic prediction-samples tibble (two samples, three sites).
new_pred_samples <- function(with_pdb_site = FALSE) {
  d <- tibble::tibble(
    sample_id   = rep(1:2, each = 3),
    i           = rep(1:3, times = 2),
    lrmsd_i_mm  = c(0.1, 0.2, 0.3, 0.15, 0.25, 0.35),
    lrmsd_i_ms  = c(0.2, 0.3, 0.4, 0.25, 0.35, 0.45),
    lrmsd_i_ma  = c(0.3, 0.4, 0.5, 0.35, 0.45, 0.55),
    lrmsd_i_msa = c(0.5, 0.6, 0.7, 0.55, 0.65, 0.75)
  )
  if (with_pdb_site) d$pdb_site <- rep(c(20L, 21L, 22L), times = 2)
  d
}

test_that("the plumbing functions validate required columns", {
  # calculate_msa_decomposition is a pure vector fn (no column validation); the
  # tibble-based plumbing fns do the validation.
  expect_error(calculate_decomposition_samples(tibble::tibble(sample_id = 1, i = 1)),
               "Missing required columns")
  expect_error(calculate_decomposition_summary(tibble::tibble(sample_id = 1, i = 1)),
               "Missing required columns")
})

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

test_that("decomposition_samples works without pdb_site (optional-branch, absent)", {
  ds <- calculate_decomposition_samples(new_pred_samples(with_pdb_site = FALSE))
  expect_false("pdb_site" %in% names(ds))
  expect_contains(names(ds), c("sample_id", "i", "phi_mut", "phi_stab", "phi_act"))

  # Values match the sequential formula applied per row.
  d <- new_pred_samples(with_pdb_site = FALSE)
  expect_equal(ds$phi_stab, d$lrmsd_i_ms - d$lrmsd_i_mm)
  expect_equal(ds$phi_act,  d$lrmsd_i_msa - d$lrmsd_i_ms)

  dsum <- calculate_decomposition_summary(ds)
  expect_false("pdb_site" %in% names(dsum))
  expect_named(dsum, c("i", "component", "mean", "sd", "median", "lower", "upper"))
  # one row per (site x component): 3 sites x 3 components
  expect_equal(nrow(dsum), 9L)
})

test_that("decomposition_samples carries pdb_site through (optional-branch, present)", {
  ds <- calculate_decomposition_samples(new_pred_samples(with_pdb_site = TRUE))
  expect_contains(names(ds), "pdb_site")

  dsum <- calculate_decomposition_summary(ds)
  expect_contains(names(dsum), "pdb_site")
  expect_equal(nrow(dsum), 9L)
})
