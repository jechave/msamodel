# Helper: a small synthetic prediction-samples tibble (two samples, three sites).
new_pred_samples <- function(with_pdb_site = FALSE) {
  d <- tibble::tibble(
    sample_id = rep(1:2, each = 3),
    i         = rep(1:3, times = 2),
    lrmsd_mm  = c(0.1, 0.2, 0.3, 0.15, 0.25, 0.35),
    lrmsd_ms  = c(0.2, 0.3, 0.4, 0.25, 0.35, 0.45),
    lrmsd_ma  = c(0.3, 0.4, 0.5, 0.35, 0.45, 0.55),
    lrmsd_msa = c(0.5, 0.6, 0.7, 0.55, 0.65, 0.75)
  )
  if (with_pdb_site) d$pdb_site <- rep(c(20L, 21L, 22L), times = 2)
  d
}

test_that("decomposition functions validate required columns", {
  expect_error(calculate_msa_decomposition(tibble::tibble(i = 1)),
               "Missing required columns")
  expect_error(calculate_decomposition_samples(tibble::tibble(sample_id = 1, i = 1)),
               "Missing required columns")
  expect_error(calculate_decomposition_summary(tibble::tibble(sample_id = 1, i = 1)),
               "Missing required columns")
})

test_that("decomposition works without pdb_site (optional-branch, absent)", {
  ds <- calculate_decomposition_samples(new_pred_samples(with_pdb_site = FALSE))
  expect_false("pdb_site" %in% names(ds))
  expect_contains(names(ds), c("sample_id", "i", "shap_mut", "shap_stab", "shap_act"))

  dsum <- calculate_decomposition_summary(ds)
  expect_false("pdb_site" %in% names(dsum))
  expect_named(dsum, c("i", "component", "mean", "sd", "median", "lower", "upper"))
  # one row per (site x component): 3 sites x 3 components
  expect_equal(nrow(dsum), 9L)
})

test_that("decomposition carries pdb_site through (optional-branch, present)", {
  ds <- calculate_decomposition_samples(new_pred_samples(with_pdb_site = TRUE))
  expect_contains(names(ds), "pdb_site")

  dsum <- calculate_decomposition_summary(ds)
  expect_contains(names(dsum), "pdb_site")
  expect_equal(nrow(dsum), 9L)
})
