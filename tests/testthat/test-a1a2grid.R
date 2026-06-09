test_that("define_selection_grid returns the default 6x9 grid", {
  g <- define_selection_grid()
  expect_named(g, c("a1", "a2"))
  expect_equal(nrow(g), 54L)            # n_a1 = 6, n_a2 = 9
  # grid is sorted; the first combination is the null model (0, 0)
  expect_equal(g$a1[1], 0)
  expect_equal(g$a2[1], 0)
})

test_that("calculate_dr2i_msa_a1a2grid output has the expected columns and alignment", {
  grid <- define_selection_grid(n_a1 = 2, n_a2 = 2)   # 4 combinations
  res <- calculate_dr2i_msa_a1a2grid(znb_spm, grid, verbose = FALSE)

  expect_contains(names(res),
                  c("a1", "a2", "i", "pdb_site",
                    "dr2_msa", "dr2_mm", "dr2_ms", "dr2_ma"))
  # one row per (grid combination x site)
  expect_equal(nrow(res), nrow(grid) * 228L)

  # pdb_site must be aligned to i via the model's site_map
  site_map <- preprocess_spm(znb_spm)$site_map
  expected_pdb <- site_map$pdb_site[match(res$i, site_map$i)]
  expect_equal(res$pdb_site, expected_pdb)
})
