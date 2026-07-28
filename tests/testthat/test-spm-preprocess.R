test_that("preprocess_spm returns energy_data, dr2_ijm, and site_map", {
  pp <- preprocess_spm(znb_spm)
  expect_named(pp, c("energy_data", "dr2_ijm", "site_map"))
})

test_that("dr2_ijm has the expected shape and carries no column names", {
  pp <- preprocess_spm(znb_spm)
  # 228 sites x 10 mutations (m > 0 rows) = 2280 rows; 228 site columns.
  expect_equal(dim(pp$dr2_ijm), c(2280L, 228L))
  # The site index is the column POSITION, carried explicitly in site_map -- NOT as
  # matrix colnames, which would leak onto every colSums-derived value vector. Guard
  # that the names stay absent (negative control: the old code set them from site).
  expect_null(colnames(pp$dr2_ijm))
})

test_that("site_map maps the internal index i to pdb_site as in the SPM", {
  pp <- preprocess_spm(znb_spm)
  expect_named(pp$site_map, c("i", "pdb_site"))
  expect_equal(nrow(pp$site_map), 228L)
  # The mapping must match the SPM's own site / pdb_site list-columns -- the whole
  # pdb_site contract depends on this being correct.
  expect_equal(pp$site_map$i, as.integer(znb_spm$site[[1]]))
  expect_equal(pp$site_map$pdb_site, as.integer(znb_spm$pdb_site[[1]]))
})

test_that("energy_data has the expected columns", {
  pp <- preprocess_spm(znb_spm)
  expect_named(pp$energy_data, c("j", "m", "ddg_jm", "ddgact_jm"))
})
