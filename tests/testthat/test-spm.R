# Contract of the assembled `spm` object (the output of generate_spm, cached as
# znb_spm). The reshaping happens inside that function; what users hold and every
# calculator consumes is this object, so the contract is pinned here on znb_spm directly.

test_that("znb_spm is a classed spm with the five expected elements", {
  expect_s3_class(znb_spm, "spm")
  expect_named(znb_spm, c("energy_data", "dr2mat_site", "dr2mat_mode", "site_map", "mode_map"))
})

test_that("dr2mat_site / dr2mat_mode have the expected shape and carry no column names", {
  # 228 sites x 10 mutations (m > 0) = 2280 mutant rows; 228 site cols, 678 mode cols.
  expect_equal(dim(znb_spm$dr2mat_site), c(2280L, 228L))
  expect_equal(dim(znb_spm$dr2mat_mode), c(2280L, 678L))
  # The response index is the column POSITION, carried explicitly (site_map for sites)
  # -- NOT as matrix colnames, which would leak onto every colSums-derived value vector.
  # Guard that the names stay absent (negative control: the old code set site colnames).
  expect_null(colnames(znb_spm$dr2mat_site))
  expect_null(colnames(znb_spm$dr2mat_mode))
})

test_that("site_map maps the internal site index to pdb_site as in the scan", {
  expect_named(znb_spm$site_map, c("site", "pdb_site"))
  expect_equal(nrow(znb_spm$site_map), 228L)
  # The mapping must match the raw scan's own site / pdb_site vectors -- the whole
  # pdb_site contract depends on this being correct.
  scan <- generate_spm_core(znb_wt, n_mutations = 1,
                            pdb_site_active = PDB_SITE_ACTIVE,
                            ensemble = 1L)
  expect_equal(znb_spm$site_map$site, as.integer(scan$site[[1]]))
  expect_equal(znb_spm$site_map$pdb_site, as.integer(scan$pdb_site[[1]]))
})

test_that("mode_map carries the mode index as in the scan", {
  expect_named(znb_spm$mode_map, "mode")
  expect_equal(nrow(znb_spm$mode_map), 678L)
  # Same shape of check as site_map above: the stored map must match the raw scan's own
  # mode vector, which is where generate_spm takes it from. Modes are not
  # residue-anchored, so there is no pdb_ counterpart -- the index is the whole map.
  scan <- generate_spm_core(znb_wt, n_mutations = 1,
                            pdb_site_active = PDB_SITE_ACTIVE,
                            ensemble = 1L)
  expect_equal(znb_spm$mode_map$mode, as.integer(scan$mode[[1]]))
  # The map must span the matrix it keys -- the invariant prepend_mode_key relies on.
  expect_equal(nrow(znb_spm$mode_map), ncol(znb_spm$dr2mat_mode))
})

test_that("energy_data has the expected columns and drops wild-type rows", {
  expect_named(znb_spm$energy_data, c("j", "m", "ddg", "ddgact"))
  # m = 0 (wild type) is dropped by the reshape.
  expect_false(any(znb_spm$energy_data$m == 0L))
})
