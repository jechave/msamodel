# Fixture-drift guard: the SPM-generation code must keep reproducing the znb_*
# fixtures in tests/testthat/fixtures/.
#
# The generation constants are SOURCED from the fixture recipe, not copied here.
# There is exactly one declaration, so this guard cannot silently test a different
# recipe than the one that built the fixture -- which is what the old duplicated
# block, carrying a "MUST match" comment and nothing enforcing it, allowed.
# Sourcing is side-effect-free: the recipe only runs when executed as a script.
source(testthat::test_path("fixtures", "make-znb-fixtures.R"))

# PDB_SITE_ACTIVE comes from helper-fixtures.R, which reads the same shared
# inst/extdata/znb_active_site.csv the vignettes and data-raw/ use.

test_that("bio3d::read.pdb reads the embedded PDB as a bio3d object", {
  pdb <- bio3d::read.pdb(system.file("extdata", "1znb_A.pdb", package = "msamodel"))
  expect_s3_class(pdb, "pdb")
})

test_that("a single mutant row reproduces from its seeded recipe (cheap coherence)", {
  # Always-on cache-coherence guard. znb_spm is a COMPUTED RESULT cached in data/
  # for speed; this confirms the committed cache still matches what the generator
  # produces -- WITHOUT the ~21 s full 2280-row regeneration (gated below).
  #
  # The scan draws each mutant via get_mutant_site(wt, site_mut, mutation, ...,
  # ensemble=SPM_ENSEMBLE): penm keys the per-mutant RNG on the hashed tuple
  # (ensemble, site_mut, mutation), so one
  # (j, m>0) cell reproduces exactly without running the whole loop. Verified to be
  # bit-exact (diff = 0), not merely within tolerance. m = 0 is the WILD TYPE and is
  # dropped in the assembled object -- this MUST use m > 0 to exercise the mutation
  # path. If single-row reproduction ever fails, that is a fail-loud drift bug.
  #
  # znb_spm is now the assembled `spm` object: the per-mutant divergences live as ROWS
  # of the dr2mat_site / dr2mat_mode matrices, aligned to energy_data row order. We locate the
  # (j=1, m=3) mutant's row k in energy_data and compare against the matrix rows.
  ed <- znb_spm$energy_data
  k  <- which(ed$j == 1L & ed$m == 3L)              # a fixed, real (j, m>0) mutant
  expect_length(k, 1L)

  mut <- penm::get_mutant_site(znb_wt, ed$j[k], ed$m[k],
                               mut_model = SPM_MODEL, mut_dl_sigma = SPM_SIGMA,
                               mut_sd_min = SPM_MIN_SD, ensemble = SPM_ENSEMBLE)

  # Same measured divergences the generation loop records for this row: the matrix row
  # compares directly against penm's per-mutant vectors (both nameless).
  expect_equal(znb_spm$dr2mat_site[k, ], penm::delta_structure_dr2i(znb_wt, mut))
  expect_equal(znb_spm$dr2mat_mode[k, ], penm::delta_structure_dr2n(znb_wt, mut))
  # Summed stability / activity energy changes (energy_data carries the summed terms).
  expect_equal(ed$ddg[k],    penm::ddg_dv(znb_wt, mut)  + penm::ddg_tds(znb_wt, mut))
  expect_equal(ed$ddgact[k], penm::ddgact_dv(znb_wt, mut, pdb_site_active = PDB_SITE_ACTIVE) +
                                penm::ddgact_tds(znb_wt, mut, pdb_site_active = PDB_SITE_ACTIVE))
})

test_that("set_enm reproduces the wild-type ENM fixture", {
  skip_if_not_full()
  # Reads the shipped PDB file, the same input make-znb-fixtures.R uses -- there is no
  # serialized bio3d object any more, because a user starts from a .pdb file too.
  pdb <- bio3d::read.pdb(system.file("extdata", "1znb_A.pdb", package = "msamodel"))
  wt  <- penm::set_enm(pdb, node = ENM_NODE, model = ENM_MODEL,
                         d_max = ENM_DMAX, frustrated = ENM_FRUST)
  expect_equal(wt, znb_wt)
})

test_that("generate_spm reproduces the embedded SPM fixture", {
  skip_if_not_full()
  # The core drift guard: SPM-generation code vs the embedded znb_spm. Heavy (~21 s
  # full 2280-row regen), so it is gated behind MSAMODEL_FULL_TESTS and does NOT run
  # in a default devtools::test().
  #
  # What actually runs it (there is no CI in this repo -- do not assume one):
  #   - .githooks/pre-commit gate 3, automatically, on any commit staging
  #     R/spm.R, the fixture recipe, or data-raw/prepare_znb_data.R;
  #   - by hand at a milestone:  MSAMODEL_FULL_TESTS=true devtools::test()
  #
  # The always-on single-mutant check above covers cache drift on every default run,
  # but it verifies ONE row of 2280 and locates that row by lookup -- so a change
  # that reorders rows or alters which are included passes it. This is the check
  # that sees that class of change.
  spm <- generate_spm(
    znb_wt,
    n_mutations     = SPM_N_MUTATIONS,
    model           = SPM_MODEL,
    sigma           = SPM_SIGMA,
    min_sd          = SPM_MIN_SD,
    pdb_site_active = PDB_SITE_ACTIVE,
    ensemble        = SPM_ENSEMBLE
  )
  expect_equal(spm, znb_spm, tolerance = 1e-8)
})

test_that("the scan's mode labels are the dense column positions of the mode axis", {
  # The scan takes `mode` from penm (penm::get_mode()), so comparing it back against
  # get_mode() would be a tautology. The invariant that survives, and that the rest of
  # the mode axis depends on, is the LABEL-vs-WIDTH agreement: the mode vector must be
  # the dense 1-based run that spans the mode axis it labels.
  #
  # The two sides come from different penm code paths -- `nma$mode` via get_mode(), and
  # the return length of delta_structure_dr2n() -- so this is not a recomputation.
  #
  # It goes red if a future penm ships a mode set that is truncated or not 1-based: the
  # whole mode axis (mode_map, and the dr2mat_mode column-position convention it keys)
  # assumes column n is mode n, so such a change would silently mislabel every
  # mode-axis result rather than erroring. Verified to fail on a 7-based mode index.
  scan <- generate_spm_core(znb_wt, n_mutations = 1, pdb_site_active = PDB_SITE_ACTIVE,
                            ensemble = SPM_ENSEMBLE)
  expect_length(scan$mode[[1]], length(scan$dr2_njm[[1]]))
  expect_equal(scan$mode[[1]], seq_len(length(scan$dr2_njm[[1]])))
})
