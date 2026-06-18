# Fixture-drift guard: the SPM-generation code must keep reproducing the embedded
# znb_* fixtures. The parameters below MUST match data-raw/prepare_znb_data.R.
ENM_NODE  <- "ca"
ENM_MODEL <- "ming_wall"
ENM_DMAX  <- 10.5
ENM_FRUST <- FALSE
SPM_N_MUTATIONS <- 10
SPM_MODEL       <- "lfenm"
SPM_SIGMA       <- 0.3
SPM_MIN_SD      <- 2
SPM_SEED        <- 1024

# 1znb_A active-site residues (PDB numbering), from data-raw/raw/dataset_1znb_A.csv.
# Must match data-raw/prepare_znb_data.R.
PDB_SITE_ACTIVE <- c(99, 101, 103, 162, 181, 184, 193, 223)

test_that("bio3d::read.pdb reads the embedded PDB as a bio3d object", {
  pdb <- bio3d::read.pdb(system.file("extdata", "1znb_A.pdb", package = "msamodel"))
  expect_s3_class(pdb, "pdb")
})

test_that("setup_enm reproduces the embedded wild-type ENM", {
  wt <- setup_enm(znb_pdb, node = ENM_NODE, model = ENM_MODEL,
                  d_max = ENM_DMAX, frustrated = ENM_FRUST)
  expect_equal(wt, znb_wt)
})

test_that("generate_spm_data reproduces the embedded SPM fixture", {
  # The core drift guard: SPM-generation code vs the embedded znb_spm.
  spm <- generate_spm_data(
    znb_wt,
    n_mutations     = SPM_N_MUTATIONS,
    model           = SPM_MODEL,
    sigma           = SPM_SIGMA,
    min_sd          = SPM_MIN_SD,
    pdb_site_active = PDB_SITE_ACTIVE,
    seed            = SPM_SEED
  )
  expect_equal(spm, znb_spm, tolerance = 1e-8)
})

test_that("SPM mode column equals penm's mode index (guards the seq_along choice)", {
  # generate_spm_data stores `mode` as seq_along(dr2n) for speed; this confirms
  # that vector is exactly penm:::get_mode(wt), i.e. the mode-axis labels are
  # the real normal-mode indices, not an off-by-one stand-in.
  expect_equal(znb_spm$mode[[1]], penm:::get_mode(znb_wt))
})

test_that("add_site_properties returns the expected site-property columns", {
  sd0 <- tibble::tibble(i = preprocess_spm(znb_spm)$site_map$i)
  props <- add_site_properties(sd0, znb_wt, PDB_SITE_ACTIVE)
  expect_contains(names(props), c("i", "dactive", "cn", "msf", "shell"))
})
