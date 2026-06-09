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

test_that("load_protein reads the embedded PDB as a bio3d object", {
  pdb <- load_protein("1znb_A", data_dir = system.file("extdata", package = "msamodel"))
  expect_s3_class(pdb, "pdb")
})

test_that("get_active_site parses the expected 1znb_A active-site residues", {
  as <- get_active_site("1znb_A", znb_dataset)
  expect_setequal(as$pdb_site_active, c(99, 101, 103, 162, 181, 184, 193, 223))
  expect_setequal(as$site_active, c(78, 80, 82, 141, 160, 163, 172, 202))
})

test_that("setup_enm reproduces the embedded wild-type ENM", {
  wt <- setup_enm(znb_pdb, node = ENM_NODE, model = ENM_MODEL,
                  d_max = ENM_DMAX, frustrated = ENM_FRUST)
  expect_equal(wt, znb_wt)
})

test_that("generate_spm_data reproduces the embedded SPM fixture", {
  # The core drift guard: SPM-generation code vs the embedded znb_spm.
  pdb_site_active <- get_active_site("1znb_A", znb_dataset)$pdb_site_active
  spm <- generate_spm_data(
    znb_wt,
    n_mutations     = SPM_N_MUTATIONS,
    model           = SPM_MODEL,
    sigma           = SPM_SIGMA,
    min_sd          = SPM_MIN_SD,
    pdb_site_active = pdb_site_active,
    seed            = SPM_SEED
  )
  expect_equal(spm, znb_spm, tolerance = 1e-8)
})

test_that("add_site_properties returns the expected site-property columns", {
  sd0 <- tibble::tibble(i = preprocess_spm(znb_spm)$site_map$i)
  pdb_site_active <- get_active_site("1znb_A", znb_dataset)$pdb_site_active
  props <- add_site_properties(sd0, znb_wt, pdb_site_active)
  expect_contains(names(props), c("i", "dactive", "cn", "msf", "shell"))
})
