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

test_that("a single mutant cell reproduces from its seeded recipe (cheap coherence)", {
  # Always-on cache-coherence guard. znb_spm is a COMPUTED RESULT cached in data/
  # for speed; this confirms the committed cache still matches what the generator
  # produces -- WITHOUT the ~21 s full 2280-row regeneration (gated below).
  #
  # generate_spm_data seeds per mutant via get_mutant_site(wt, site_mut, mutation,
  # ..., seed=SPM_SEED): the per-mutant RNG is keyed on (seed, site_mut, mutation),
  # so one (j, m>0) cell reproduces exactly without running the whole loop. Verified
  # to be bit-exact (diff = 0), not merely within tolerance. m = 0 would be the WILD
  # TYPE (no mutation) and is preprocessed away -- this MUST use m > 0 to exercise
  # the mutation path. If single-cell reproduction ever fails, that is a fail-loud
  # drift bug to investigate, not to relax.
  idx <- which(znb_spm$j == 1L & znb_spm$m == 3L)   # a fixed, real (j, m>0) mutant
  expect_length(idx, 1L)

  mut <- get_mutant_site(znb_wt, znb_spm$j[idx], znb_spm$m[idx],
                         mut_model = SPM_MODEL, mut_dl_sigma = SPM_SIGMA,
                         mut_sd_min = SPM_MIN_SD, seed = SPM_SEED)

  # Same measured effects the generation loop records for this cell.
  expect_equal(penm::delta_structure_dr2i(znb_wt, mut), znb_spm$dr2_ijm[[idx]])
  expect_equal(penm::delta_structure_dr2n(znb_wt, mut), znb_spm$dr2_njm[[idx]])
  expect_equal(ddg_dv(znb_wt, mut),  znb_spm$ddg_dv_jm[idx])
  expect_equal(ddg_tds(znb_wt, mut), znb_spm$ddg_tds_jm[idx])
  expect_equal(ddgact_dv(znb_wt, mut, pdb_site_active = PDB_SITE_ACTIVE),
               znb_spm$ddgact_dv_jm[idx])
  expect_equal(ddgact_tds(znb_wt, mut, pdb_site_active = PDB_SITE_ACTIVE),
               znb_spm$ddgact_tds_jm[idx])
})

test_that("setup_enm reproduces the embedded wild-type ENM", {
  skip_if_not_full()
  wt <- setup_enm(znb_pdb, node = ENM_NODE, model = ENM_MODEL,
                  d_max = ENM_DMAX, frustrated = ENM_FRUST)
  expect_equal(wt, znb_wt)
})

test_that("generate_spm_data reproduces the embedded SPM fixture", {
  skip_if_not_full()
  # The core drift guard: SPM-generation code vs the embedded znb_spm. Heavy (~21 s
  # full 2280-row regen); gated to milestone / CI via MSAMODEL_FULL_TESTS. The
  # always-on single-mutant check above covers cache drift on every default run.
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
