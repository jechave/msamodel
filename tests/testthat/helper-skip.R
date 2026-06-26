# Gate for the heavy FULL SPM regeneration (~21 s: re-runs the entire 2280-row
# generate_spm_data + setup_enm to prove the embedded znb_* fixtures reproduce from
# their data-raw recipe). The default suite runs a CHEAP single-mutant coherence
# check instead (see test-spm-generate.R); this full guard runs at milestone / CI.
# Full run:  MSAMODEL_FULL_TESTS=true devtools::test()
# Default:   devtools::test()  (skips only the full-regen blocks)
skip_if_not_full <- function() {
  testthat::skip_if_not(
    identical(Sys.getenv("MSAMODEL_FULL_TESTS"), "true"),
    "full SPM regeneration (set MSAMODEL_FULL_TESTS=true to run)"
  )
}
