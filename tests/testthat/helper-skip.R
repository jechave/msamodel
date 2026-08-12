# Gate for the heavy FULL SPM regeneration (~21 s: re-runs the entire 2280-row
# generate_spm + set_enm to prove the embedded znb_* fixtures reproduce from
# their recipe). The default suite runs a CHEAP single-mutant coherence check
# instead (see test-spm-generate.R).
#
# There is NO CI in this repo. What runs the full guard is .githooks/pre-commit
# gate 3, automatically, on commits that stage SPM/ENM generation code -- plus a
# by-hand run at a milestone.
# Full run:  MSAMODEL_FULL_TESTS=true devtools::test()
# Default:   devtools::test()  (skips only the full-regen blocks)
skip_if_not_full <- function() {
  testthat::skip_if_not(
    identical(Sys.getenv("MSAMODEL_FULL_TESTS"), "true"),
    "full SPM regeneration (set MSAMODEL_FULL_TESTS=true to run)"
  )
}
