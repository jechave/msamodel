# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this package is

`msamodel` implements the Mutation-Stability-Activity (MSA) model of neutral
structural evolution: it predicts site-specific structural divergence in protein
evolution from single-point-mutation scans (SPMs). Structure/active-site input is a
bio3d pdb object + a plain `pdb_site_active` integer vector.

v0.1 and v0.2 (API cleanup) shipped; the 0.3.0 site+mode structural suite is done
(predict + fit, both axes); the inference rework is in flight. The package computes
structure-divergence profiles (structure × site and × mode) and fits them to data:
SPM generation, divergence profiles, parameter estimation (ML / AGQ / MCMC), and a
site-level decomposition (the `phi_*` columns). **Read the `<!-- NOW -->` block in
`dev/LOG.md` first each session** — it holds the live state; `dev/plan.md` is the
roadmap.

It depends on `penm` (an `Imports:` dependency) for the ENM machinery — `set_enm()`,
mutation scans, response matrices. **penm is a dependency, never a source**: call
`penm::fn()`; never copy or wrap its code, and never rename a penm function to fit a
local convention.

## Planning & cadence (how work flows)

Planning lives in **two files**: `dev/plan.md` (coarse, durable roadmap — what each
version *is*; edited only when a version's *goal* changes) and `dev/LOG.md`
(append-only history + a delimited `<!-- NOW -->` … `<!-- /NOW -->` block at the top
holding the live agenda). Read the NOW block first; keep it current.

Work is reviewed at the **work-item level**, not per commit:

1. **Plan up front.** Enter plan mode, read the code the item touches, sketch the
   approach, and get the user's approval before implementing. Skip the plan only for
   changes you could describe in one sentence (typo, log line, rename).
2. **Execute autonomously.** Run the inner loop (`load_all()` → targeted
   `test_file()`) across as many atomic commits as the item needs. Commits are
   logical units, but they are **not** individual review gates — don't stop the user
   at each one. The test suite is the between-touchpoint verification loop, so keep
   it strong (see the rigor gate below).
3. **Review the diff.** At the item's end, show the user the finished diff plus a
   verification artifact sized to the change (a scratchpad table/PNG for code; a
   `dev/preview/<name>.html` render for a vignette). Nothing is pushed until the user
   has reviewed and said go.

If a change goes wrong mid-way, `git restore` the uncommitted tree and re-approach —
nothing is committed mid-unit.

### The rigor gate — decided by the DIFF, not by how big the change felt

Before committing, whether the full `devtools::test()` runs is decided by one rule:

- **Diff touches `R/`, `data/`, `data-raw/`, roxygen, `NAMESPACE`, or any test
  snapshot value → run ONE full `devtools::test()`.** This is the rigor gate.
- **Otherwise the diff is non-testable → SKIP `test()`.** No test exercises those
  bytes, so the run is vacuous. This covers every docs-only diff *including a whole
  vignette commit* — `.Rmd.orig`, the re-knit `.Rmd`, and `_files/` figures are
  documentation, not code/data/roxygen. Verification is still sized to the change:
  confirm `git diff --stat` shows only docs/vignette bytes moved. `check()` re-knits
  + runs the full suite at milestones — that is where vignette-touching code gets its
  test coverage, not the per-commit gate.

While iterating, run only the relevant `test_file()` / `test(filter=)`. The full
`test()` is a **gate action fired once** before a code/data/roxygen commit — never a
mid-work "did it pass" reflex. `check()` runs at milestones only.

### At the work-item milestone

Run the **`/done` skill** (the command-output reconciliation gate: grep stale
cross-refs across code + memory, confirm the LOG NOW block and a dated history entry
are current, confirm git is clean + pushed). The principle: satellite state (LOG,
memory, cross-references) goes stale silently, so it must show up in *command
output*, not be asserted. Run `check()` at the milestone too.

## Test discipline

A permanent test must be able to **fail for a real reason**. Before trusting any new
invariant/value test, run a **negative control**: feed a deliberately-wrong input and
confirm the assertion goes RED. If it can't go red, it's a tautology — drop it. A
"numbers didn't change" refactor-invariance check is a one-time scratchpad artifact,
**not** a permanent suite test. Invoke the **`/test-review` skill** when writing or
reviewing tests — it holds the full checklist (permanence filter, negative control,
anti-patterns, loop discipline).

## Migration (essentially done — read only when resuming)

`tmp_src/` is a **read-only frozen snapshot** of the source project that is *not*
part of this package (it is `.Rbuildignore`d). **Never edit inside `tmp_src/`** —
copy out of it only. Most planned migration is done; new work is developed directly.
When migration resumes for a specific piece, **invoke the `/migration` skill** — it
holds the mining procedure (`dev/find-source.sh`, hidden dirs), the copy rules, and
the restructure-and-verify-vs-archive discipline. Provenance honesty is the global
honesty rule applied to code origin: never claim code is new / migrated-from-X /
has-no-source without having *just* searched to confirm.

## House conventions

- **Naming:** `snake_case`, following the source.
- **`dr2`-family index-signature convention.** Every `dr2`-family name msamodel
  *creates* is `dr2_<indices>`: one underscore, then exactly the free indices the
  object spans, in order **response index (`i` site / `n` mode), mutated site `j`,
  mutation `m`**. A reduction over an axis drops that axis's letter (`dr2_ijm` →
  `dr2_i`). Transform prefixes (`l`=log, `n`=normalized) and source labels
  (`_msa`/`_obs`) are kept, index still underscore-set: `lrmsd_i`, `nlrmsd_i_msa`,
  `dr2_i_msa`. This governs names msamodel **creates**, not what it **calls** — call
  penm's own names (`penm::delta_structure_dr2i`) directly; never wrap or rename a
  dependency. Applies to future motion-arm names (`dh_*`, `nh_*`) too.
- **Roxygen on every function.** `#' @export` for public API, `#' @noRd` for internal
  helpers. Document `@param`, `@return`/`@returns`, add `@family` tags to group
  related functions. `@examples` in `\dontrun{}` when they need real data.
- **Organize `R/` files by function family/role** (`@family`: spm / model / objective
  / fitting / decomposition / setup), not by source filename.
- **Imports live in one place:** `R/msamodel-package.R` carries the `"_PACKAGE"` doc
  and the `@importFrom` directives. Re-export `%>%` via `#' @importFrom magrittr %>%`.
- **Tibbles** (not data.frames) for tabular returns; tidyverse for data manipulation.

## Development commands

```bash
# INNER LOOP (constant) — load_all + targeted tests; this is the working loop
Rscript -e "devtools::load_all()"
Rscript -e "testthat::test_file('tests/testthat/test-msa-mcmc.R')"   # one file
# (or devtools::test(filter='msa-mcmc') for a name-filtered subset)

# Document (regenerate NAMESPACE + man/) — ONLY after a roxygen/@importFrom change
Rscript -e "devtools::document()"

# BEFORE A CODE/DATA/ROXYGEN COMMIT (rigor gate) — one full run of the whole suite
Rscript -e "devtools::test()"

# AT A MILESTONE ONLY — full check (v0.1 baseline: 0 errors, 1 warning, 2 notes,
# deliberately accepted, GitHub-only not CRAN; see dev/LOG.md).
Rscript -e "devtools::check()"

# Install locally (only when a vignette needs the working tree's new functions)
Rscript -e "devtools::install()"

# Preview a pre-rendered vignette as standalone HTML, WITHOUT touching the repo
Rscript dev/preview-vignette.R mode-analysis   # -> dev/preview/<name>.html (git-ignored)
```

NAMESPACE and everything under `man/` are roxygen-generated — edit the roxygen
comments and run `document()`, never hand-edit them.

### Vignettes — HARD RULES

- **Never commit/push `vignettes/` without the user's explicit HTML approval.** Do
  NOT commit any change touching `vignettes/` (`.Rmd`, `.Rmd.orig`, `_files/`) until
  the **user** has said they previewed the rendered HTML and approved it. Re-knitting,
  "tests pass", or *you* glancing at a preview do not satisfy this — only the user's
  explicit OK does. When vignettes change: regenerate **every** affected preview, tell
  the user exactly which `dev/preview/<name>.html` to open, then STOP and wait. The
  `.githooks/pre-commit` hook enforces this mechanically.
- **`.Rmd.orig` precompute pattern.** Edit `<name>.Rmd.orig` → knit to `<name>.Rmd`
  (`Rscript -e 'knitr::knit("vignettes/<name>.Rmd.orig", output="vignettes/<name>.Rmd")'`)
  → preview with `dev/preview-vignette.R` → user approves → commit. Install the working
  tree first if the vignette calls new functions (`library(msamodel)` loads the
  *installed* package). **Never render the committed `vignettes/<name>.Rmd` in place** —
  `html_vignette` renders `--self-contained`, which base64-embeds figures and then
  DELETES the `<name>_files/` dir the shipped `.Rmd` references. Use
  `dev/preview-vignette.R` (renders a temp-dir copy; leaves repo figures intact).

Embedded `znb_*` datasets are (re)generated by `data-raw/prepare_znb_data.R`. The
fixture is **frozen** — regenerate it intentionally, not incidentally.
`test-spm-generate.R` catches drift between the SPM-generation code and the fixture,
so the data-prep script and that test must stay in sync.

## Working style

- Don't hand-wave uncertainty — flag anything unverified. Only say
  "confirmed/verified/fixed" when you actually checked.
- Own contradictions across turns instead of smoothing them over.
- Push back on bad ideas rather than silently implementing them.
- Terse is good; real uncertainty should be stated explicitly.
- **Tool choice:** read files with `Read` (use `offset`/`limit` for big files). Avoid
  `sed`/`head`/`cat`/`tail`/`awk`/`echo` in Bash for reading — they trigger permission
  prompts here.

When this file and `dev/plan.md` disagree, the roadmap wins — update this file or flag
the contradiction.
