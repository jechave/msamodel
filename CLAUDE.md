# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this package is

`msamodel` implements the Mutation-Stability-Activity (MSA) model of neutral
structural evolution: it predicts site-specific structural divergence in protein
evolution from single-point-mutation scans (SPMs). **v0.1 scope** (re-scoped
2026-06-04, see `dev/plan.md` §1) is "compute divergence profiles + fit them to
data": SPM generation, structure-divergence profiles, a1/a2 grid exploration,
Bayesian parameter estimation (MCMC), and site-level Shapley decomposition.

Deferred to later versions (real code, not fossils — see `dev/plan.md` §1/§12):
the model-fit *assessment* layer (goodness-of-fit / MSA-vs-MM/MS/MA comparison
and its LOESS + error-metric utilities), proportional-*allotment* decomposition,
protein-level Shapley decomposition, the trajectory/star-tree route, and
visualization.

It depends on the `penm` package (an `Imports:` dependency) for the ENM
machinery — `set_enm()`, mutation scans, response matrices. `penm` is consumed
as an installed library; its source is not part of this repo and is not edited
here.

## Status: this is a migration in progress, not a greenfield package

The code is being migrated **out of** a paper-specific research project into a
clean package. As of this writing `R/` is essentially empty. The work is
governed by:

- **`dev/plan.md`** — the single normative blueprint (13 sections): scope, final
  file layout, DESCRIPTION, refactor rules, the data-prep script, the export
  list, the test matrix, migration order, known bugs to fix, resolved decisions.
  **Read it before doing anything.**
- **`dev/PROGRESS.md`** — human-readable history + an at-a-glance checklist of
  the §9 migration order (what's done, what's left). Read it to see where we
  are. (The precise "what's next" pointer lives in Claude project memory.)

When the plan and this file disagree, **the plan wins** — update this file or
flag the contradiction, don't silently follow the stale one.

### The source project is READ-ONLY

The migration reads from a **local frozen snapshot** of the source project that
lives inside this package:

```
tmp_src/
```

referred to as `<SRC>` in the plan. This snapshot is a one-time copy of the live
project (`.../lab/projects/MSA-neutral-structure-evolution`), made so the
migration cannot touch the original — the live path is intentionally **not** a
working/additional directory in this session, so it is unreachable from here.

`tmp_src/` is disposable and is excluded from the package build via
`.Rbuildignore`. It is a snapshot, so it goes stale: if the live source changes
during the migration, re-copy it deliberately. **Treat `tmp_src/` as read-only**
— copy *out of* it, never edit *in* it.

## Migration rules (from `dev/plan.md` §4 — minimal change, "stay close")

This is a migration, not a refactor. When copying a file from `<SRC>/R/` into
`R/`:

1. Remove `library(...)` calls — both top-level and inline-in-functions. All
   package usage goes through `Imports:` + `@importFrom`.
2. Keep function signatures **identical**. Do not rename functions or variables
   or "improve" logic in v0.1.
3. Namespacing: **Option A** by default — package-level `@importFrom` + bare
   calls (smallest diff). **Exception:** ggplot2 is in `Suggests:`, so in
   `loess_compare.R` qualify calls as `ggplot2::...` behind a `requireNamespace`
   guard (you cannot `@importFrom` a Suggested package).
4. Keep `%>%` (magrittr, re-exported).
5. **Strip the non-standard `@requires` roxygen tags** (every source file has
   them; roxygen2 doesn't recognise them). Keep standard tags.
6. Add `@export` to top-level functions per the export list in plan §6; helpers
   stay unexported (`#' @noRd`). Note `plot_loess_comparison` is **not** exported
   in v0.1 (ggplot-dependent; deferred to v0.3).
7. Fix the bugs catalogued in **plan §11** as you reach each file (duplicate
   `define_selection_grid`, inline `library()` calls, the `"lower,Quite upper"`
   roxygen typo, the `paste1`→`paste0` typo in the data-prep script).

The plan (`dev/plan.md` §4/§11) is the authoritative, verified version of these
rules — this is a summary; defer to it on specifics.

Several files are **renamed** during the copy (e.g. `msa_bayesian_workflow.R` →
`msa_workflow.R`); the full rename map is in plan §2 and the copy order is in
plan §9. No public function name changes in v0.1.

## House conventions

- **Naming:** `snake_case`, following the source. v0.1 renames nothing
  (plan §4/§6), so the existing source names are the convention.
- **Roxygen on every function.** `#' @export` for public API, `#' @noRd` for
  internal helpers. Document `@param`, `@return`/`@returns`, and add `@family`
  tags to group related functions. `@examples` wrapped in `\dontrun{}` when they
  need real data.
- **Imports live in one place:** `R/msamodel-package.R` carries the
  `"_PACKAGE"` doc and the `@importFrom` directives for the whole package.
  Re-export `%>%` via `#' @importFrom magrittr %>%`.
- **Tibbles** (not data.frames) for tabular returns; tidyverse for data
  manipulation.

## Tracking progress

**After completing any plan §9 step (or sub-step), update `dev/PROGRESS.md`
immediately:** tick the checklist item and add a one-line entry under the
current date in the log. Do this as soon as the step is done — it is tied to
*finishing a step*, not to ending the session (there is no reliable
end-of-session moment). Keeping it current means the worst case is one unticked
item, not a lost session's worth of history.

## Development commands

```bash
# Document (regenerate NAMESPACE + man/ from roxygen) — run after any roxygen change
Rscript -e "devtools::document()"

# Run all tests / a single test file
Rscript -e "devtools::test()"
Rscript -e "testthat::test_file('tests/testthat/test-msa-mcmc.R')"

# Full check — must pass clean (zero errors/warnings/notes) before v0.1 ships
Rscript -e "devtools::check()"

# Install locally
Rscript -e "devtools::install()"
```

Embedded `znb_*` datasets are (re)generated by running
`data-raw/prepare_znb_data.R` — see plan §5. The fixture is **frozen**:
regenerate it intentionally, not incidentally. `test-spm-generate.R` exists to
catch drift between the SPM-generation code and the embedded fixture, so the
parameters in the data-prep script and that test must stay in sync (plan §5).

NAMESPACE and everything under `man/` are roxygen-generated — edit the roxygen
comments and run `document()`, never hand-edit them.

## Working style (from project memory)

- Don't hand-wave uncertainty — flag anything unverified. Only say
  "confirmed/verified/fixed" when you actually checked.
- Own contradictions across turns instead of smoothing them over.
- Push back on bad ideas rather than silently implementing them.
- Terse is good; real uncertainty should be stated explicitly.
- **Tool choice:** read files with `Read` (use `offset`/`limit` for big files).
  Avoid `sed`/`head`/`cat`/`tail`/`awk`/`echo` in Bash for reading — they
  trigger permission prompts in this environment.
