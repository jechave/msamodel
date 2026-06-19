# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this package is

`msamodel` implements the Mutation-Stability-Activity (MSA) model of neutral
structural evolution: it predicts site-specific structural divergence in protein
evolution from single-point-mutation scans (SPMs).

**v0.1 is shipped; v0.2 (API cleanup) is in flight.** The package computes
structure-divergence profiles (structure × site) and fits them to data: SPM
generation, structure-divergence profiles, Bayesian parameter estimation (MCMC),
and a site-level decomposition (the `phi_*` columns — renamed from `shap_*` in
v0.2, since "Shapley" was a misnomer). v0.2 also drops the a1/a2 grid API and
switches structure/active-site input to a bio3d pdb object + a plain
`pdb_site_active` integer vector; see the `dev/plan.md` roadmap.

The package was migrated out of a paper-specific research project, and **much
capability is still unmigrated** — most notably the model's motion/mode arm (only
the structure × site cell of a larger divergence grid is implemented). The
`dev/plan.md` roadmap tracks what each future version adds.

It depends on the `penm` package (an `Imports:` dependency) for the ENM
machinery — `set_enm()`, mutation scans, response matrices. `penm` is consumed as
an installed library; its source is not part of this repo and is not edited here.

## How planning works here (two tiers)

The v0.1 plan was a single 575-line blueprint that tried to be both the durable
spec and the step-by-step execution script — so every execution surprise forced an
edit to a huge document. We deliberately split that:

- **`dev/plan.md`** — the **coarse, durable roadmap**: what each version *is* (goal
  + rough scope) and findings that must not be re-discovered (e.g. the
  precomputation property). High-level and stable; edited only when a version's
  *goal* changes. Read it for orientation.
- **`dev/PROGRESS.md`** — the checklist for the **version currently in flight**.
  Written from that version's detailed plan when the version starts; dormant
  between versions.
- **`dev/LOG.md`** — append-only history (sessions, decisions, dead ends). Holds
  the "why" behind the current state and the v0.1 detail the roadmap no longer
  carries.

**Detailed planning is done per-version, at execution time** — not up front. When a
version starts: enter plan mode, read the specific code that version touches, write
the executable detail then (and write `dev/PROGRESS.md` from it). Deep code
exploration belongs where it's actionable.

**Record decisions before coding.** This still holds — it now means *write the
per-version detailed plan before coding*, and touch the roadmap only if the
version's goal changed. Don't fold surprises into ad-hoc code; plan them first.

When this file and `dev/plan.md` disagree, the roadmap wins — update this file or
flag the contradiction.

## The migration source is READ-ONLY

Continuing migration reads from a **local frozen snapshot** of the source project
that lives inside this package:

```
tmp_src/
```

A one-time copy of the live project (`.../lab/projects/MSA-neutral-structure-evolution`),
made so the migration cannot touch the original — the live path is intentionally
**not** reachable from this session. `tmp_src/` is `.Rbuildignore`d. It is a
snapshot, so it goes stale: re-copy deliberately if the live source changes.
**Treat `tmp_src/` as read-only** — copy *out of* it, never edit *in* it.

Much unmigrated code lives in **`tmp_src/archive/`** (the model's motion/mode arm,
rates, etc.). This was previously `tmp_src/.archive/` (hidden) — renamed
2026-06-18 because the dot-prefix made `grep -r` silently skip its ~18 `.R` files,
which once caused a false "no source exists" conclusion. **Always search all of
`tmp_src/` — including any hidden dirs — before concluding a function has no
migration source.** Use `dev/find-source.sh '<pattern>'` (greps `tmp_src/`
hidden-dirs-included) as the canonical check.

## Migration rules (apply whenever migration continues)

This is a migration, not a refactor — stay close to the source. When copying a
function from `tmp_src/R/` into `R/`:

1. Remove `library(...)` calls (top-level and inline). All package usage goes
   through `Imports:` + `@importFrom` (or `pkg::fn()`).
2. Keep function signatures identical. Do not rename functions or variables or
   "improve" logic as part of a copy. (Deliberate, recorded API changes — like the
   v0.2 `phi_*` rename — are separate, planned work, not folded into a copy.)
3. Namespacing: **Option A** by default — package-level `@importFrom` + bare calls
   (smallest diff). **Exception:** a Suggested package (e.g. ggplot2) cannot be
   `@importFrom`'d — qualify calls as `ggplot2::...` behind a `requireNamespace`
   guard.
4. Keep `%>%` (magrittr, re-exported).
5. Strip the non-standard `@requires` roxygen tags (every source file has them;
   roxygen2 doesn't recognise them). Keep standard tags.
6. Add `@export` to public top-level functions; helpers stay unexported
   (`#' @noRd`).
7. Fix catalogued bugs as you reach each file — record them in the version's
   detailed plan (e.g. the tree route's `p_act`→`p_ma` bug, the `akima`→`interp`
   swap; see `dev/plan.md` findings).

**Organize `R/` files by function family/role, not by source filename**
(decided 2026-06-18, retiring the earlier "keep `tmp_src` filenames" rule). A
migrated function goes in the `R/` file for its `@family` (spm / model /
objective / fitting / decomposition / setup), regardless of which `tmp_src/` file
it came from; functions with no `tmp_src` source go in the matching family file
too. Traceability to the source is provided by `dev/find-source.sh` (content
search), **not** by filename. This does **not** loosen the provenance HARD RULES
below, nor migration rule 2 — function *names* and *signatures* still don't change
in a copy; only file placement and grouping do.

### Provenance & verification discipline (HARD RULES — a 2026-06-18 failure)

A bad "migration" of the mode functions, followed by repeatedly asserting false
provenance, motivates these. They are not optional.

- **Migration sources are ONLY `tmp_src/`** (incl. `tmp_src/archive/` and any
  other hidden dirs). Before claiming any function or logic is "new" /
  "package-native" / "written from scratch", you MUST search for it:
  `dev/find-source.sh '<fn name>'` **and** `dev/find-source.sh '<core formula>'`.
  A zero-hit search on the *name* is not proof of no source — names get changed
  in migration; the **math** does not. Search the formula.
- **penm is a DEPENDENCY, never a migration source.** Call `penm::fn()`. Never
  copy or "migrate" code out of penm — not the installed library, not the local
  `../penm` source. Reading penm to check a *signature you are calling* is fine;
  reading it to source/justify package code is not. If a computation seems to
  need penm internals, that means a `penm::` call, not a copy.
- **Migration = restructure + VERIFY, not clone-and-rename.** When migrating a
  property the user has NOT already refactored, do not assume that axis-swapping
  an existing function (e.g. site → mode) is equivalent. Reproduce the OLD
  (archive) version's numbers and assert equality — machine precision for
  deterministic code, seeded for stochastic — **before** trusting the new form.
  "Looks structurally identical" is not verification; run the comparison.
- **Provenance honesty.** Never state where code came from ("I wrote it" / "it's
  migrated from X" / "no source exists") unless you have *just* searched/read to
  confirm it. These are factual claims, governed by the same rule as
  "verified"/"fixed": no claim without a check.

## House conventions

- **Naming:** `snake_case`, following the source.
- **`dr2`-family index-signature convention (decided 2026-06-18).** Every
  `dr2`-family name msamodel *creates* is `dr2_<indices>`: one underscore between
  `dr2` and the index block, then exactly the free indices the object spans in its
  representation, letters joined, in order **response index (`i` site / `n` mode),
  then mutated site `j`, then mutation `m`**. A reduction (sum/mean/weighted-mean
  over an axis) drops that axis's letter. Examples: a per-mutant SPM list-column
  whose cell is an `i`-vector (with `j`,`m` as sibling columns) is `dr2_ijm`; its
  `[mutant × i]` matrix form keeps the same name `dr2_ijm`; the jm-averaged profile
  is `dr2_i`. No `mat`/`msa`/provenance words *in the index part*. Transform
  prefixes (`l`=log, `n`=normalized) and source labels (`_msa`/`_obs`) are kept,
  with the index still underscore-set: `lrmsd_i`, `nlrmsd_i_msa`, `nlrmsd_n_obs`,
  `dr2_i_msa`. **Scope: this governs names msamodel CREATES, not what it CALLS** —
  penm's own (no-underscore) names like `penm::delta_structure_dr2i` /
  `delta_structure_dr2n` are called directly by their own names; never wrap or
  rename a dependency function to fit this convention. This convention also applies
  to future motion-arm names (`dh_*`, `nh_*`) and all later migration.
- **Roxygen on every function.** `#' @export` for public API, `#' @noRd` for
  internal helpers. Document `@param`, `@return`/`@returns`, and add `@family` tags
  to group related functions. `@examples` wrapped in `\dontrun{}` when they need
  real data.
- **Imports live in one place:** `R/msamodel-package.R` carries the `"_PACKAGE"`
  doc and the `@importFrom` directives for the whole package. Re-export `%>%` via
  `#' @importFrom magrittr %>%`.
- **Tibbles** (not data.frames) for tabular returns; tidyverse for data
  manipulation.

## Tracking progress

- **`dev/PROGRESS.md`** — checklist for the in-flight version (bare substep titles,
  `[ ]`/`[x]`). Rewritten from that version's detailed plan when it starts.
- **`dev/LOG.md`** — append-only history, newest first. Reverts, scope changes,
  decisions, dead ends live here.

**After finishing any substep:** (a) tick it in `dev/PROGRESS.md`, and (b) add a
one-line dated entry to `dev/LOG.md`. Do this as the substep finishes — tied to
*finishing*, not to ending a session.

### Definition of Done (run this reconciliation pass before declaring any plan/substep done)

Verifying the *deliverable* (tests pass, NAMESPACE clean) is not enough — a change
also has **satellite state** that describes the codebase and silently goes stale.
A plan's own verification section never covers this, so it drifts and only surfaces
later ("get ready for exit" keeps finding it). After the deliverable is verified,
run this pass against the **whole project**, not just the files you edited:

1. **Stale cross-references.** If you renamed/moved/deleted a file, function, or
   version label, grep for the OLD name across `R/`, `tests/`, `dev/`, `man/`,
   `vignettes/*.orig`, `DESCRIPTION`, and memory — fix every dangling reference.
   (E.g. after the R/ reorg, `dev/PROGRESS.md` still named the deleted files.)
2. **`dev/PROGRESS.md`** reflects reality — ticked to match what's done; set
   **dormant** if no version is in flight; never describing a superseded state.
3. **`dev/LOG.md`** has an entry for what just happened (and, at session end, the
   commit hashes + the next step).
4. **Memory** (`MEMORY.md` + files): no entry now contradicts the change; the
   "next session" pointer is current; new non-obvious decisions are recorded.
5. **Git:** intended changes committed, nothing stray staged, pushed to `main`
   (solo repo), `git status` clean — unless the user asked you to hold off.

If any item turns up work, it is **part of this task, not a later cleanup** — do it
now. Only report done once this pass is clean.

## Development commands

```bash
# Document (regenerate NAMESPACE + man/ from roxygen) — run after any roxygen change
Rscript -e "devtools::document()"

# Run all tests / a single test file
Rscript -e "devtools::test()"
Rscript -e "testthat::test_file('tests/testthat/test-msa-mcmc.R')"

# Full check (v0.1 baseline: 0 errors, 1 warning, 2 notes — all deliberately
# accepted, GitHub-only not CRAN; see dev/LOG.md)
Rscript -e "devtools::check()"

# Install locally
Rscript -e "devtools::install()"

# Preview a pre-rendered vignette as standalone HTML, WITHOUT touching the repo
Rscript dev/preview-vignette.R dr2n-analysis   # -> dev/preview/<name>.html (git-ignored)
```

**HARD RULE — never commit/push `vignettes/` without the user's explicit HTML
approval.** Do NOT `git commit` or `git push` any change touching `vignettes/`
(`.Rmd`, `.Rmd.orig`, or `_files/`) until the **user** has said they previewed the
rendered HTML and approved it. Re-knitting it, "tests pass", "check() is clean",
or *you* glancing at a preview do **not** satisfy this — only the user's explicit
OK on the HTML does. When vignettes change: regenerate **every** affected preview
(not just one), tell the user exactly which `dev/preview/<name>.html` files to
open, and then STOP and wait. (Violated 2026-06-18 — committed re-knit vignettes
the user never saw; this gate exists because of that.) The
`.githooks/pre-commit` hook enforces this mechanically (see below).

**Vignettes use the `.Rmd.orig` precompute pattern.** Workflow when changing a
vignette: edit `<name>.Rmd.orig` → knit it to `<name>.Rmd`
(`Rscript -e 'knitr::knit("vignettes/<name>.Rmd.orig", output="vignettes/<name>.Rmd")'`)
→ **preview** with `dev/preview-vignette.R` → **user reviews the HTML and approves**
→ only then commit. Knitting requires the working tree to
be **installed first** if the vignette calls new package functions
(`library(msamodel)` in the chunk loads the *installed* package). The
`.githooks/pre-commit` guard blocks committing an edited `.orig` without its
re-knit `.Rmd`, and (see the approval gate above) blocks committing `vignettes/`
changes without a recorded user approval. **Never render the committed
`vignettes/<name>.Rmd` in place** — `html_vignette` renders `--self-contained`,
which base64-embeds the figures and then DELETES the `<name>_files/` dir the
shipped `.Rmd` references. Use `dev/preview-vignette.R` (renders a temp-dir copy;
leaves the repo figures intact).

Embedded `znb_*` datasets are (re)generated by running
`data-raw/prepare_znb_data.R`. The fixture is **frozen**: regenerate it
intentionally, not incidentally. `test-spm-generate.R` catches drift between the
SPM-generation code and the embedded fixture, so the parameters in the data-prep
script and that test must stay in sync.

NAMESPACE and everything under `man/` are roxygen-generated — edit the roxygen
comments and run `document()`, never hand-edit them.

## Working style (from project memory)

- Don't hand-wave uncertainty — flag anything unverified. Only say
  "confirmed/verified/fixed" when you actually checked.
- Own contradictions across turns instead of smoothing them over.
- Push back on bad ideas rather than silently implementing them.
- Terse is good; real uncertainty should be stated explicitly.
- **Tool choice:** read files with `Read` (use `offset`/`limit` for big files).
  Avoid `sed`/`head`/`cat`/`tail`/`awk`/`echo` in Bash for reading — they trigger
  permission prompts in this environment.
