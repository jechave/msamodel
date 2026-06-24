# Log — msamodel migration

Append-only history of what was done and attempted (newest first), including
dead ends and decisions. The durable roadmap is `dev/plan.md`; `dev/PROGRESS.md`
is the checklist for the version currently in flight. This file is the history
those two don't keep.

One short entry per working session.

### 2026-06-24 — Slice 2: ML fit arm (`fit_msa_ml`) + MCMC-vs-ML vignette section

Second of the three 0.3.0 fit slices. Added a maximum-likelihood **point**
estimator parallel to the MCMC, then a vignette section comparing the two arms and
the observed profile.

- **`R/fitting.R` — `fit_msa_ml()`** (new file, `@family fitting`). Maximises the
  slice-1-corrected `calculate_loglik_msa` over the *same* coords as the MCMC
  (`a1`, `b = log2(a2+1)`) and box bounds. L-BFGS-B, deterministic coarse grid-max
  start (25×25), `optimHess` covariance at the optimum, delta-method `se_a2`,
  **fail-loud** on a singular / non-PD Hessian (no silent NA SEs). Returns a list
  (a1, a2, logLik, sigma_hat, cov, se_a1, se_a2, convergence, par_fit). Added
  `optim`/`optimHess` to `@importFrom stats` and two NSE globals
  (`nlrmsd_obs`, `nlrmsd_i_msa`).
- **Provenance** (find-source.sh, name + formula): NO `tmp_src` source for
  ML-fitting `(a1,a2)` — source project fit only by MCMC / grid; the lone `optim`
  (`archive/R_backup/model_rates.R`) is the unmigrated rate arm, a different object.
  ⇒ new code, not a migration. Recorded in `dev/plan.md` §v0.3 (B).
- **Cross-checks on znb.** optim optimum (logLik −138.546) ≥ 81×81 grid max
  (−138.632), same basin → refines off-grid. ML point (0.458, 42.30) lies inside
  the seeded (seed 2024, 4000/1000) MCMC posterior, near its mean; ML SEs
  (0.122, 10.19) ≈ posterior sd (0.131, 11.3).
- **Tests** `test-fit-ml.R` (33 assertions): list shape, grid agreement,
  MCMC-mode proximity, frozen regression, cov/SE shape, fail-loud paths, pdb_site
  contract parity. Suite 116/0F (+33); `check()` at v0.1 baseline (0E/1W/2N).
- **Vignette `msamodel-intro` §5.2.** Compares MCMC fit, ML fit, observed profile.
  *Design correction after user review:* the §5.1 MCMC profile is the posterior
  *predictive mean* (per-site prediction averaged over the posterior); pairing that
  against an ML *point* prediction conflates the estimate difference with the Jensen
  gap. So §5.2 predicts BOTH arms at a single `(a1,a2)` — posterior mean for MCMC,
  point estimate for ML — making the only difference the estimate. Adds a param
  table (both R²=0.605), a 3-line profile overlay, a by-method faceted scatter, and
  a **direct MCMC-vs-ML scatter** (method-vs-method R² = 0.9999 → arms agree
  site-for-site; Jensen gap negligible on znb). User approved the rendered HTML;
  committed under `VIGNETTE_APPROVED=1`.
- **Knit-cwd footgun:** first knit from repo root dumped figures to a stray
  root-level `msamodel-intro_files/` (relative `fig.path` resolved against cwd).
  Re-knit with `setwd("vignettes")` so figures land in `vignettes/…`; removed the
  stray dir. The shipped `.Rmd` had correct refs; only the PNG location was wrong.
- Committed `aa868a3`, pushed `b7251be..aa868a3`. **Next:** slice 3 — `lrmsd_n`
  (mode) fit on seeded synthetic observed data (plan it when it starts).

### 2026-06-24 — Architecture decisions + start σ correctness fix (slice 1)

Planning session, then began the first of three 0.3.0 *fit* slices.

Decisions recorded (in `dev/plan.md` §v0.3, split the old "observed profiles + fit"
bullet):
- **Observed-pattern computation (lrmsd from homologs + alignment) is OUT of
  `msamodel`** → a future separate "protein-evolution-patterns" package (seq /
  structure / motion), serving other models too. `msamodel` keeps "observed data =
  a vector". No homolog/alignment source exists in `tmp_src` to migrate; the EC2024
  pipeline and the user's earlier empirical-model project kept it as a separate
  upstream stage. (I wrongly claimed "no source" early in the session from a
  name-only `.R` grep — corrected after the user noted `.Rmd`/other-project code; it
  doesn't change the architecture conclusion.)
- **New ML fit arm wanted**, agile scope: estimator only (no S3, no `method=` flag,
  separate fn from MCMC). Speed + `lm`/`gam`-style direction. Shared S3 convention
  deferred to a later pass.
- **`lrmsd_n` fit** will bootstrap on **seeded synthetic** observed data until the
  patterns package exists; mark synthetic on the data, not via runtime warning.
- Slice order: **σ fix → ML `lrmsd_i` → `lrmsd_n` synthetic fit.** Each slice
  planned in detail only when it starts (two-tier rule).

Slice 1 (σ fix) — in flight:
- **The bug:** `calculate_loglik_msa` (`R/objective.R:56`) profiles σ out but used
  `sd(residuals)` (divisor `n−1`); the profiled-Gaussian MLE is
  `sqrt(mean(residuals^2))` (divisor `n`). Wrong constant in the likelihood the MCMC
  has been sampling.
- **Math verified (not asserted):** numeric maximiser of `sum dnorm(r,0,σ,log=T)`
  matches `sqrt(mean(r²))` (not `sd(r)`); and argmax over θ is identical for SSR,
  loglik-n, loglik-(n−1) — so the **`(a1,a2)` point estimate is unchanged**; only
  σ-derived quantities (logLik, SEs, posterior width) move.
- **Verified outcome (real znb data):**
  - Deterministic grid check (a1×a2): OLD/NEW argmax IDENTICAL (0.5, 31); logLik
    change is a **constant +0.00111 shift** across the whole surface (Spearman cor
    = 1) — because σ profiled out makes loglik = const − (n/2)log(σ̂²), and the
    n vs n−1 divisor is an (a1,a2)-independent constant.
  - Seeded MCMC OLD vs NEW (seed 2024, 4000/1000): posterior **bit-identical** —
    `a1`/`a2` samples `identical`, max summary diff = 0. A constant logLik shift
    leaves the Metropolis acceptance ratio (differences only) unchanged, so the
    posterior cannot move. ⇒ no PAUSE, nothing downstream stale.
  - Blast radius was NOT zero on tests (my earlier prediction was wrong): 3 tests
    froze the *absolute* logLik value and moved by the constant +0.0011. Updated as
    corrected-expected (not drift): `test-msa-evaluate.R` literal → −184.3230779142;
    `test-contract.R` manual route `sd(res)`→`sqrt(mean(res^2))` (keeps it a true
    independent route); `profile-invariance` snapshot accepted. Suite 83/0F;
    check() at v0.1 baseline (0E/1W/2N). Intro vignette unaffected (no MCMC, no
    posterior numbers printed) — commit-gate not triggered.
- Files touched: `R/objective.R`, `tests/testthat/{test-msa-evaluate,test-contract}.R`,
  `tests/testthat/_snaps/profile-invariance.md`, plus `dev/{plan,PROGRESS,LOG}.md`.
- NOT committed yet — awaiting user review (commit slice 1 alone, then ML arm).

### 2026-06-19 — Mode nested-models fn + dr2n vignette parity with intro

Brought the `dr2n-analysis` vignette to parity with the intro vignette's site §3.
- New `calculate_lrmsd_n_nested_models(spm_pp, a1, a2)` in `R/model.R` — mode
  counterpart of `calculate_lrmsd_i_nested_models`: returns `n` + `lrmsd_n_mm/ms/
  ma/msa`, NO `pdb_site` (modes aren't residue-anchored). The pure
  `calculate_msa_decomposition` is reused unchanged on the four `lrmsd_n` vectors
  (the payoff of the axis-agnostic refactor). Globals get `lrmsd_n_*` + `n`. New
  test in `test-msa-mode.R` pins the recipe (independent route). Suite 83/0F,
  check() at v0.1 baseline.
- Vignette restructured to mirror the intro flow: profile → **nested models +
  decomposition (new)** → sweeps (moved after decomposition). y axes show
  `lrmsd_n = log(sqrt(dr2_n))` everywhere (computed inline; `calculate_dr2_n_msa`
  unchanged); all plots show ALL modes (removed the `n<=40` sweep filter).
  Removed the "Consistency with the site form" section (per user). Four-model plot
  uses the site model colours; mode decomposition plot uses **free_y** panels
  (user: shared scale hid the small stab/act components) + paper component colours.
- User reviewed `dev/preview/dr2n-analysis.html` and approved.

### 2026-06-19 — REFACTOR: pure-vector decomposition + nested-models fn + `_i_` rename

Triggered by a vignette request (decompose the profile at fixed `(a1,a2)=(2,500)`).
Surfaced that `calculate_msa_decomposition` was coupled to hard-coded column names.
Per user's design principle (see [[pure-functions-over-param-flexibility]]): make
the math a **pure 4-vector function** rather than adding an axis/key argument
(false flexibility). Changes:

- `calculate_msa_decomposition(mm, ms, ma, msa)` → returns a named list of the
  three phi vectors. Context-free (no tibble/columns/axis); serves site `i` and
  future mode `n` unchanged. `ma` kept as 4th arg (unused by sequential, reserved
  for Shapley). `calculate_decomposition_samples` keeps its tibble interface and
  calls the vector fn on its columns; `_summary` unchanged.
- New `calculate_lrmsd_i_nested_models(spm_pp, a1, a2)` in `R/model.R` — single
  source of truth for the four-variant lrmsd recipe (MM/MS/MA/MSA). Refactored
  `calculate_prediction_samples` to call it per posterior sample (verified
  numerically identical to the old four-block+pivot version — equivalence check,
  bit-for-bit as plain data.frames).
- Full rename `lrmsd_mm/ms/ma/msa` → `lrmsd_i_mm/...` everywhere (prediction_samples
  cols, objective.R `lrmsd_i_msa`/`nlrmsd_i_msa`, prediction_summary `variable`
  values, package globals, tests). Parallels `dr2_i`/`dr2_n`; leaves room for the
  mode arm's `lrmsd_n_*`.
- Tests: new pure-vector unit test, nested-models recipe test; updated fixtures to
  `lrmsd_i_*`. Full suite 77/0F. `check()` at v0.1 baseline (0E/1W/2N).
- Vignette: new §3 "Decomposing the profile at given (a1,a2)" (fixed-(a1,a2),
  before fitting) with a four-model line plot + the decomposition plot; §3-6
  renumbered. Paper component colours applied to both decomposition plots
  (mut #FF8C00 / stab #0000CD / act #C41E3A, from tmp_src shap_decomposition
  chunks); four-model plot MM=orange/MS=blue/MA=green #1B9E77/MSA=red. **Vignette
  HELD from commit pending user HTML approval** (commit gate).

### 2026-06-19 — BUGFIX: site decomposition was the wrong formula (Shapley → sequential)

User caught that the migrated `phi_*` site decomposition shipped the **wrong
formula**. `tmp_src` has two decompositions of `lrmsd_msa`: the symmetric
**Shapley** form (`tmp_src/R/msa_decomposition.R`, `methods.md:101-103`) and the
**sequential/additive** M0→MM→MS→MSA form (`methods.md:106-111`,
`shared_data_preparation.R:48-51`, the paper-analysis one literally named `phi_*`).
The migration copied Shapley and then renamed its columns `shap_*`→`phi_*` —
attaching the `phi_*` name to the wrong math. Both telescope to `lrmsd_msa`, so a
sum check passes either way → bug invisible, survived migration + rename + a memory
that recorded the Shapley formula as "what it actually is."

Fix: `R/msa_decomposition.R` now computes `phi_mut=lrmsd_mm`,
`phi_stab=lrmsd_ms-lrmsd_mm`, `phi_act=lrmsd_msa-lrmsd_ms`. Per user, `lrmsd_ma`
kept as a required input (unused by sequential) reserved for a future `method`
flag (sequential vs shapley); MA still computed in `calculate_prediction_samples`.
Added a value-pinning test (sample-1-site-1: seq stab=.1/act=.3 vs Shapley
.15/.25) that fails under the old formula. Updated roxygen+man, NEWS, and the
`decomposition-not-shapley` memory. Full suite 68/0F; `check()` at v0.1 baseline.

Vignette `msamodel-intro` (§5 formulas+prose) re-edited & re-knit (LaTeX now
sequential, Shapley "average of marginal effect" prose dropped); preview at
`dev/preview/msamodel-intro.html`. **Vignette HELD from commit pending user HTML
approval** (commit gate); code committed separately.

### 2026-06-18 — VERSIONING FIX: retire the lettered sub-version scheme

The "v0.3a / v0.3b / v0.3c" lettered sub-version scheme was a mistake and is
retired. **Versions are plain semver in `DESCRIPTION`** (`MAJOR.MINOR.PATCH`,
`.9000` while in development); a version is a *release*, not a work step. **Work is
named, not sub-versioned** — the roadmap lists whole versions and describes each
release's work as named bullets (no letters, no "phase N"). v0.1/v0.2/v0.4/v0.5
were already lettered-free; only 0.3 had drifted.

Worst damage (the letters leaking into `DESCRIPTION`/`NEWS`, the package's version
namespace) was already fixed earlier this session: `DESCRIPTION` is `0.3.0.9000`
and `NEWS.md` is one `# (development version)` entry grouped by change type. This
pass cleaned the rest: `dev/plan.md` (v0.3 block reworded to named work items),
`dev/PROGRESS.md` (letters removed + "version in flight" vocabulary → "work item"),
`CLAUDE.md` (convention dated, not labelled), one stray `NEWS.md` line, and memory.

**This LOG is append-only history — the 24 older `v0.3a/b` mentions below are left
intact** (they were accurate with the label used at the time). Not rewriting git
history either: the pushed commits whose messages say "v0.3b" stay as-is
(rewording them needs force-push). Root cause for the record: I adopted a
non-standard versioning scheme the user floated without questioning it, and should
have flagged it as non-standard up front.

### 2026-06-18 — v0.3b SHIPPED: dr2* naming convention applied package-wide

Executed the rename below in full. Order: locked the profile-invariance guard
FIRST (commit `f068a45`), then renamed code+fixture (`6115b7d`), then
vignettes/docs/version/persistence (final commit this session).

**Result — every verification gate green:**
- Profile-invariance snapshots (captured pre-rename) reproduce **bit-for-bit**
  across BOTH the code rename AND the `znb_spm` regen — no numeric value moved.
- `znb_spm` regen: Validation OK vs tmp_src original (tol 1e-8) — data unchanged,
  names only.
- Full suite 62/0F. `check()` = v0.1 baseline (0E/1W/2N) after removing stray
  top-level `*_files` (knit must run from within `vignettes/`, not repo root, or
  `fig.path` pollutes the package root — gotcha worth remembering).
- Code surface (R/ tests/ data-raw/ vignettes) clean of old tokens; only penm::
  calls + intentional NEWS history references remain.
- DESCRIPTION 0.2.0.9000 → 0.3.0.9000. Convention recorded in CLAUDE.md,
  dev/plan.md, memory ([[dr2-naming-convention]]).

Notable: `data-raw/prepare_znb_data.R` validation needed `names(orig)` `dr2`→
`dr2_ijm` remap so the shared-column check compares CONTENTS not the old name.
NEWS v0.3a section + one v0.2.0 line had stale old names → fixed.

### 2026-06-18 — v0.3b decision detail: dr2* naming convention + profile guard

Long design discussion settled the `dr2*` naming mess. **Decision:** msamodel adopts
one convention for every `dr2`-family name it *creates* — `dr2_<indices>` (one
underscore, then free indices `i`/`n`,`j`,`m` joined; reductions drop the averaged
letter). penm is the only internally-consistent codebase but uses no-underscore
(`dr2i`), and msamodel's transform/source layers (`l`/`n`, `_msa`/`_obs`) make
penm's run-together form ambiguous (`nlrmsdi_msa`) — so we diverge deliberately. Key
rule: **the convention governs names msamodel CREATES, not what it CALLS** —
`penm::delta_structure_dr2i/n` are called directly by their own names, NO wrappers
(rejected pass-through adapters as indirection-without-behavior). Local
`delta_structure_dr2` is a verified-bit-identical penm duplicate → deleted. Exported
fns normalized too: `calculate_dr2i_msa`→`calculate_dr2_i_msa`,
`calculate_dr2n_msa`→`calculate_dr2_n_msa` (breaking; GitHub-only/solo/pre-1.0).
`_msa`/`_obs` source labels stay (not index sigs). Profile cols `dr2_i`/`dr2_n`
already correct.

Per project rule (fix-plan-first), updated `dev/plan.md` v0.3b bullet BEFORE coding.
dev/plan.md quantity table (line ~163) already used `dr2_ijm`/`dr2_njm` — consistent.

**Step 0 done (guard locked before any rename):** added
`tests/testthat/test-profile-invariance.R` snapshotting the `dr2_i`/`dr2_n` profiles
+ loglik at (a1=2,a2=5) via `expect_snapshot_value(style="serialize")` — verified
json2 is NOT lossless (~3e-15) so used serialize. Recorded
`_snaps/profile-invariance.md` from CURRENT pre-rename code; stable on re-run (PASS
3, no "Adding"). This is the numeric net: any rename that moves a profile value will
fail. Committing the guard before touching names.

### 2026-06-18 — SESSION END: reorg + roadmap committed & pushed; next = v0.3b

Committed & pushed both pieces of this session to `main`:
- `56dac7d` — the R/ reorg (rule retired, files grouped by `@family`, evaluation
  split into model+objective; see the entry below for detail).
- `992e771` — roadmap reshape in `dev/plan.md`: **v0.3 is now the complete
  site+mode structural-divergence suite** (model AND observed, predict AND fit),
  NOT motion. v0.3a SHIPPED; **v0.3b = consistency tidy** (`dr2→dr2i`; use
  `penm::delta_structure_dr2i`); v0.3c+ = TBD step-by-step (observed dr2i/dr2n from
  homologous structures + alignment; fit dr2n; joint fit). **Motion promoted to its
  own v0.4**; tree route renumbered v0.4→**v0.5**.

Tree clean, in sync with origin/main. PROGRESS.md set dormant. **Next session:
start v0.3b** — enter plan mode, read `R/spm.R` + `data-raw/prepare_znb_data.R` +
the drift test, write the detailed plan (the `dr2→dr2i` rename + `penm::` helper
swap need a deliberate `znb_spm` fixture regen). Memory `project-next-session`
points here.

### 2026-06-18 — File-rename rule RETIRED; R/ reorganized by @family; evaluation split

**Decision (with user): retire the "keep `tmp_src` filenames" migration rule.**
New rule — organize `R/` files by function `@family`/role, not by source filename;
sourceless functions go in the matching family file; `dev/find-source.sh` (content
search) provides traceability instead of filenames. Updated CLAUDE.md (migration
rules), dev/plan.md (resolved decisions), and the `feedback-no-renames` memory
(now states the *opposite* of its old self). Provenance HARD RULES and the
no-function-rename rule are untouched — only *file placement* is freed.

**Naming fix (with user): `@family evaluation` was misleading** ("evaluation" reads
as *assessing*; these are the model's *forward* computation). Split into:
- `@family model` — `calculate_dr2i_msa` / `calculate_dr2n_msa`: the pfix model's
  forward map (predicted divergence profile at given `(a1,a2)`). NOT `prediction`
  — that word is taken by the fit-level `calculate_prediction_*` (posterior/ML
  profile).
- `@family objective` — `calculate_loglik_msa`: one fitting *criterion* (future
  RMSE/robust/other-likelihoods join here). Single-member family for now, by
  design.

**Reorg (pure moves + tag edits, no logic change):**
- `R/spm.R` ← `generate_spm_data` + `delta_structure_*` (`@noRd`) + `preprocess_spm`
  + `preprocess_spm_mode` (consolidates old `generate_spm_data.R` +
  `msa_bayesian_data_preparation.R`, both deleted).
- `R/model.R` ← `calculate_dr2i_msa`, `calculate_dr2n_msa`.
- `R/objective.R` ← `calculate_loglik_msa`. (Old `msa_model_evaluation.R` deleted.)
- `msa_bayesian_workflow.R` / `msa_bayesian_analysis.R` kept (genuinely MCMC,
  `@family fitting`).
- **Boundary recorded so it isn't re-litigated:** `spm` = per-mutant *raw*
  divergence (the data); `model` = pfix-weighted *profile* (the forward map).

**Verification (all green):** baseline captured first (16 exports, 23 man pages,
59 tests pass). After: `document()` → NAMESPACE unchanged, man/ file set unchanged
(only the 6 moved pages' content changed: provenance line + the evaluation trio's
`\concept`/`\seealso` retag `model`/`objective` + first-function title-header
text). `test()` 59 pass / 0 fail (= baseline). `check()` 0 errors / 1 warning /
2 notes (= accepted v0.1 baseline, no new diagnostics). `grep @family evaluation
R/` = 0. `find-source.sh 'preprocess_spm'` and `'calculate_dr2i_msa'` still resolve
to `tmp_src/` despite the new filenames — proves the rule-change premise.
**Visual vignette regression:** installed the refactored tree, re-knit both
`.Rmd.orig`; the only diff was the `sessionInfo()` version banner
(`0.2.0` → `0.2.0.9000`, a precompute artifact) — zero numerical/figure change, so
reverted the version churn. Previews in `dev/preview/` for eyeballing.

### 2026-06-18 — SESSION END / NEXT-SESSION AGENDA

**Pick up here next session.** Open topic the user wants to discuss (NOT yet
decided — do not implement, discuss first):

- **Revisit the "don't rename files when migrating" rule** (CLAUDE.md migration
  rule + [[no file renames]] memory). The rule's original purpose was
  traceability: keep migrated files matching `tmp_src/` filenames 1:1. But this
  session exposed two cracks: (a) the migration ALREADY routinely moves functions
  into files that don't match their `tmp_src` origin, so strict file-to-file
  correspondence is already broken; (b) it produced a genuine mis-housing —
  `preprocess_spm` / `preprocess_spm_mode` (general SPM reshaping, `@family spm`)
  and `calculate_dr2*_msa` (per-point prediction, `@family evaluation`) sit in
  `R/msa_bayesian_data_preparation.R` / live under `msa_bayesian_*` names, even
  though they are NOT Bayesian-specific. The `@family` roxygen tags already
  encode the *real* grouping (spm / evaluation / fitting / decomposition /
  setup / datasets), independent of filenames.
- **The proposal to weigh:** shift the mental model from "migrate whole FILES
  (keep names)" to "migrate FUNCTIONS into well-named, function-grouped files"
  (e.g. a neutral `spm_preprocess.R`, a pipeline-neutral evaluation file, leaving
  `msa_bayesian_*` for genuinely-MCMC code). Trade-off: lose easy
  file-to-tmp_src diffing (mitigated now by `dev/find-source.sh`, which finds a
  source regardless of filename) vs. gain a layout that reflects architecture.
- Decide: keep the rule, retire it, or replace it with "group by @family /
  function role; filenames need not match tmp_src." If retired, a follow-up
  reorg task (pure moves + document() + tests) would realign existing files.

**State at session end:** working tree clean, in sync with origin/main
(`a0dba7e`). Tests 59 pass. v0.3a mode code (`calculate_dr2n_msa`,
`preprocess_spm_mode`) is STAYING — user confirmed the dr2n-analysis vignette
HTML looks correct, so no further verification/revert. The migration guardrails
(unhide archive, find-source.sh, provenance rules) are in place.

### 2026-06-18 — PROCESS FIX: unhide archive + migration-discipline guardrails

After a serious failure (see below), added guardrails so it can't recur.

- **The failure:** Claude "migrated" the v0.3a mode functions poorly and then
  asserted false provenance repeatedly — claimed `calculate_dr2n_msa` was
  "package-native, no tmp_src source" when the source
  (`calculate_diff_mode_msa.R`) was sitting in `tmp_src/.archive/`. Root causes:
  (1) `.archive/` is a HIDDEN dir, so `grep -r tmp_src/` silently skipped its 18
  `.R` files and returned 0 hits → false "no source" conclusion; (2) confused
  penm (a dependency to CALL) with a migration source — even read penm internals
  to reason about provenance; (3) cloned the user's already-refactored `dr2_i`
  function, swapped `i`→`n`, and shipped it WITHOUT verifying against the
  archive's `dr2_n` — skipping the refactor-and-check the user would have done.
- **Fix 1 — unhid the archive:** `mv tmp_src/.archive tmp_src/archive`. Now
  `grep -rln dr2_njm tmp_src/` returns 3 files (was 0). Pure local rename
  (tmp_src is git-ignored + Rbuildignored). Updated all `.archive` refs in
  plan.md/LOG.md.
- **Fix 2 — `dev/find-source.sh`:** mechanical "does this have a migration
  source?" check; greps all of tmp_src/ incl. hidden dirs. Run before claiming
  any code is new.
- **Fix 3 — CLAUDE.md migration rules + memory:** migration sources are ONLY
  tmp_src/; penm is a dependency never a source; migration = restructure + VERIFY
  (reproduce old numbers), not clone-and-rename; never state provenance unchecked.
- **penm access:** kept (rule-only guardrail, no hard block) — user's call.
- **STILL OPEN (not done here):** the v0.3a mode functions on `main` remain
  UNVERIFIED against the archive. Verifying / re-migrating them, and deciding
  whether to revert the v0.3a commit, are separate user-directed follow-ups.

### 2026-06-18 — v0.3a started: mode-form structural divergence (`dr2_njm`)

- Split the roadmap's large v0.3 into slices; **v0.3a** is the first/smallest:
  add the mode-form of `dr2` only (`dr2_njm` / `dr2_n`), introducing the mode
  response axis with no new quantity and no motion. Predict-only.
- Design forks resolved in conversation before planning (the substance of the
  session): (1) keep the **slow per-mutant loop**, not penm's fast `smrs`
  (the `C·F` matrix batch) or `amrs` (analytical) — the loop is the general
  engine that survives into motion (v0.3b) and trees (v0.4) and is the
  benchmark; `dr2n` is just one more cheap reduction of the same `dr`.
  (2) **Raw SPM stays one "all-effects" object** — `dr2n` is another SPM
  column next to energies + `dr2`; the user's framing: energy changes are
  mutation effects too, and all effects of one mutation share the `(j,m)` key.
  (3) **Site/mode separation lives at the preprocess/evaluate layer**, which is
  a reshape-for-purpose step, not in the raw object: new `preprocess_spm_mode()`
  + `calculate_dr2n_msa()`; a mode analysis never touches `dr2mat`/`site_map`.
  (4) **Predict-only** — no fit (no observed mode data), v0.1 site fit untouched.
- Vignette: new **`dr2n-analysis`** topic vignette (best practice = vignettes by
  task, not feature accumulation; intro stays the narrow on-ramp). v0.3a content
  = predict-only `a1`/`a2` sweep of the `dr2_n` profile; a fitting section is
  added there later when alignment-derived `dr2_n` + a mode fit exist.
- Verified pre-planning (penm local source at `../penm`):
  `delta_structure_dr2n(wt,mut)` returns 678 finite per-mode values on `znb_wt`
  (678 = 3·228−6), same call shape as `delta_structure_dr2`; `get_mode` is NOT
  exported (use `penm:::get_mode`, as the archive did). Memory written:
  `project-penm-mutscan-tiers`, `project-msamodel-v03-mode-axis`.

### 2026-06-18 — vignette preview tooling + .Rbuildignore fix (post-v0.3a)

- **Problem:** rendering a committed `vignettes/<name>.Rmd` to HTML in place
  DESTROYS its figures — `html_vignette` uses pandoc `--self-contained`, which
  base64-embeds the figures then deletes the `<name>_files/` dir the shipped
  `.Rmd` references. Hit during the v0.3a commit (figures vanished, nearly
  committed missing).
- **Fix:** `dev/preview-vignette.R` — renders a COPY of the `.Rmd` + `_files/`
  in a temp dir and copies only the self-contained `.html` back to
  `vignettes/<name>.html` (git-ignored). The repo's `vignettes/` is never
  written to. Verified: figures present & `git status` clean before and after.
  Usage `Rscript dev/preview-vignette.R <name>`. Documented in CLAUDE.md.
  (Output goes to `dev/preview/<name>.html`, a git-ignored scratch dir — NOT
  `vignettes/`, which should hold only files that ship.)
- **Also fixed a latent v0.3a bug:** `.Rbuildignore` had a filename-specific
  `^vignettes/msamodel-intro\.Rmd\.orig$`, so the new `dr2n-analysis.Rmd.orig`
  would have SHIPPED in the build tarball. Generalized to
  `^vignettes/.*\.Rmd\.orig$` (ignores all `.orig` sources; rendered `.Rmd`
  still ships). Verified both `.orig` ignored, both `.Rmd` not.

### 2026-06-18 — v0.3a implemented (steps 1-8 + verification)

- `generate_spm_data()` now emits `mode` + `dr2n` SPM list-columns
  (`dr2n <- penm::delta_structure_dr2n(wt, mut)` per mutant; `mode <-
  seq_along(dr2n)`). Added `preprocess_spm_mode()` (builds `dr2nmat`, no
  `pdb_site` map) and `calculate_dr2n_msa()` (predict-only; identical
  mutant-axis reweighting). `get_mode` stays `penm:::` (not exported). Exports
  14 → 16.
- Regenerated `znb_spm`. The one-time `tmp_src` validation now compares only
  shared columns (`znb_spm[names(orig)]`) — the migrated SPM legitimately has
  extra columns the 2021 original lacks; shared columns still match at 1e-8.
- Tests: added `mode == penm:::get_mode(znb_wt)` guard to test-spm-generate;
  new `test-msa-mode.R` (shape; a **Parseval** site/mode total-equality check —
  basis invariance, a real invariant not a recomputation; and an
  independent-route reweighting check). Suite 45 → 59 pass.
- New `dr2n-analysis` vignette (predict-only `a1`/`a2` sweep of the per-mode
  profile; consistency check shows site_total == mode_total). NOTE: knitting
  required `devtools::install()` first — the vignette's `library(msamodel)`
  loads the INSTALLED package, which lacked the new functions until installed.
  (Same will apply to any future vignette using new code: install, then knit.)
- Generalized `.githooks/pre-commit` to guard every `vignettes/*.Rmd.orig`
  pair, not just the intro.
- DESCRIPTION → `0.2.0.9000`; NEWS development section added.
- `check()` = 0E/1W/2N (baseline). **Finding (not acted on):** the data-size
  NOTE grew — `znb_spm` 12.5 → 23.6 MB, because `dr2n` (678 modes/row) is ~3x
  wider than `dr2` (228 sites/row). Each future v0.3 quantity adds another wide
  column to the same fixture, so fixture-size management (e.g. dropping the
  redundant `dr` raw-vector column, or shipping a smaller example) is worth a
  deliberate decision before v0.3b. Not fixed here — out of v0.3a scope.
- v0.3a COMPLETE; NOT committed (awaiting user review).

### 2026-06-11 — vignette: switch to the `.Rmd.orig` precompute pattern

- **Why:** the intro vignette had drifted and violated "fail loud". Heavy compute
  lived in a *separate* `vignettes/precompute.R` → `vignette_cache.rds`; the `.Rmd`
  showed hand-typed `eval=FALSE` snippets *resembling* that script and read the
  cache. Two hand-maintained copies of the pipeline code with nothing comparing
  them → they diverged: the `sweep-code` chunk showed only `sweep_a1` while the
  cache (and a rendered figure) also had `sweep_a2`, so a figure was produced by
  code the reader never saw. Worse, every figure chunk was gated
  `eval = have_cache && has_ggplot`, so a missing/stale cache *silently skipped*
  figures instead of erroring.
- **Fix (rOpenSci `.Rmd.orig` pattern):** `vignettes/msamodel-intro.Rmd.orig` is
  now the single source of truth — every pipeline chunk is real `eval=TRUE` and
  computes inline; the shipped `msamodel-intro.Rmd` is *pre-rendered* from it via
  `knitr::knit("vignettes/msamodel-intro.Rmd.orig", output="vignettes/msamodel-intro.Rmd")`.
  Shown code == executed code (one copy), so the drift class is gone structurally.
  Deleted `precompute.R` and `vignette_cache.rds`. The MCMC fit runs at 500 iter /
  100 burn-in (fast render). The only `eval=FALSE` chunk is the SPM-scan recipe in
  §1, explicitly labelled as the code that produced the shipped `znb_spm` (already
  drift-guarded by `test-spm-generate.R`).
- **No drift test (deliberate):** a "recompute-and-assert-equal" test (the guard
  `znb_spm` uses) is wrong for a vignette — a vignette is *meant* to change across
  versions, so such a test fails on every intentional edit. The `.Rmd.orig` pattern
  already removes the root cause (two copies). `set.seed(1024)` stays at the top of
  `.Rmd.orig` only to keep the committed render stable, not as an asserted invariant.
- **Stale-render guard = a pre-commit hook** (`.githooks/pre-commit`, enabled via
  `git config core.hooksPath .githooks`). It blocks a commit that stages
  `msamodel-intro.Rmd.orig` *without* also staging the re-rendered
  `msamodel-intro.Rmd`. It does NOT recompute/diff content (no brittleness, doesn't
  fight intentional edits) — it only enforces "re-knit was staged alongside the
  source edit." `--no-verify` bypasses for genuine no-render edits. Tested: blocks
  when only `.orig` staged; passes when both staged or `.orig` untouched. (A buried
  comment alone was no safeguard — it can't tell you that you forgot.) `.githooks`
  is `.Rbuildignore`d. The re-knit command is also in a comment atop `.Rmd.orig` for
  reference.
- **Wiring:** `.Rbuildignore` now ignores `msamodel-intro.Rmd.orig` (was
  `precompute.R`). `fig.path = "msamodel-intro_files/figure-html/"` so figures don't
  land in a bare `figure/` dir (which `R CMD check` flags as a knitr leftover —
  that was the one new NOTE, now cleared). Figures are checked in (un-ignored in
  `vignettes/.gitignore`).
- **Verify:** knit ran all chunks incl. the live MCMC (no skips); rendered `.Rmd`
  carries embedded results (`a1≈0.455`, `a2≈44.4`) and references both sweeps;
  `tools::pkgVignettes()` sees exactly one vignette; `check()` back at baseline
  (0 errors, 1 warning, 2 notes); the built `doc/msamodel-intro.html` embeds all
  10 figures as base64 (0 broken path refs). NOTE: the 2026-06-11 v0.2 entry below
  describes the now-removed `precompute.R`/cache mechanism — superseded by this
  entry, not edited (append-only).

### 2026-06-11 — v0.2 COMPLETE: all four API changes executed

- **All four v0.2 changes done; `check()` = 0 errors, 1 warning, 2 notes — same
  as the v0.1 baseline (no new issues).** Tests 55 → 45 pass (the drop removed
  the grid test file + the `get_active_site` parse test + 3 contract grid lines;
  no coverage lost — the dropped assertions were redundant with the manual i-join
  test and the prediction/decomposition pdb_site checks). Exports 18 → 14.
  DESCRIPTION bumped to 0.2.0; NEWS.md written (changes flagged breaking).
- **Step 1 (rename):** `shap_* → phi_*` in `R/msa_decomposition.R` (3 columns +
  both validation `required_cols` + `starts_with("phi_")`), `test-decomposition.R`,
  vignette prose + `labs_comp`, roxygen prose in `msa_bayesian_workflow.R`,
  `DESCRIPTION` Description field (also dropped its "grid-based exploration"
  clause), `CLAUDE.md`. NO `data/*.rda` regen. `man/*.Rd` regenerated. `tmp_src/`
  untouched.
- **Step 2 (bio3d object):** deleted `load_protein` (one-line `read.pdb`
  wrapper). Added `inherits(pdb,"pdb")` guard to `setup_enm` with a clear error
  pointing at `read.pdb`/`read.cif`. Updated callers (test-spm-generate, vignette
  precompute + Rmd) to `bio3d::read.pdb`. **mmCIF verified empirically:**
  `bio3d::get.pdb("1znb", format="cif")` → `read.cif` (class "pdb", 458 sites for
  the full multi-chain CIF) → `setup_enm` OK, `get_msf_site` works. Caveats
  recorded in NEWS: bio3d flags `read.cif` as beta; it warned helix/sheet records
  unparsed (`set_enm` doesn't need them).
- **Step 3 (active-site vector):** deleted `get_active_site` → `R/pdb_utils.R`
  emptied → `git rm`'d the file. Inlined `c(99,101,103,162,181,184,193,223)`
  (from `data-raw/raw/dataset_1znb_A.csv`) in `data-raw/prepare_znb_data.R`,
  `test-spm-generate.R` (`PDB_SITE_ACTIVE` const), `vignettes/precompute.R` +
  `.Rmd`. Kept `znb_dataset` shipped; rewrote its data-doc as *illustrative*
  (its `@seealso get_active_site` would've broken check). Moved `bio3d`
  Imports→Suggests and dropped `stringr` from Imports + the two `@importFrom`
  lines (both unused in `R/` after the deletions).
- **Step 4 (drop grid):** `git rm` `R/msa_a1a2grid_workflow.R` +
  `tests/testthat/test-a1a2grid.R`. Removed the grid lines from
  `test-contract.R` (pdb_site-on-output already covered by the prediction/
  decomposition checks; pdb_site alignment by the independent manual-join test).
  Replaced vignette §6 with a `preprocess_spm` once → `map_dfr` over
  `calculate_dr2i_msa` example (makes reshape-once explicit). Reworded
  `site_properties.R` `@param site_data`. globalVariables unchanged (`dr2_msa`
  still used in `calculate_loglik_msa`; `dr2_mm/ms/ma` were never listed).
- **Verify:** `document()` clean (14 exports, 4 dropped `.Rd` deleted, no stale
  bio3d/stringr in NAMESPACE); `test()` 45 pass / 0 fail / 0 warn; reran
  `precompute.R` (~100s MCMC) — cache now holds `phi_*` (confirmed no `shap` left
  anywhere in the `.rds`), so the built vignette matches the renamed columns;
  `check()` clean vs baseline. The MCMC is seeded (`precompute.R:20`
  `set.seed(1024)`, base-R RNG), so the rerun is deterministic: every numeric
  value in the cache is bit-for-bit identical to the committed one (verified
  parameter_summary + all.equal(prediction_summary) TRUE); only the `component`
  factor labels changed `shap_*`→`phi_*`. (No drift — an earlier claim of
  "slight drift" was wrong.)

### 2026-06-11 — v0.2 started: API decisions (plan-docs-first)

- **v0.2 entered plan mode; four API decisions made with the user, recorded in
  `dev/plan.md` BEFORE any code** (fix-plan-first). No new model capability;
  exports 18 → 14; all four are breaking; `tmp_src/` not touched.
  1. **`shap_*` → `phi_*`** rename (misnomer; φ in the paper). Load-bearing edits
     confined to `R/msa_decomposition.R` (cols + validation + `starts_with`),
     `test-decomposition.R`, vignette labels/prose, roxygen prose. No `data/*.rda`
     regen (no `shap_` columns there) — BUT verified `vignette_cache.rds` carries
     `shap_act/mut/stab` as *values* in its `component` column, so
     `precompute.R` must be rerun (~100s MCMC; demo R²/summaries will drift
     slightly — cosmetic, expected).
  2. **Structure input = bio3d pdb object, NOT a path** (changed from the
     roadmap's earlier "accept a `.pdb` path"). Rationale (user-driven): the real
     boundary is already `setup_enm(pdb, ...)`; `load_protein` is a one-line
     `read.pdb` wrapper forcing an ID→filename layout. Drop it; user calls
     `bio3d::read.pdb()`/`read.cif()`. Decouples from file I/O and gets **mmCIF
     for free** (bio3d reads mmCIF via `read.cif`; the package no longer sniffs
     extensions/tracks formats). Add `inherits(pdb,"pdb")` guard to `setup_enm`.
     CAVEAT to verify empirically (penm not in repo): does `set_enm` work on a
     `read.cif` object, or need a `.pdb`-only field? If the latter, dial back the
     mmCIF docs claim.
  3. **Active-site input = plain `pdb_site_active` integer vector.** Drop
     `get_active_site` (only thing forcing the `dataset_ec2024.csv` lookup);
     downstream (`generate_spm_data`, `add_site_properties`) already take the
     vector. Inline 1znb_A's vector `c(99,101,103,162,181,184,193,223)` in
     data-raw/tests/vignette. **Keep `znb_dataset`** shipped as *illustrative* of
     the source CSV format (no longer consumed by package code).
  4. **Drop the grid API** (`define_selection_grid` +
     `calculate_dr2i_msa_a1a2grid`). No real consumer — fit + vignette sweeps use
     `calculate_dr2i_msa(pp, ...)` directly; the grid fn also re-ran
     `preprocess_spm` internally (violates reshape-once). A sweep is a 3-line
     `map_dfr` in user code. Delete `R/msa_a1a2grid_workflow.R` +
     `test-a1a2grid.R`; re-prove the `pdb_site`-alignment property in
     `test-contract.R` via `calculate_dr2i_msa` + `site_map` (not the grid route).
- **Also tightened the durable "precomputation" finding in `dev/plan.md`** — it
  conflated three stages the user wants kept distinct: SPM *physics*
  (`generate_spm_data`, the actual precomputation, once) vs *reshaping*
  (`preprocess_spm`, cheap, not a computation) vs *reweighting*
  (`calculate_dr2i_msa`, the only a1/a2-dependent step).
- `dev/PROGRESS.md` rewritten as the v0.2 checklist (Step 0 done). No v0.2 code
  written yet.

### 2026-06-10 — planning restructured + full tmp_src survey

- **Restructured planning into two tiers (user-driven).** Lesson from v0.1: the
  575-line `dev/plan.md` was both blueprint and execution script, so every
  execution surprise forced a plan edit. Now: `dev/plan.md` = COARSE durable
  roadmap (what each version is) + durable findings; detailed planning is done
  per-version at execution time (in plan mode, against the actual code) and tracked
  in `dev/PROGRESS.md`. Rewrote `dev/plan.md` (dropped the 13 v0.1 sections — that
  job is done; history is here), repurposed `dev/PROGRESS.md` as the in-flight-
  version checklist (currently dormant), and rewrote `CLAUDE.md` for the post-v0.1
  reality + the two-tier model.
- **Full read of `tmp_src/` (.R + .Rmd) — much more unmigrated than §12 tracked.**
  The old §12 only listed unmigrated `R/*.R` files. Surveyed all four regions:
  (A) deferred R/ model code [assessment 3-file bundle, allotment, protein decomp,
  viz]; (B) `Rmd/` paper analyses [figure library, AlphaFold2 validation,
  sequence-vs-structure divergence — new feature signals]; (C) `someday_maybe/tree/`
  trajectory route; (D) **`archive/` — the model's MOTION/MODE arm**, which the
  user wants for future projects and which was buried/never analyzed.
  *(2026-06-18: this dir was renamed `.archive/`→`archive/` — unhidden so
  `grep -r` finds it; path updated here to stay valid.)*
- **KEY FINDING (durable, now in plan.md): the precomputation property generalizes
  for the SPM-mean route but NOT the tree.** The model is a 2×2 grid {structure
  dr2, motion dh} × {site i, mode n} (+ nh_n); v0.1 implements only structure×site.
  All archive quantities are `X = sum(pfix_jm * X_jm)/sum(pfix_jm)` with a1/a2 only
  in per-mutant pfix → the v0.1 `[mutant×site]` reweighting matrix generalizes to
  more value-matrices reweighted by the same weights. So motion/mode (D) ≈ "more
  SPM columns + reweighting", mechanically close to v0.1. The TREE (C) uses
  sequential substitutions → weights don't factorize → it interpolates over an
  a1/a2 grid (akima) instead; no clean precomputation. Hence roadmap order
  v0.3=motion/mode (D), v0.4=tree (C).
- **Roadmap set (user decisions):** v0.2 = API-only (phi rename, PDB path input,
  pdb_site-vector active sites, grid keep/demote/drop); v0.3 = motion/mode (D),
  read rest of `archive/` incl. model_rates.R at its start; v0.4 = tree (C);
  later/unordered = assessment+allotment+protein-decomp, viz/figure-lib, AF2,
  sequence divergence. No v0.2 code written yet.

### 2026-06-10 (earlier)
- **§9 step 10 done — `devtools::check()` as clean as v0.1 gets (GitHub-only, not
  CRAN).** Final: **0 errors, 1 warning, 2 notes**, tests pass. Fixed the two
  genuinely-fixable items: (N2) `.Rbuildignore`d top-level `CLAUDE.md`; (N4)
  silenced the package-wide tidyverse-NSE "no visible binding" NOTE via
  `utils::globalVariables(c(...))` in `R/msamodel-package.R` (28 data-column /
  grouping names — not real globals). **Accepted (deliberate, user is GitHub-only):**
  LazyData WARNING + installed-size NOTE (both the known `znb_spm`~13 MB/`znb_wt`~6 MB
  data-size story) and the spurious future-timestamp NOTE (clock/network skew).
  **Dead end:** tried removing `LazyData: true` to clear the WARNING — broke 19
  tests (`object 'znb_spm' not found`: tests use the `znb_*` datasets bare, relying
  on lazy auto-load). Reverted; not worth editing the whole test suite to chase a
  cosmetic warning. **v0.1 migration COMPLETE** — all 10 plan §9 steps done. Next:
  the v0.2 backlog in plan §12 (rename "Shapley"→phi_*, reconsider the grid API,
  the deferred assessment/allotment/protein-decomp/trajectory/viz layers).
- **Vignette revisions (user review of step 9).** Per user feedback on the rendered
  HTML: (3) added active-site vlines to the a1/a2 sweep plots; (4) moved the R²
  annotation to top-left on the obs-vs-fitted profile (was overlapping); added a
  "four model variants" table (MM/MS/MA/MSA) — they were referenced but never
  defined; (5) rewrote the decomposition section — **dropped "Shapley" and the
  wrong "baseline" gloss** (the latter was my invention, grounded in nothing — the
  `shap_mut` term is the mutation/no-selection component, which still varies among
  sites; a baseline would be residue-independent), added the explicit φ formulas
  (φ_mut, φ_stab, φ_act with the ½-averaging) + the additivity identity; (6)
  demoted the grid section to a brief mention. Flagged both v0.2 items in plan §12.
  Re-rendered + re-inspected the changed figures; cache unchanged (only plot/prose
  edited).

### 2026-06-10
- **§9 step 9 done — vignette `msamodel-intro.Rmd` written, builds clean.** Full
  plotted walkthrough (supersedes the old §8 "no plots, print() only" spec —
  rewrote §8 first, fix-plan-first). **Build mode: precompute + cache** —
  `vignettes/precompute.R` runs the heavy pipeline once (SPM preprocess, profiles,
  a1/a2 sweeps, a 4000-iter/1000-burn MCMC fit, decomposition) and saves
  `vignettes/vignette_cache.rds` (96 KB); the `.Rmd` shows the pipeline code as
  `eval=FALSE` and draws every plot off the cache. Keeps check fast, no `penm`/MCMC
  at build time. `precompute.R` is `.Rbuildignore`d; the cache ships. ggplot2 is
  Suggested, so every plot chunk is guarded by `requireNamespace("ggplot2")` +
  `eval = have_cache && has_ggplot`.
  **Sections:** (1) setup from PDB (load_protein→get_active_site→setup_enm→
  generate_spm_data, eval=FALSE, noted as reproducing the shipped znb_wt/znb_spm);
  (2) profile at a demo (a1=2,a2=500) — raw lrmsd vs residue (active sites marked),
  vs lrmsf, vs dactive (+geom_smooth); (3) a1-sweep (a2=0) and a2-sweep (a1=0),
  whole raw lrmsd; (4) Bayesian fit — parameter summary + observed-vs-fitted on
  **mean-centered nlrmsd** (what the likelihood compares) as profile/scatter/vs-
  lrmsf/vs-dactive, each annotated with **R² (≈0.605)** computed inline (v0.1 has
  no gof fns — §1); (5) Shapley decomposition (shap_mut/stab/act, active sites
  marked — activity-selection troughs visibly coincide); (6) grid API note.
  **Decisions (user):** gof = R²; centering differs by section (raw for profile/
  sweep, centered for fit-vs-obs); full PDB setup path shown. **Verify:** ran
  precompute (~99 s, MCMC dominates) → cache written; visually inspected all 10
  figures as standalone PNGs (Read tool) and iterated (added active-site vlines to
  the decomposition); `rmarkdown::render` → 1 MB self-contained HTML, all 46 chunks
  ran; **`devtools::check()` = 0 errors**, 1 warning + 4 notes, **all preexisting /
  step-10** (LazyData + installed-size = the known znb_spm/znb_wt data-size issue;
  top-level CLAUDE.md; future-timestamp; the package-wide tidyverse NSE
  "no visible binding"). Confirmed via git that this work touched only
  `.Rbuildignore`, `dev/plan.md`, `vignettes/` — DESCRIPTION/data untouched, so the
  warning/notes are not vignette-introduced. Only step 10 (`check()` clean) remains;
  the data-size warning/note is the main thing to resolve there.

### 2026-06-09 (later)
- **§9 step 8 done — test suite written (8 files, 55 tests, all pass).** Scaffolded
  `tests/testthat.R` + `tests/testthat/`. Files: spm-generate (the fixture-drift
  guard — reproduces znb_wt/znb_spm from the §5 params, ~20s, dominates runtime),
  spm-preprocess (dr2mat shape + the site_map mapping the pdb_site contract relies
  on), msa-evaluate (dr2i finiteness + a FROZEN loglik literal -184.3241923285),
  msa-mcmc (prior-range validation errors + deterministic fix_a1/fix_a2 constancy
  + shape), decomposition (validation stop()s + the pdb_site optional branch both
  ways + summary shape — additivity NOT asserted, it's tautological), a1a2grid
  (54-row grid + output cols + pdb_site alignment), workflow-endtoend (7-element
  return, real pipeline), contract (unknown-pdb_site error, subset-coverage finite
  loglik, independent-route mapping equality, pdb_site on outputs). Suite: 55 pass
  / 0 fail / 0 warn, ~29s. **Failure-sanity check:** perturbed the frozen loglik
  and the dr2mat shape → the relevant tests FAILED, then reverted → clean. Confirms
  assertions aren't vacuous.
- **Plan §7 reconciled for step 8 (tests).** Updated the test matrix to the
  pdb_site contract and to vetted assertions (the prior matrix had tautological /
  untestable rows). Key changes: dropped the Shapley-additivity assertion (it is
  algebraically the implementation — can't fail except by editing the formula);
  replaced the loglik "regression" with a FROZEN literal (captured once, not
  recomputed in-test); dropped the "acceptance rate > 0" claim (`accepts` isn't
  returned); added deterministic `fix_a1/fix_a2` and validation-error tests for
  MCMC; made decomposition test validation + the pdb_site optional-branch + shape.
  Added **`test-contract.R`** (8th file) for the pdb_site contract: unknown-pdb_site
  error, subset-coverage finite loglik, independent-route mapping equality (non-
  circular), pdb_site on outputs. v0.1 test count 7 → 8. Done before writing tests
  (fix-plan-first). See memory `feedback-test-quality`.

### 2026-06-09
- **Fixed the `observed_data` site-key contract (pdb_site, not i) + documented
  datasets (step 7).** Surfaced while prepping step 7: `observed_data` is
  *user-provided* input (public arg of `run_msa_bayesian_analysis`), but the
  contract required the user to supply `i` — the model's internal 1..N response-
  site renumbering — which a user can't reliably produce, and the fit joined on it
  with no validation (silent-misalignment risk). For 1znb_A it worked by luck
  (verified 225/225 pdb_site agreement), not by construction.
  **Decision (user): full fix overriding CLAUDE.md §4 "rename nothing"** — justified
  as a public-API correctness fix. New contract: `observed_data = {pdb_site,
  lrmsd_obs}`; package maps pdb_site→i internally via a new `site_map` returned by
  `preprocess_spm`, erroring on any unknown pdb_site. All six site-keyed objects
  now also carry `pdb_site` alongside `i` (prediction_samples/summary,
  decomposition_samples/summary, observed_data passthrough, grid output). Internal
  i/j convention (i=response site, j=perturbed site) unchanged. Files:
  `R/msa_bayesian_data_preparation.R` (site_map), `R/msa_model_evaluation.R`
  (map+validate+join), `R/msa_bayesian_analysis.R`, `R/msa_decomposition.R`
  (pdb_site optional — synthetic test input without it still works),
  `R/msa_a1a2grid_workflow.R`, `R/msa_bayesian_workflow.R` (docs).
  **`znb_profile` regenerated** as `{pdb_site, lrmsd_obs}` (dropped i; dropped the
  derivable dactive/lrmsf — verified reproducible from znb_wt to machine
  precision: dactive=get_dactive, lrmsf=log(sqrt(get_msf_site))). SPM one-time
  validation still passed on re-run.
  **Step 7:** wrote `R/data-doc.R` documenting all 5 datasets (`@family datasets`,
  accurate `@format`; znb_profile framed as the {pdb_site, lrmsd_obs} fit target).
  **Verify:** correctness gate — new contract gives loglik identical to the old
  i-join reference for 1znb_A (-184.3242, all.equal TRUE); negative test — unknown
  pdb_site errors clearly; `document()` clean, NAMESPACE still 18 exports;
  `load_all()` clean; full workflow returns 7 elems with pdb_site on every
  site-keyed output. Updated `dev/plan.md` §2/§5 (contract note + deliberate §4
  exception) BEFORE coding. v0.2 maybe: make pdb_site the canonical key throughout.

### 2026-06-06
- **§9 step 6 done — embedded `znb_*` data, SPM generated + validated.** Vendored
  the small inputs (`inst/extdata/1znb_A.pdb`; `data-raw/raw/dataset_1znb_A.csv`
  1 row; `data-raw/raw/profiles_1znb_A.csv` 225 rows), wrote
  `data-raw/prepare_znb_data.R` (generates `znb_spm` via the package's own
  `generate_spm_data`, params as named constants flagged to match the step-8
  test), ran it. **One-time validation PASSED:** `all.equal(znb_spm,
  readRDS(tmp_src/.../1znb_A_spm.rds), tol=1e-8)` is TRUE — the migrated
  `generate_spm_data` reproduces the source SPM (migration-correctness proof;
  block fenced as removable). Saved 5 `.rda` via `use_data(compress="xz")`.
  Sanity: znb_spm 2508 rows, znb_profile 225, 8 active sites, dr2mat 2280×228.
  Added `^data-raw$` to `.Rbuildignore`; wrote `data-raw/source_files.txt`. Scan
  runtime ~20 s. **Flag for step 10:** `znb_spm.rda` ≈ 13 MB and `znb_wt.rda`
  ≈ 6 MB even xz-compressed — will trigger an `R CMD check` data-size NOTE; needs
  a decision later (subsample? fewer mutations? accept the NOTE?). Not a step-6
  blocker.
- **Plan change — §5 embedded-data approach: generate the SPM, don't copy it.**
  Rejected the original §5 (copy `tmp_src/.../1znb_A_spm.rds`, point `SRC` at
  `tmp_src/`/the live project): not self-consistent — once `tmp_src/` is deleted,
  `data-raw/` would have no in-repo recipe to rebuild `data/`, and pointing at the
  live project isn't reproducibility. **New approach:** vendor only the *small*
  inputs (PDB → `inst/extdata/`, the 1znb_A rows of the two CSVs →
  `data-raw/raw/`), and **generate** `znb_spm` with the package's own
  `generate_spm_data()`. Validate the generated SPM **once** against the
  `tmp_src/` original via `stopifnot(all.equal(...))` (migration-correctness
  proof, fenced as removable) — this is NOT the permanent test; the
  regenerate-and-compare drift guard remains step 8's `test-spm-generate.R`.
  **Feasibility verified before deciding:** recovered the generation params from
  `tmp_src/scripts/03_scan_mutants_all_cases.R` (n_mutations=10, model="lfenm",
  sigma=0.3, min_sd=2, seed=1024; ENM node="ca", ming_wall, d_max=10.5,
  frustrated=FALSE); confirmed the scan is seed-deterministic — penm's
  `get_mutant_site_lfenm` does `set.seed(seed + site_mut*mutation)` per mutant
  (tested: same seed → identical mutant xyz, different seed → different;
  penm 0.2.0.9000, R 4.2). The .rds carries no provenance, so the params are
  hard-coded in the script and must match the step-8 test (option a). **Edited
  `dev/plan.md` §5 (rewrite), §2 (layout: znb_spm "generated"; data-raw adds
  `raw/`), §9 step 6, §13 decision 5.** Done as a standalone turn *before* writing
  any step-6 code, so a fresh context reads the corrected plan (per
  `CLAUDE.md` "read plan.md before doing anything"). Step-6 code execution is the
  next turn.

### 2026-06-05
- **§9 step 4 done (all sub-steps 4.1–4.7), in one batch + step 5.** Copied the 7
  v0.1 source files from `tmp_src/R/` into `R/` keeping source filenames (per the
  no-rename change below), in dependency order. Applied §4 rules: stripped all
  `@requires` tags; added `@export` to the 18 public fns + `@family` groups; the
  2 helpers (`delta_structure_dr`, `delta_structure_dr2`) marked `@noRd`. No
  inline `library()` existed to strip. **Bug fixes:** §11.1 — deleted the
  duplicate `define_selection_grid` (the `a1_range, a2_range` override at source
  lines 106–119), keeping the real defaulted one; §11.3 — fixed the
  `"lower,Quite upper"` roxygen typo in `calculate_prediction_summary`'s
  `@return`. **4.7 assessment removal:** dropped the two
  `calculate_model_comparison_samples/summary` blocks + their return-list entries
  from `run_msa_bayesian_analysis` → 7-element return; updated its `@return`.
  **Imports:** added `transmute` + `group_split` to the dplyr `@importFrom` in
  `msamodel-package.R` (used by `msa_decomposition.R`). **Verify:**
  `devtools::document()` clean (no roxygen warnings; 18 `.Rd` written, none for
  the `@noRd` helpers); NAMESPACE has exactly 18 exports, `delta_structure_*`
  absent; `devtools::load_all()` clean. Behavioral/numeric checks deferred to
  steps 6–8 (no data/tests yet).
- **Plan change — no file renames in v0.1.** Dropped the §2 rename list (10
  renames) before starting §9 step 4. Files are now copied from `tmp_src/R/`
  keeping their **source filenames**. Rationale: source names keep provenance
  obvious during the migration, and some rename rationales were v0.1-scoped (e.g.
  `msa_bayesian_workflow.R` → `msa_workflow.R` "because only one workflow in
  v0.1") — an assumption v0.2 breaks, which would force a *second* rename (churn
  the migration explicitly avoids). File naming/layout is now an explicit later
  refactor, not folded into the copy. Edited `dev/plan.md` §2 (rename block + `R/`
  layout tree) and §9 step-4 sub-steps; regenerated the §9-step-4 labels in
  `dev/PROGRESS.md` to match. Deferred-file renames mentioned elsewhere in the
  plan are left as advisory (decided when those files are actually migrated).

### 2026-06-04
- **Docs split:** `dev/PROGRESS.md` was doing two jobs (append-only log + a
  checklist meant to mirror plan.md §9); the re-scope made the checklist
  numbering diverge from the plan and lose meaning. Split into three files:
  `dev/plan.md` (normative, unchanged role; §9 step-4 re-labelled `4.1`–`4.7`),
  `dev/PROGRESS.md` (pure checklist, bare 1:1 mirror of §9), and this
  `dev/LOG.md` (history). Updated the `## Tracking progress` convention in
  `CLAUDE.md` and the cross-references.
- **v0.1 RE-SCOPED (major):** narrowed v0.1 to "compute divergence profiles + fit
  to data". Traced the real call graph from `run_msa_bayesian_analysis` + the grid
  entry points to decide scope. **Kept:** setup, SPM, profiles/eval, grid, MCMC
  fitting + the four model profiles, site-level Shapley decomposition (wired into
  the workflow). **Deferred to vNext:** the model-fit *assessment* layer
  (`model_comparison_functions.R` + its sole-consumer `loess_compare.R` + the
  `utils.R` error metrics), *allotment* (unsettled), protein-level decomposition,
  trajectory route, visualization. Decisions: keep site decomposition; defer
  assessment (NOT the profiles — those are kept); pull `utils.R` + `loess_compare.R`
  back out. **Reverted** R/utils.R + R/loess_compare.R + their man pages + 7
  exports (forward commit, no history rewrite); `document()`/`load_all()` clean,
  NAMESPACE now 0 exports, R/ = just `msamodel-package.R`. **Rewrote `dev/plan.md`
  §1/§2/§6/§7/§8/§9/§11/§12** to match (exports 39 → 18; tests 9 → 7; workflow
  edited to drop assessment calls → 7-elem return). Reconciled `CLAUDE.md` and
  the progress docs. The §9 step-4 sub-step list is the new dependency order.
  Next: resume migration at §9 step 4.1 (`pdb_utils.R`, `enm_setup.R`,
  `site_properties.R`) against the re-scoped plan.
- **§9 step 4.2 done** (later REVERTED in the re-scope above): copied
  `compare_loess_fits.R` → `R/loess_compare.R` (renamed per §2).
  `compare_loess_fits` + `compare_loess_rmse` exported with
  `@family loess comparison`; `plot_loess_comparison` kept `@noRd` (deferred to
  v0.3), inline `library(ggplot2)` stripped, every ggplot2 call qualified
  `ggplot2::` behind the existing `requireNamespace` guard (§4 rule 2/8).
  Dropped the orphaned `@param evaluation_points` (not a real arg — doc-only).
  `gridExtra` left out of Suggests and no `globalVariables()` shim added (both
  decided: minimal change, the only consumer is the deferred internal plot fn).
  `document()` wrote 2 `.Rd` (none for the `@noRd` fn) + NAMESPACE (then 7
  exports); `load_all()` clean; smoke test of both exported fns returned finite
  values with the expected result names.
- **Toolchain decision (not a §9 step):** investigated whether to update the R
  stack first. Finding: R itself is old (4.2.1; latest 4.6.0) but the packages
  are current (dplyr 1.1.4, ggplot2 4.0.1, etc.). **Decision: do NOT upgrade R
  now** — 4.2.1 satisfies the plan (`R >= 3.5`); revisit after v0.1. When done,
  upgrade R **side-by-side via `rig`** (keep 4.2 for the MSA paper project),
  reinstall packages into the new library, rebuild `penm` from source — that's
  the safe path (single shared library; all 436 pkgs built under 4.2.x). The
  only risk to penm/paper is upgrading *R*, not packages.
  **Updated roxygen2 7.2.3 → 8.0.0** (compiles under R 4.2; `Depends: R >= 4.1`):
  resolves the version-mismatch warning. 8.0.0 records its version as the new
  `Config/roxygen2/version: 8.0.0` field and removed the old `RoxygenNote` line.
  8.0.0 also enforces **one-line `@importFrom`** — flattened the wrapped dplyr/
  stats/penm directives in `msamodel-package.R`. `document()` + `load_all()`
  clean; 5 utils exports intact (those exports later reverted in the re-scope).
- **§9 step 4.1 done** (later REVERTED in the re-scope above): copied `utils.R`
  (`rmse_trend`, `rmse`, `r2`, `mrr`, `mrr_trend`). Stripped the `@requires`
  tags, added `@export` + `@family error metrics` to all 5, kept
  signatures/logic identical (no `library()` calls in source). `document()`
  wrote 5 `.Rd` + NAMESPACE (5 exports); `load_all()` clean. Bare
  `loess`/`predict`/`cor` resolve via the existing `@importFrom stats`.
- **§9 step 3 done:** added `R/msamodel-package.R` with the package-wide
  `@importFrom` directives (Option A namespacing — imports declared once;
  copied function files use bare calls). Inventory verified against the 16
  migrated source files; ggplot2 excluded (Suggests). `document()` populated
  NAMESPACE (incl. `importFrom(rlang,"!!")`, `importFrom(magrittr,"%>%")`, no
  exports yet); `load_all()` clean.
- Set up source isolation: froze a snapshot of the source project to
  `tmp_src/` (build- and git-ignored), removed the live source from the
  session's additional directories, repointed all docs at `tmp_src/`.
- Wrote `CLAUDE.md`.
- Reviewed `dev/plan.md` against the actual source and reconciled it: fixed the
  `load_protein` test, made `plot_loess_comparison` internal (ggplot/Suggests),
  added a rule to strip non-standard `@requires` roxygen tags, corrected the
  `library()` inventory, added the `R/data-doc.R` step, folded the resolved
  pre-flight questions into §13.
- Reorganized: `dev/plan.md` is now the single normative spec; deleted
  `dev/RESUME.md`; added `dev/PROGRESS.md`.
- `git init` on `main`. **Initial commit done** (`020ac2b`): skeleton +
  `dev/plan.md`, `dev/PROGRESS.md`, `CLAUDE.md`. `tmp_src/` and
  `.claude/settings.local.json` excluded.
- **Remote set up:** public repo `github.com/jechave/msamodel` (SSH origin),
  `main` pushed and tracking `origin/main`. Verified `tmp_src/` not pushed.

### 2026-04-22
- Created the package skeleton: `DESCRIPTION` (plan §3), `LICENSE` + `LICENSE.md`
  (MIT), empty roxygen-managed `NAMESPACE`, `.Rbuildignore`. `R/` left empty.
- Moved `plan.md` into `dev/`.
- Source-verification pass; answered the pre-flight questions.
