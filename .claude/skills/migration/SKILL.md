---
name: migration
description: Procedure and provenance discipline for migrating code from tmp_src/ into R/. Invoke when resuming migration of a specific function or property from the frozen source snapshot (motion/mode arm, rates, trees, or any unmigrated capability).
---

# Migrating code from `tmp_src/` into the package

Most planned migration is done; new work is developed directly. Use this skill only
when the user says "we're migrating X" — a specific function or property is being
brought over from the frozen source snapshot.

## The source and the boundary

- Migration sources are **ONLY `tmp_src/`** (a read-only frozen snapshot of the
  original project, `.Rbuildignore`d). **Never edit inside `tmp_src/`** — copy out of
  it only.
- Much unmigrated code lives in **`tmp_src/archive/`** (motion/mode arm, rates, etc.).
  A plain `grep -r` silently skips hidden dirs and once produced a false "no source
  exists" conclusion — always search *all* of `tmp_src/`.
- **penm is a DEPENDENCY, never a migration source.** Call `penm::fn()`. Never copy or
  "migrate" code out of penm — not the installed library, not the local `../penm`
  source. Reading penm to check a *signature you are calling* is fine; reading it to
  source/justify package code is not. If a computation seems to need penm internals,
  that means a `penm::` call, not a copy.

## Step 0 — find the source (canonical check)

Before claiming any function or logic is "new" / "package-native" / "written from
scratch", search for it:

```sh
dev/find-source.sh '<fn name>'        # e.g. calculate_dr2n_msa
dev/find-source.sh '<core formula>'   # e.g. 'sum\(pfix_jm \* dr2_njm\)'
```

A zero-hit search on the **name** is not proof of no source — names get changed in
migration; the **math does not**. Search the formula too. The script searches all of
`tmp_src/` including hidden dirs (`.R`/`.Rmd`) and reports whether a source exists.

## Copy rules (this is a migration, not a refactor — stay close to the source)

When copying a function from `tmp_src/R/` into `R/`:

1. Remove `library(...)` calls (top-level and inline). All package usage goes through
   `Imports:` + `@importFrom` (or `pkg::fn()`).
2. **Keep function signatures identical.** Do not rename functions or variables or
   "improve" logic as part of a copy. (Deliberate, recorded API changes — like the
   v0.2 `phi_*` rename — are separate, planned work, not folded into a copy.)
3. Namespacing: **Option A** by default — package-level `@importFrom` + bare calls
   (smallest diff). **Exception:** a Suggested package (e.g. ggplot2) cannot be
   `@importFrom`'d — qualify as `ggplot2::...` behind a `requireNamespace` guard.
4. Keep `%>%` (magrittr, re-exported).
5. Strip the non-standard `@requires` roxygen tags (every source file has them;
   roxygen2 doesn't recognise them). Keep standard tags.
6. Add `@export` to public top-level functions; helpers stay unexported (`#' @noRd`).
7. Fix catalogued bugs as you reach each file — record them in the version's plan
   (e.g. the tree route's `p_act`→`p_ma` bug, the `akima`→`interp` swap; see
   `dev/plan.md` findings).

**Place the migrated function in the `R/` file for its `@family`** (spm / model /
objective / fitting / decomposition / setup), regardless of which `tmp_src/` file it
came from. Traceability is provided by `dev/find-source.sh` (content search), not by
filename. This does not loosen rule 2 — names and signatures still don't change; only
file placement does.

## Restructure + VERIFY — not clone-and-rename

Migration = restructure + verify. When migrating a property the user has NOT already
refactored, **do not assume axis-swapping an existing function (e.g. site → mode) is
equivalent.** Reproduce the OLD (archive) version's numbers and assert equality —
machine precision for deterministic code, seeded for stochastic — **before** trusting
the new form. "Looks structurally identical" is not verification; run the comparison.
(This verify step is a one-time migration artifact → it belongs in a scratchpad
script, not the permanent suite. See the `/test-review` skill.)

## Provenance honesty (the global honesty rule, applied to code origin)

Never state where code came from ("I wrote it" / "it's migrated from X" / "no source
exists") unless you have *just* searched/read to confirm it. These are factual claims,
governed by the same rule as "verified"/"fixed": no claim without a check. A bad
"migration" of the mode functions followed by repeatedly asserting false provenance is
why this rule exists.
