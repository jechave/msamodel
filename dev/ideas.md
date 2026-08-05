# Ideas / maybe later

A parking lot for things worth doing eventually. **Not a plan and not a status
file** — nothing here is scheduled, ordered, or claimed to be current. Not read at
session start; open it when you want to pick something up.

Delete an entry when it is done, or when it stops being a good idea.

---

## Terse band-machinery names

`spm_hmat`, `grad_theta`, `var_param_delta`, `as_theta` (all in
`R/predict_band_helpers.R`, some also used in `R/api.R`) are named in a terse style
the author rejects — as are the internal `M`/`G`/`Gc`/`h`/`w` inside `spm_hmat`.

The `t` → `theta` coordinate piece of this was done (2026-08-04); the rest was not.

Coupled to the notation in `dev/reports/spm_band_validation.Rmd`, which calls
`spm_hmat` by name in 7 places under `load_all()` and mirrors the terse notation
throughout. Rename the code first so the validated behaviour is preserved, then
follow the PDF notation to match.

## "band" → "se" vocabulary

The uncertainty vocabulary is inconsistent: `R/` carries ~28 uses of
band/bands/banded/banding, while the public API returns `_se` columns. Worth
settling on one word. Touches `R/predict_band_helpers.R` (including its filename),
roxygen, and the validation report.

## Rethink the `R/` family / file structure

The `CLAUDE.md` rule fixing how `R/` files map to `@family` tags was **removed
2026-08-05** pending a rethink, deliberately not replaced with a corrected version.

Two concrete inputs for that thinking:

- `api` had become a seventh family the rule never listed (`R/api.R:90,133,185,238`)
  — and it holds the main public surface.
- `R/msa_decomposition.R` carries no `@family` tag at all (`decompose_nested`,
  line 28).

Fixing the missing tag is an `R/` + roxygen change, so it needs `document()` and the
full-suite gate; it was deliberately left out of the docs-only cleanup.

## `predict_decomposition` errors on its own default

`predict_decomposition(fit, spm)` cannot succeed as written: `match.arg` defaults
`which` to `"lrmsd"`, and that branch is an explicit `stop()` saying the uncentred
component band is not yet derived (`R/api.R:246-251`). Either derive the uncentred
component SE, or make `"nlrmsd"` the default so a bare call works.

## `dev/tmp2/` scratch

~14 MB of untracked throwaway scratch. Delete when convenient.

## v0.4 — motion arm

The motion/fluctuation quantities (`dh_ijm`, `dh_njm`, `nh_njm`) are unimplemented;
see `dev/findings.md` for the grid and why the existing reweighting machinery should
extend to them. Available in penm's slow per-mutant loop only. Not scoped.
