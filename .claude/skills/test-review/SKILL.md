---
name: test-review
description: Checklist and mechanical gate for writing or reviewing R package tests — reject tautologies and circular re-derivations, force a negative control before trusting a test, and decide permanent-vs-disposable. Invoke whenever writing, adding, or reviewing a test in tests/testthat/.
---

# Test review — write tests that can FAIL, not tests that pass by construction

The root failure this skill prevents: **writing the assertion to PASS, not to FAIL** —
optimizing for the green checkmark instead of for catching a bug. Every weak-test
symptom (tautology, circular re-derivation, permanent refactor-invariance check,
full-suite reflex) flows from that. The fix is mechanical steps that produce output
the green checkmark cannot fake.

## 1. Permanence filter — ask BEFORE adding to `tests/`

Three questions, all must pass:

- **Can it fail for a real reason?** (If not → §2 will catch it; drop or downgrade.)
- **One-time or forever?** A **refactor-invariance / "the numbers didn't change" check
  is one-time** — it belongs in the **scratchpad verify script**, NEVER the suite. A
  test earns a `tests/` slot only if it guards a *future* regression risk.
- **Already covered?** If an existing test (positivity, an independent-route nested
  check, a frozen-loglik snapshot) already guards this, adding another is churn, not
  rigor. More green tests ≠ more rigor.

## 2. Negative control — REQUIRED, and it produces output

Before trusting any new invariant/value test, **feed a deliberately-wrong input**
(swapped `(a1,a2)`, a perturbed value, a shuffled vector) and confirm the assertion
goes **RED**. Do this **before**, not after.

- If it can't go red, it is a tautology — delete it or downgrade to a smoke test.
- This is the falsifiable step: run it, see the red, then revert to the correct input
  and see the green. The red output is the artifact.

## 3. Anti-patterns — reject on sight

- **True by construction / re-derives the implementation.** e.g. asserting decomposition
  additivity `phi_mut + phi_stab + phi_act == lrmsd_msa` when that identity IS the
  formula — cannot fail except by editing the formula. A change-detector, not a
  correctness check.
- **Circular reference.** The "expected" side must NOT call the function under test (or
  a shared helper the function also calls). If both sides call `weights_jm_spm(...)`,
  garbage weights match garbage weights. For a refactor-invariance check the correct
  reference is the **OLD formula recomputed independently** — never the new fn.
- **Recomputed regression value.** A fixed-value/regression assertion is meaningful
  only if the expected is a **FROZEN literal or a `testthat` snapshot** — not
  recomputed by the code under test at test time.

Identity/algebraic facts evident from the code add ~no value: drop them, or keep as a
labelled smoke test at most.

## 4. Loop discipline

- While iterating: `testthat::test_file(...)` / `devtools::test(filter=)` — the ONE
  relevant target. Stop there.
- The full `devtools::test()` is a **gate action fired once** before a
  code/data/roxygen commit (per the diff-rule in CLAUDE.md) — **never** a mid-work
  "did it pass" reflex after a green `test_file`.
- If you notice the impulse "run all tests to be sure," that impulse itself is the
  signal you are optimizing for the appearance of done. Targeted check, then continue.

## Prefer testing (things that can actually regress)

Real branching logic (validation `stop()`s, optional-column branches), integration /
wiring (does the real pipeline output flow through), error paths (`expect_error`), and
behavior that could plausibly break — over identities the code guarantees.

## Optional: adversarial second pass

For a batch of new tests, have a fresh-context sub-agent read only the new test diff and
flag any tautology/circularity/recomputed-value it sees. A reviewer in fresh context
isn't biased by the reasoning that produced the tests.
