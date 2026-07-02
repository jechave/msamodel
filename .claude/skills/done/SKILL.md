---
name: done
description: Definition-of-Done reconciliation gate for a finished work item — a run-commands-and-paste-output check that satellite state (LOG, memory, cross-references, git) is current, not just that tests pass. Invoke at the work-item milestone before reporting done.
---

# Definition of Done — a COMMAND-OUTPUT gate, not a self-attestation

Verifying the *deliverable* (tests pass, NAMESPACE clean) is not enough. A change also
has **satellite state** — the `dev/LOG.md` NOW block, memory files, cross-references —
that goes stale silently. A failing test stops you cold; a stale memory file does not,
it just sits there looking fine. So the "looks done" signal arrives at `git push` with
the easy part done (the LOG entry) and the rest skipped.

**The fix: this gate is RUN COMMANDS + PASTE OUTPUT, not "I reconciled."** "I checked
for stale refs" is unfalsifiable; a grep that prints the dangling line is not. Stale
state must *show up in output*, not be reasoned away.

**When it runs:** at the **work-item milestone**, before reporting done / pushing. Run
against the whole project, not just edited files.

```bash
# 1. STALE CROSS-REFS. For EVERY name you renamed/closed/moved this work item (old
#    file, function, version label, status word like "in flight"/"next"), grep it
#    across code AND satellite state. Substitute the real old term(s):
git grep -n 'OLD_NAME'                         # tracked: R/ tests/ dev/ man/ vignettes DESCRIPTION
grep -rn 'OLD_NAME' ~/.claude/projects/-Users-julianechave-Library-Mobile-Documents-com-apple-CloudDocs-lab-Rpackages-msamodel/memory/
#    -> every hit is either correct or a fix you make NOW. Zero unexplained hits to pass.

# 2. LOG + NOW reflect reality (read, don't assert):
#    open dev/LOG.md, read the <!-- NOW -->…<!-- /NOW --> block: is the live item
#    right? closed work crossed off? + a dated dev/LOG.md history entry for what just
#    happened exists.

# 3. MEMORY currency — the next-session pointer is the one loaded first each session:
#    read the first lines of the memory MEMORY.md index — is that pointer current?
#    + grep (step 1) already proved no memory file contradicts the change.

# 4. GIT clean + pushed:
git status -sb        # clean tree, '## main...origin/main' (no ahead/behind), nothing stray
```

(Use `Read` for steps 2–3 — read the LOG NOW block and the MEMORY.md head directly;
avoid `sed`/`cat` per the tool-choice convention.)

**Pass condition:** steps 1–4 have each been *run this session* and their output is
clean. If any turns up work, it is **part of this task, not a later cleanup** — do it
now, then re-run. Do not report done on "tests pass + committed + pushed" alone — that
is the easy part; the grep output is the gate.
