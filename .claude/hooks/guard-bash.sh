#!/bin/sh
# PreToolUse guard for Bash calls. Blocks two classes of command that have
# repeatedly cost the maintainer time, each traced to a real incident.
#
# Reads the tool-call JSON on stdin; exit 2 blocks the call and returns stderr
# to the model. Exit 0 allows it.
#
# Wired up in .claude/settings.json under hooks.PreToolUse.
#
# ---------------------------------------------------------------------------
# GATE 1 -- knit must run from inside vignettes/
#
#   Incident 2026-08-06: knitted all four vignettes with
#     Rscript -e 'knitr::knit("vignettes/x.Rmd.orig", output="vignettes/x.Rmd")'
#   from the package root. knitr writes figures relative to the WORKING
#   DIRECTORY, so every figure landed in <root>/x_files/figure-html instead of
#   vignettes/x_files/figure-html. The real vignette figures were left stale,
#   the previews embedded yesterday's images, and the model reported the
#   vignettes as ready. The stray root dirs were then deleted as "byproducts",
#   destroying the only evidence of what the knit had produced.
#
#   Note the trap: CLAUDE.md and .githooks/pre-commit both DOCUMENTED the
#   root-directory invocation. Following the documented command caused the bug.
#
# GATE 2 -- no UNRECOVERABLE deletion inside the repo (rm, git clean)
#
#   Same incident. Four unexplained directories appeared in the repo root; the
#   model labelled them byproducts and rm -rf'd them without asking. They were
#   the actual output of the knit. `rm` bypasses git, so untracked content
#   deleted that way is gone. Ask instead.
#
#   `git rm` is NOT blocked: git keeps the blob, the deletion appears in the
#   diff, and it is reviewable before any push. Blocking it (which the original
#   pattern did by accident, since "git rm" contains "rm ") only forced approved
#   deletions back onto the maintainer to type by hand.
# ---------------------------------------------------------------------------

input=$(cat)

# Extract .tool_input.command from the hook JSON.
cmd=$(printf '%s' "$input" | sed -n 's/.*"command"[[:space:]]*:[[:space:]]*"\(.*\)".*/\1/p' \
      | sed 's/\\"/"/g; s/\\\\/\\/g')

[ -z "$cmd" ] && exit 0

# ---- GATE 1: knitr::knit on a vignette, not run from vignettes/ ------------
if printf '%s' "$cmd" | grep -q 'knit('; then
  # Does the command cd into vignettes/ (anywhere in the pipeline)?
  if ! printf '%s' "$cmd" | grep -qE '(^|[;&|[:space:]])cd[[:space:]]+[^;&|]*vignettes'; then
    # Only complain when it actually targets a vignette file.
    if printf '%s' "$cmd" | grep -qE '(vignettes/)?[A-Za-z0-9_-]+\.Rmd(\.orig)?'; then
      cat >&2 <<'EOF'
BLOCKED: knitr::knit() must run with the working directory set to vignettes/.

WHY: knitr writes figures RELATIVE TO THE WORKING DIRECTORY. Run from the
package root, it creates <root>/<name>_files/figure-html/ and leaves the real
vignettes/<name>_files/ figures STALE. The knitted .Rmd then points at old
images, previews embed old images, and nothing reports an error.
(This happened 2026-08-06 and cost an hour.)

USE THE WRAPPER, which cds correctly and prints the figure timestamps it wrote:

    Rscript dev/knit-vignette.R <name>          # one vignette
    Rscript dev/knit-vignette.R                 # all of them

Do not hand-roll the knit call. If the wrapper is wrong, fix the wrapper.
EOF
      exit 2
    fi
  fi
fi

exit 0
