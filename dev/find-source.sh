#!/bin/sh
# find-source.sh — "does this have a migration source in tmp_src/?"
#
# THE canonical check before claiming any function/logic is new or
# package-native. Searches ALL of tmp_src/ recursively, INCLUDING hidden
# directories — which a plain `grep -r` silently skips. (On 2026-06-18 a plain
# grep missed tmp_src/.archive/ entirely and produced a false "no source"
# conclusion; this script exists so the answer never depends on remembering to
# unhide / --include.)
#
# Usage:
#   dev/find-source.sh '<pattern>'      # ERE pattern (function name OR a core formula)
#
# Search BOTH the name and the math: names get renamed during migration, the
# formula does not. Examples:
#   dev/find-source.sh 'calculate_dr2n_msa'
#   dev/find-source.sh 'sum\(pfix_jm \* dr2_njm\)'
#
# Exit status: 0 if matches found, 1 if none.

if [ "$#" -ne 1 ] || [ -z "$1" ]; then
  echo "usage: dev/find-source.sh '<pattern>'" >&2
  exit 2
fi
pattern=$1

# Resolve tmp_src relative to the script's package root (works from anywhere).
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
src="$root/tmp_src"

if [ ! -d "$src" ]; then
  echo "tmp_src/ not found at $src (snapshot removed?)." >&2
  exit 2
fi

# `find ... -name '*.R*'` descends into HIDDEN dirs (.archive etc.); the default
# grep -r does not. Match .R and .Rmd. Use grep -nH for file:line output.
matches=$(find "$src" \( -name '*.R' -o -name '*.Rmd' \) -type f \
  -exec grep -nHE -- "$pattern" {} + 2>/dev/null)

if [ -n "$matches" ]; then
  # Strip the long absolute root prefix for readable, clickable paths.
  printf '%s\n' "$matches" | sed "s|$root/||"
  n=$(printf '%s\n' "$matches" | wc -l | tr -d ' ')
  files=$(printf '%s\n' "$matches" | cut -d: -f1 | sort -u | wc -l | tr -d ' ')
  echo "---"
  echo "FOUND: $n line(s) in $files file(s) under tmp_src/ — a migration source EXISTS."
  exit 0
else
  echo "NO migration source found in tmp_src/ for: $pattern"
  echo "(Searched all of tmp_src/ incl. hidden dirs. Try the core FORMULA, not just the name.)"
  exit 1
fi
