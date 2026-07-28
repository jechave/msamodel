#!/usr/bin/env Rscript
# dev/callgraph.R — regenerate the package's function call graph.
#
#   Rscript dev/callgraph.R
#
# Produces (git-ignored, in dev/preview/):
#   callgraph.pdf   — a 2-page vector document, one arm per page, each sized to
#                     print legibly on a single Letter sheet (portrait):
#                       p1  site branch   (site fns + everything they call)
#                       p2  mode branch   (mode fns + everything they call)
#
# HOW IT WORKS (no hand-drawn edges): mvbutils::foodweb does static analysis of the
# loaded package's function bodies to get the real caller -> callee edges; those are
# written to Graphviz DOT and rendered with `dot`. Every arrow means "A calls B". Each
# page is an induced subgraph: an arm's own functions plus their transitive callees (so
# the shared spine each arm reaches is shown, coloured grey). Re-run after any change to
# R/ to refresh the picture.
#
# PRINT-LEGIBLE by design: high-contrast (black text, black arrows, saturated fills),
# and each page is scaled UP to fill a Letter sheet at a readable font. The whole-graph
# view was dropped — at 47 nodes it cannot be both on one sheet and legible; the two arm
# pages together already cover every function.
#
# Requires: mvbutils + qpdf (R pkgs) and graphviz's `dot` on PATH
# (`brew install graphviz`).

suppressMessages({
  library(mvbutils)
})

# --- locate the package root (works whether run from root or dev/) -----------------
here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
root <- if (!is.null(here)) normalizePath(file.path(here, "..")) else getwd()
if (!file.exists(file.path(root, "DESCRIPTION"))) root <- getwd()

# --- fail loud if a dependency is missing (don't silently produce nothing) ----------
if (nzchar(Sys.which("dot")) == FALSE) {
  stop("graphviz `dot` not found on PATH. Install it first: brew install graphviz")
}
if (!requireNamespace("qpdf", quietly = TRUE)) {
  stop("R package `qpdf` not installed (needed to merge the pages). install.packages('qpdf')")
}

# --- load the package under development ---------------------------------------------
suppressMessages(devtools::load_all(root, quiet = TRUE))
ns   <- getNamespace("msamodel")
objs <- ls(ns)
fns  <- objs[vapply(objs, function(o) is.function(get(o, envir = ns)), logical(1))]

# --- static call analysis -> caller x callee matrix (m[a, b] == 1 means a calls b) --
# NB: mvbutils::foodweb opens a graphics device even with plotting = FALSE, which
# dumps a stray Rplots.pdf in the working dir. Route it to a temp file and clean up.
.tmpdev <- tempfile(fileext = ".pdf")
grDevices::pdf(.tmpdev)
fw <- foodweb(where = ns, funs = fns, plotting = FALSE)
grDevices::dev.off()
unlink(.tmpdev)
m  <- fw$funmat
nm <- rownames(m)

# --- classify each node by axis; export status decides the border -------------------
exported <- getNamespaceExports("msamodel")
axis <- ifelse(grepl("_i(_|$)", nm), "site",
        ifelse(grepl("_n(_|$)", nm), "mode", "agnostic"))
axis[nm == "preprocess_spm_mode"] <- "mode"   # named mode, no _n token

# --- high-contrast, print-first palette --------------------------------------------
# Saturated arm fills, BLACK node text and BLACK arrows so it survives a b/w laser and
# holds contrast in colour. Shared (non-arm) nodes are a solid mid-grey on branch pages.
# dashed border = @noRd internal; solid = exported.
ARM <- list(
  site = list(fill = "#bfe0d8", border = "#0f4b41"),   # teal
  mode = list(fill = "#f2dca6", border = "#7a4e00")     # gold
)
GREY <- list(fill = "#d9d7d2", border = "#3a3f4b")      # shared spine

# One node-declaration line. `grey` overrides the arm colour for shared nodes.
node_line <- function(x, arm_axis, grey = FALSE, indent = "  ") {
  is_exp <- x %in% exported
  col    <- if (grey) GREY else ARM[[arm_axis]]
  sprintf(paste0("%s\"%s\" [shape=box, style=\"filled,%s\", color=\"%s\", penwidth=1.6, ",
                 "fillcolor=\"%s\", fontcolor=\"black\"];"),
          indent, x, if (is_exp) "solid" else "dashed", col$border, col$fill)
}

# Transitive closure: all functions reachable from `roots` via the call matrix.
reachable <- function(roots) {
  seen <- character(0)
  frontier <- intersect(roots, nm)
  while (length(frontier)) {
    seen <- union(seen, frontier)
    nxt  <- character(0)
    for (a in frontier) nxt <- union(nxt, colnames(m)[m[a, ] > 0])
    frontier <- setdiff(nxt, seen)
  }
  seen
}

# Induced-subgraph DOT for one arm, sized to fill a Letter portrait sheet.
# `roots` = the arm's own functions; nodes NOT in `roots` are shared spine, drawn grey.
arm_dot <- function(roots, title, arm_axis) {
  keep <- reachable(roots)
  L <- c(
    "digraph arm {",
    # These arm graphs are inherently WIDE and SHALLOW (many roots across the top, only
    # a few ranks deep), so left alone dot renders a thin horizontal strip that scales
    # to a tiny band on the sheet. ratio="fill" with an explicit size FORCES the drawing
    # to the box's proportions -- it stretches the few ranks apart vertically to fill a
    # landscape-Letter sheet, so nodes and text are large on paper. nodesep/ranksep set
    # generous minimums so the fill has room to work.
    "  size=\"10,7.5\"; ratio=fill; margin=0;",
    "  labelloc=\"t\";",
    sprintf("  label=\"%s\";", title),
    "  fontname=\"Helvetica-Bold\"; fontsize=16;",
    "  rankdir=TB;",
    "  graph [splines=true, nodesep=0.5, ranksep=1.0];",
    "  node  [fontname=\"Courier-Bold\", fontsize=15, margin=\"0.18,0.10\"];",
    "  edge  [color=\"black\", penwidth=1.2, arrowsize=0.9];"
  )
  for (x in keep) L <- c(L, node_line(x, arm_axis, grey = !(x %in% roots)))
  for (a in keep) for (b in intersect(colnames(m)[m[a, ] > 0], keep)) {
    L <- c(L, sprintf("  \"%s\" -> \"%s\";", a, b))
  }
  c(L, "}")
}

site_roots <- nm[axis == "site" | nm == "preprocess_spm"]
mode_roots <- nm[axis == "mode"]

# --- render each page to its own temp PDF, then merge --------------------------------
outdir <- file.path(root, "dev", "preview")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

render_pdf <- function(dot_lines) {
  dotfile <- tempfile(fileext = ".dot")
  pdffile <- tempfile(fileext = ".pdf")
  writeLines(dot_lines, dotfile)
  # -Gsize is already in the DOT; force Letter media box so the sheet is standard.
  ok <- system2("dot", c("-Tpdf", shQuote(dotfile), "-o", shQuote(pdffile)))
  unlink(dotfile)
  if (ok != 0) stop("dot failed to render a page")
  pdffile
}

pages <- c(
  render_pdf(arm_dot(site_roots,
                     "msamodel — site arm (i) + everything it calls   (arrow: A calls B; teal = site, grey = shared; dashed = @noRd)",
                     "site")),
  render_pdf(arm_dot(mode_roots,
                     "msamodel — mode arm (n) + everything it calls   (arrow: A calls B; gold = mode, grey = shared; dashed = @noRd)",
                     "mode"))
)
on.exit(unlink(pages), add = TRUE)

pdf_out <- file.path(outdir, "callgraph.pdf")
qpdf::pdf_combine(pages, output = pdf_out)

cat(sprintf("call graph: %d functions, %d edges — 2-page print-legible PDF\n  %s\n",
            length(nm), sum(m > 0), pdf_out))
