#!/usr/bin/env Rscript
# dev/callgraph.R — regenerate the package's function call graph.
#
#   Rscript dev/callgraph.R
#
# Produces (git-ignored, in dev/preview/):
#   callgraph.png   — open/embed anywhere
#   callgraph.pdf   — vector, for printing
#
# HOW IT WORKS (no hand-drawn edges): mvbutils::foodweb does static analysis of the
# loaded package's function bodies to get the real caller -> callee edges; those are
# written to a temporary Graphviz DOT file and rendered with `dot`. Every arrow means
# "A calls B". Re-run after any change to R/ to refresh the picture.
#
# Requires: mvbutils (R pkg) and graphviz's `dot` on PATH (`brew install graphviz`).

suppressMessages({
  library(mvbutils)
})

# --- locate the package root (works whether run from root or dev/) -----------------
here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
root <- if (!is.null(here)) normalizePath(file.path(here, "..")) else getwd()
if (!file.exists(file.path(root, "DESCRIPTION"))) root <- getwd()

# --- fail loud if graphviz is missing (don't silently produce nothing) -------------
if (nzchar(Sys.which("dot")) == FALSE) {
  stop("graphviz `dot` not found on PATH. Install it first: brew install graphviz")
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

# --- node styling: colour = axis, border = export status ----------------------------
exported <- getNamespaceExports("msamodel")
axis <- ifelse(grepl("_i(_|$)", nm), "site",
        ifelse(grepl("_n(_|$)", nm), "mode", "agnostic"))
axis[nm == "preprocess_spm_mode"] <- "mode"   # named mode, no _n token

pen  <- c(site = "#216b62", mode = "#9a6f14", agnostic = "#3a3f4b")[axis]
fill <- c(site = "#e6f1ef", mode = "#f6ecd6", agnostic = "#eceae6")[axis]

# --- assemble the DOT source --------------------------------------------------------
L <- c(
  "digraph msamodel {",
  "  labelloc=\"t\";",
  "  label=\"msamodel — function call graph  (arrow: A calls B; row = call depth; ",
  "          colour = axis: teal site / amber mode / grey agnostic; dashed box = @noRd internal)\";",
  "  fontname=\"Helvetica\"; fontsize=13;",
  "  rankdir=TB;",
  "  graph [splines=true, nodesep=0.35, ranksep=0.75];",
  "  node  [fontname=\"Courier\", fontsize=11, margin=\"0.12,0.06\"];",
  "  edge  [color=\"#9aa0aa\", arrowsize=0.7];"
)
for (i in seq_along(nm)) {
  is_exp <- nm[i] %in% exported
  L <- c(L, sprintf(
    "  \"%s\" [shape=box, style=\"filled,%s\", color=\"%s\", fillcolor=\"%s\", fontcolor=\"%s\"];",
    nm[i], if (is_exp) "solid" else "dashed",
    pen[i], if (is_exp) fill[i] else "#f2f1ee", pen[i]))
}
for (a in nm) for (b in colnames(m)[m[a, ] > 0]) {
  L <- c(L, sprintf("  \"%s\" -> \"%s\";", a, b))
}
L <- c(L, "}")

# --- render: DOT is a throwaway temp; the PNG + PDF are the deliverables ------------
outdir <- file.path(root, "dev", "preview")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dotfile <- tempfile(fileext = ".dot")
on.exit(unlink(dotfile), add = TRUE)
writeLines(L, dotfile)

png_out <- file.path(outdir, "callgraph.png")
pdf_out <- file.path(outdir, "callgraph.pdf")
stopifnot(system2("dot", c("-Tpng", "-Gdpi=160", shQuote(dotfile), "-o", shQuote(png_out))) == 0)
stopifnot(system2("dot", c("-Tpdf",              shQuote(dotfile), "-o", shQuote(pdf_out))) == 0)

cat(sprintf("call graph: %d functions, %d edges\n  %s\n  %s\n",
            length(nm), sum(m > 0), png_out, pdf_out))
