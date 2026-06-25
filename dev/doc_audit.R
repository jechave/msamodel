#!/usr/bin/env Rscript
# dev/doc_audit.R — deterministic documentation-quality checker (Half I).
#
# Reads the GENERATED man/*.Rd (never R/), classifies each page, and runs the
# mechanical rubric checks from the doc-quality plan. Emits a per-page scorecard
# and an overall PASS/FAIL gate. The qualitative criteria (purpose-vs-mechanism,
# dangling article, jargon nuance, example realism) are judged separately by a
# subagent (Half II); this script owns only what a regex can settle reliably.
#
# Usage:  Rscript dev/doc_audit.R
# Exit:   0 if no page has a blocking (hard_fail) defect, 1 otherwise.

suppressWarnings(suppressMessages({
  library(tools)
  # Load the package (source tree) so we can read real formals() for the
  # param-coverage check. Fall back to the installed package if load_all fails.
  ok <- tryCatch({ pkgload::load_all(".", quiet = TRUE); TRUE },
                 error = function(e) FALSE)
  if (!ok) requireNamespace("msamodel", quietly = TRUE)
}))

man_dir <- "man"
stopifnot(dir.exists(man_dir))

# ---- helpers ---------------------------------------------------------------

# Flatten an Rd tag's content to plain text.
rd_text <- function(x) {
  if (is.null(x)) return("")
  paste(rapply(x, as.character, how = "unlist"), collapse = "")
}

# Pull the first top-level node of a given Rd tag from a parsed Rd object.
get_tag <- function(rd, tag) {
  tags <- vapply(rd, function(node) attr(node, "Rd_tag") %||% "", character(1))
  hit <- which(tags == tag)
  if (!length(hit)) return(NULL)
  rd[[hit[1]]]
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# All nodes of a tag (e.g. multiple \alias).
get_tags <- function(rd, tag) {
  tags <- vapply(rd, function(node) attr(node, "Rd_tag") %||% "", character(1))
  rd[hit <- which(tags == tag)]
}

# Banned developer-internal phrases that should not appear in user-facing prose.
BANNED <- c(
  "single source of truth", "MCMC path", "carries no method token",
  "no method token", "axis-agnostic", "forward-compatibility",
  "forward compatibility", "future method switch", "future `method`",
  "Predict-only", "predict-only", "building block reused", "recipe",
  "non-breaking addition"
)

# Mechanism words that, in a TITLE, signal implementation over purpose. Soft —
# flagged for the human/judge, not auto-blocking.
TITLE_MECH <- c("using matrix operations", "matrix operations", "one by one")

# ---- per-page audit --------------------------------------------------------

rd_files <- list.files(man_dir, pattern = "[.]Rd$", full.names = TRUE)

exports <- tryCatch(getNamespaceExports("msamodel"), error = function(e) character(0))

audit_one <- function(f) {
  rd <- parse_Rd(f)

  name    <- trimws(rd_text(get_tag(rd, "\\name")))
  title   <- trimws(rd_text(get_tag(rd, "\\title")))
  desc    <- trimws(rd_text(get_tag(rd, "\\description")))
  details <- trimws(rd_text(get_tag(rd, "\\details")))
  value   <- get_tag(rd, "\\value")
  examples<- get_tag(rd, "\\examples")
  seealso <- get_tag(rd, "\\seealso")
  format  <- get_tag(rd, "\\format")
  usage   <- rd_text(get_tag(rd, "\\usage"))
  docType <- trimws(rd_text(get_tag(rd, "\\docType")))
  keywords<- vapply(get_tags(rd, "\\keyword"), rd_text, character(1))
  concepts<- vapply(get_tags(rd, "\\concept"), rd_text, character(1))

  # Page type.
  is_data     <- docType == "data"
  is_internal <- any(trimws(keywords) == "internal")
  type <- if (is_data) "dataset" else if (is_internal) "internal" else "function"

  # Whole-page rendered text (what the user sees) for render-bug detection.
  page_txt <- paste(title, desc, details, rd_text(value), rd_text(examples),
                    rd_text(seealso), rd_text(format), collapse = "\n")

  # ---- mechanical booleans ----
  literal_bracket <- grepl("\\[[A-Za-z_][A-Za-z0-9_.]*\\(\\)\\]", page_txt) ||
                     grepl("\\[[^]]+\\]\\([^)]+\\)", page_txt)
  pipe_table      <- grepl("(^|\\n)\\s*\\|?[-:| ]{3,}\\|", page_txt) ||
                     grepl("\\|[^|\\n]+\\|[^|\\n]+\\|", page_txt)
  title_period    <- grepl("[.]$", title)
  title_repeats   <- nzchar(name) && grepl(name, title, fixed = TRUE)
  title_len       <- lengths(strsplit(title, "\\s+"))
  title_too_long  <- title_len > 12
  desc_eq_title   <- nzchar(desc) && identical(tolower(gsub("\\s+", " ", title)),
                                               tolower(gsub("\\s+", " ", desc)))
  banned_hit      <- BANNED[vapply(BANNED, function(p)
                       grepl(p, page_txt, ignore.case = TRUE), logical(1))]
  title_mech_hit  <- TITLE_MECH[vapply(TITLE_MECH, function(p)
                       grepl(p, title, ignore.case = TRUE), logical(1))]

  # Params: compare the function's real formals against the documented
  # \arguments. Names come from the live function (formals()), never from a
  # regex over \usage (defaults like c(0, 10) break paren-matching).
  arg_node <- get_tag(rd, "\\arguments")
  documented <- character(0)
  if (!is.null(arg_node)) {
    child_tags <- vapply(arg_node, function(n) attr(n, "Rd_tag") %||% "", character(1))
    for (k in which(child_tags == "\\item")) {
      it <- arg_node[[k]]                     # \item = list(name-node, desc-node)
      nm <- rd_text(it[[1]])                  # the name block, e.g. "x, y"
      nm <- unlist(strsplit(nm, "[,[:space:]]+"))
      documented <- c(documented, trimws(nm))
    }
  }
  documented <- documented[nzchar(documented)]

  formals_found <- character(0)
  if (type == "function" && nzchar(name) && exists(name, envir = asNamespace("msamodel"))) {
    fn <- get(name, envir = asNamespace("msamodel"))
    if (is.function(fn)) {
      formals_found <- names(formals(fn))
      formals_found <- formals_found[nzchar(formals_found) & formals_found != "..."]
    }
  }
  undocumented_param <- setdiff(formals_found, documented)

  missing_return   <- type == "function" && is.null(value)
  missing_examples <- type == "function" && is.null(examples)
  missing_format   <- type == "dataset"  && is.null(format)
  missing_seealso  <- type == "function" && is.null(seealso) && length(concepts) == 0
  # return shape: a \value that says only "tibble"/"list" with no columns/items.
  value_txt <- rd_text(value)
  return_no_columns <- type == "function" && !is.null(value) &&
    !grepl("\\\\item|`[A-Za-z]|column|:", value_txt) &&
    grepl("tibble|list|data\\.?frame", value_txt, ignore.case = TRUE)

  # ---- blocking (hard) fails ----
  hard <- c()
  if (literal_bracket)             hard <- c(hard, "literal_bracket")
  if (pipe_table)                  hard <- c(hard, "pipe_table")
  if (desc_eq_title)               hard <- c(hard, "desc==title")
  if (length(undocumented_param))  hard <- c(hard, paste0("undocumented_param:",
                                              paste(undocumented_param, collapse="/")))
  if (missing_return)              hard <- c(hard, "missing_return")
  if (missing_examples)            hard <- c(hard, "missing_examples")
  if (missing_format)              hard <- c(hard, "missing_format")

  # ---- soft/minor + judgment-deferred flags ----
  soft <- c()
  if (title_period)        soft <- c(soft, "title_period")
  if (title_repeats)       soft <- c(soft, "title_repeats_name")
  if (title_too_long)      soft <- c(soft, paste0("title_len:", title_len))
  if (length(banned_hit))  soft <- c(soft, paste0("jargon:", paste(banned_hit, collapse="/")))
  if (length(title_mech_hit)) soft <- c(soft, paste0("title_mechanism:",
                                          paste(title_mech_hit, collapse="/")))
  if (missing_seealso)     soft <- c(soft, "missing_seealso")
  if (return_no_columns)   soft <- c(soft, "return_no_columns")

  # Jargon is a category-3 FAIL, not merely cosmetic: treat as hard for scoring.
  if (length(banned_hit))  hard <- c(hard, paste0("jargon:", paste(banned_hit, collapse="/")))

  # ---- score ----
  # 5 = no hard, no soft. 4 = no hard, only cosmetic soft. cap 2 on any hard.
  score <- if (length(hard)) 2L
           else if (length(soft)) 4L
           else 5L

  list(
    page = name, type = type, score = score,
    hard = paste(hard, collapse = "; "),
    soft = paste(soft, collapse = "; ")
  )
}

results <- lapply(rd_files, function(f) {
  tryCatch(audit_one(f), error = function(e)
    list(page = basename(f), type = "ERROR", score = NA_integer_,
         hard = paste("parse error:", conditionMessage(e)), soft = ""))
})

# ---- report ----------------------------------------------------------------

ord <- order(vapply(results, function(r) r$score %||% -1L, integer(1)),
             vapply(results, function(r) r$page, character(1)))
results <- results[ord]

cat(sprintf("\n%-34s %-9s %-5s %s\n", "PAGE", "TYPE", "SCORE", "FAILS"))
cat(strrep("-", 100), "\n")
for (r in results) {
  fails <- paste(c(r$hard, r$soft)[nzchar(c(r$hard, r$soft))], collapse = " | ")
  cat(sprintf("%-34s %-9s %-5s %s\n", r$page, r$type,
              ifelse(is.na(r$score), "ERR", r$score), fails))
}

n_total <- length(results)
n_5     <- sum(vapply(results, function(r) isTRUE(r$score == 5L), logical(1)))
n_hard  <- sum(vapply(results, function(r) nzchar(r$hard), logical(1)))

cat(strrep("-", 100), "\n")
cat(sprintf("TOTAL %d pages | %d at 5/5 | %d with blocking fails\n",
            n_total, n_5, n_hard))

# Coverage cross-check: every export must have a page.
documented_names <- vapply(results, function(r) r$page, character(1))
missing_pages <- setdiff(exports, documented_names)
if (length(missing_pages)) {
  cat("WARNING: exported but no man page: ", paste(missing_pages, collapse = ", "), "\n")
}

quit(status = if (n_hard == 0) 0 else 1)
