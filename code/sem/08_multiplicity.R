# ====================================================================
# Multiplicity. Run last, after H1, H2 and H3.
# Applies the pre-specified corrections and records the checksum of every input
# file, so that the corrected values can be traced to the exact result files they
# came from. H1 is a family of two under Bonferroni; H2 and H3 are separate
# families of twelve under Benjamini-Hochberg, because they ask different
# questions and a finding about cognition should not change status according to
# how many lifestyle variables happened to be measured.
# ====================================================================

# =====================================================================
# Multiplicity correction across the hypothesis families
#
# Reads the estimate files written by the hypothesis scripts, applies the
# correction each family calls for, and writes one table that states every
# confirmatory decision in the project. Nothing is estimated here.
#
# HOW TO RUN
#   Rscript code/sem/08_multiplicity.R      (or source())
#   Honours FUNCTIONAL_TAG, so it corrects the run you actually used.
#
# WHY THE TWO FAMILIES ARE CORRECTED DIFFERENTLY
#   H1 asks one question of two modalities: does the gap relate to
#   cognition? Two tests, and a false positive on either would be claimed
#   as support, so the correction is Bonferroni at alpha = .025 per test.
#   With two tests Bonferroni and Holm coincide, and the stricter control
#   is affordable because the family is small.
#
#   H2 asks six questions of two modalities: twelve tests, screening a set
#   of exposures rather than confirming a single claim. Bonferroni at
#   .05/12 = .004 would be the wrong instrument - it controls the chance of
#   any false positive at the cost of missing most true ones, which is not
#   the trade a screening analysis wants. Benjamini-Hochberg controls the
#   expected proportion of false discoveries among the rejections, which is
#   the quantity that matters when the output is a shortlist.
#
#   Both modalities stay in one H2 family. They are two measurements of one
#   hypothesis, not two hypotheses, and splitting them into two families of
#   six would loosen the correction without a substantive reason.
#
# WHAT COUNTS AS A TEST
#   H1  the g ~ BAG coefficient of the PRIMARY specification, per modality.
#       The supplementary specification is descriptive and is not corrected.
#   H2  the twelve rows flagged is_test = TRUE in the single-exposure
#       estimates, which are the six exposures x two modalities. The
#       mutually adjusted model in 05b is a sensitivity analysis and carries
#       no correction, so it is not read here.
#   H3  the same twelve, from 06_h3_single, with latent cognition as the
#       outcome. 06b is the sensitivity and is not read here either.
#
#   H2 and H3 are SEPARATE families. H2 asks what relates to the gap; H3
#   asks what relates to cognition given the gap. A result in one carries
#   no implication for the other, so correcting all twenty-four together
#   would be correcting across two questions instead of within one.
#
#   Covariates are never part of a family. Age, sex and education are in
#   every model as controls and no hypothesis is stated about them.
#
# ON READING THE OUTPUT
#   A q value is not a p value. q = .04 does not mean "4% chance this is a
#   false positive"; it means that if every result at or below this q were
#   declared a discovery, about 4% of those discoveries would be expected
#   to be false. It is a property of the list, not of the row.
#
#   A test that fails correction is reported with its uncorrected p and its
#   q, and described as not supported. It is not described as a trend.
#
# INPUT   results/sem/h1_estimates_<TAG>.csv
#         results/sem/h2_single_estimates_<TAG>.csv
#         results/sem/h3_estimates_<TAG>.csv
# OUTPUT  results/sem/multiplicity_<TAG>.csv
#         results/sem/multiplicity_log_<TAG>.txt
# =====================================================================

TAG  <- Sys.getenv("FUNCTIONAL_TAG", "pca")
ROOT <- "C:/Users/inaye/Desktop/MASTERTHESIS"
RES  <- file.path(ROOT, "results", "sem")

ALPHA_H1 <- 0.05   # split across the two modalities by Bonferroni -> .025 each
Q_H2     <- 0.05   # Benjamini-Hochberg level for H2, twelve tests
Q_H3     <- 0.05   # Benjamini-Hochberg level for H3, twelve tests

if (interactive())
  message("NOTE: run with Rscript or source(), not line by line.")

while (sink.number() > 0) sink()
try(sink(type = "message"), silent = TRUE)

con <- file(file.path(RES, sprintf("multiplicity_log_%s.txt", TAG)), open = "wt")
sink(con, split = TRUE); sink(con, type = "message")

log_closed <- FALSE
close_log <- function() {
  if (log_closed) return(invisible(NULL))
  log_closed <<- TRUE
  try(sink(type = "message"), silent = TRUE)
  while (sink.number() > 0) sink()
  try(close(con), silent = TRUE)
  invisible(NULL)
}
on.exit(close_log(), add = TRUE)


# ---------------------------------------------------------------------
# 1. Provenance
# ---------------------------------------------------------------------
this_script <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(normalizePath(f[1], winslash = "/", mustWork = FALSE))
  fr <- sys.frames()
  for (i in rev(seq_along(fr)))
    if (!is.null(fr[[i]]$ofile))
      return(normalizePath(fr[[i]]$ofile, winslash = "/", mustWork = FALSE))
  normalizePath(file.path(ROOT, "code", "sem", "08_multiplicity.R"),
                winslash = "/", mustWork = FALSE)
})

cat("=== provenance ===\n")
cat(sprintf("  %-18s %s\n", "run_time",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
cat(sprintf("  %-18s %s\n", "tag", TAG))
cat(sprintf("  %-18s %s\n", "r_version", R.version.string))
cat(sprintf("  %-18s %s\n", "script_md5",
            if (file.exists(this_script))
              unname(tools::md5sum(this_script)) else NA_character_))

# The inputs are checksummed too. A corrected table that cannot be traced to
# the estimate files it corrected is not worth more than the uncorrected one.
read_est <- function(file, what) {
  f <- file.path(RES, file)
  if (!file.exists(f)) {
    cat(sprintf("\n  %s not found: %s\n", what, basename(f)))
    return(NULL)
  }
  cat(sprintf("\n  %s: %s\n    md5 %s\n", what, basename(f),
              unname(tools::md5sum(f))))
  read.csv(f, stringsAsFactors = FALSE)
}

cat("\n=== inputs ===")
h1 <- read_est(sprintf("h1_estimates_%s.csv", TAG), "H1")
h2 <- read_est(sprintf("h2_single_estimates_%s.csv", TAG), "H2 (single-exposure)")
h3 <- read_est(sprintf("h3_estimates_%s.csv", TAG), "H3")

rows <- list()


# ---------------------------------------------------------------------
# 2. H1: two tests, Bonferroni
# ---------------------------------------------------------------------
# The coefficient is taken from the standardized rows so that the reported
# effect is on one scale across modalities. The decision does not depend on
# that choice: p is identical in the two row types up to rounding, because
# standardising rescales the estimate and its standard error together.

if (!is.null(h1)) {
  k <- h1$label == "b_bag" & h1$spec == "primary" & h1$type == "standardized"
  f <- h1[k, , drop = FALSE]

  if (!nrow(f)) {
    cat("\nH1: no primary b_bag rows found; nothing corrected\n")
  } else {
    f$p_adj  <- p.adjust(f$pvalue, method = "bonferroni")
    f$crit   <- ALPHA_H1 / nrow(f)
    f$decision <- ifelse(f$p_adj < ALPHA_H1, "supported", "not supported")

    cat(sprintf("\n\n=== H1: %d tests, Bonferroni at alpha = %.3f (%.4f per test) ===\n",
                nrow(f), ALPHA_H1, ALPHA_H1 / nrow(f)))
    print(data.frame(modality = f$modality, beta = round(f$est, 4),
                     se = round(f$se, 4), p = signif(f$pvalue, 3),
                     p_bonferroni = signif(f$p_adj, 3),
                     decision = f$decision),
          row.names = FALSE)

    rows[["H1"]] <- data.frame(
      family = "H1", method = "bonferroni", level = ALPHA_H1,
      modality = f$modality, term = f$rhs, label = "BAG -> cognition",
      estimate = f$est, se = f$se, ci.lower = f$ci.lower, ci.upper = f$ci.upper,
      p = f$pvalue, p_adjusted = f$p_adj, decision = f$decision,
      stringsAsFactors = FALSE)
  }
}


# ---------------------------------------------------------------------
# 3. H2: twelve tests, Benjamini-Hochberg
# ---------------------------------------------------------------------
# The unstandardized rows are the ones to report: the gap is the outcome and
# stays in years, so the estimate reads as years of gap per SD of exposure.

if (!is.null(h2)) {
  k <- h2$is_test %in% c(TRUE, "TRUE", "True") &
       h2$type == "unstandardized_years"
  f <- h2[k, , drop = FALSE]

  if (nrow(f) != 12)
    cat(sprintf("\nH2 WARNING: expected 12 test rows, found %d. The family is\n",
                nrow(f)),
        "  whatever is in the file, and the correction below reflects that.\n",
        sep = "")

  if (!nrow(f)) {
    cat("\nH2: no test rows found; nothing corrected\n")
  } else {
    f$p_adj <- p.adjust(f$pvalue, method = "BH")
    f$decision <- ifelse(f$p_adj < Q_H2, "supported", "not supported")
    o <- order(f$pvalue)
    f <- f[o, , drop = FALSE]

    cat(sprintf("\n\n=== H2: %d tests, Benjamini-Hochberg at q = %.2f ===\n",
                nrow(f), Q_H2))
    print(data.frame(modality = f$modality, exposure = f$exposure,
                     years_per_SD = round(f$est, 3),
                     ci = sprintf("[%.2f, %.2f]", f$ci.lower, f$ci.upper),
                     p = signif(f$pvalue, 3), q = signif(f$p_adj, 3),
                     decision = f$decision),
          row.names = FALSE)

    n_sig <- sum(f$decision == "supported")
    cat(sprintf("\n  %d of %d tests survive correction.\n", n_sig, nrow(f)))
    if (n_sig)
      cat("  Expected number of false discoveries among them: about",
          sprintf("%.1f", n_sig * Q_H2), "\n")

    rows[["H2"]] <- data.frame(
      family = "H2", method = "BH", level = Q_H2,
      modality = f$modality, term = f$rhs, label = f$exposure,
      estimate = f$est, se = f$se, ci.lower = f$ci.lower, ci.upper = f$ci.upper,
      p = f$pvalue, p_adjusted = f$p_adj, decision = f$decision,
      stringsAsFactors = FALSE)
  }
}


# ---------------------------------------------------------------------
# 4. H3, when it exists
# ---------------------------------------------------------------------
# Twelve tests, as for H2, and corrected the same way and for the same
# reason. They are a separate family: H2 asks what relates to the gap, H3
# asks what relates to cognition given the gap. A test that fails in one
# says nothing about the other, so pooling the twenty-four into one family
# would correct across two questions rather than within one.
#
# The standardized rows are the ones to read here: the outcome is latent
# cognition, so there is no natural unit and the unstandardized estimate
# depends on the scaling of g.

if (!is.null(h3)) {
  k <- h3$is_test %in% c(TRUE, "TRUE", "True") & h3$type == "standardized"
  f <- h3[k, , drop = FALSE]

  if (nrow(f) != 12)
    cat(sprintf("\nH3 WARNING: expected 12 test rows, found %d. The family is\n",
                nrow(f)),
        "  whatever is in the file, and the correction below reflects that.\n",
        sep = "")

  if (!nrow(f)) {
    cat("\nH3: no test rows found; nothing corrected\n")
  } else {
    f$p_adj <- p.adjust(f$pvalue, method = "BH")
    f$decision <- ifelse(f$p_adj < Q_H3, "supported", "not supported")
    f <- f[order(f$pvalue), , drop = FALSE]

    cat(sprintf("\n\n=== H3: %d tests, Benjamini-Hochberg at q = %.2f ===\n",
                nrow(f), Q_H3))
    print(data.frame(modality = f$modality, exposure = f$exposure,
                     beta = round(f$est, 3),
                     ci = sprintf("[%.2f, %.2f]", f$ci.lower, f$ci.upper),
                     p = signif(f$pvalue, 3), q = signif(f$p_adj, 3),
                     decision = f$decision),
          row.names = FALSE)

    n_sig <- sum(f$decision == "supported")
    cat(sprintf("\n  %d of %d tests survive correction.\n", n_sig, nrow(f)))
    if (n_sig)
      cat("  Expected number of false discoveries among them: about",
          sprintf("%.1f", n_sig * Q_H3), "\n")

    rows[["H3"]] <- data.frame(
      family = "H3", method = "BH", level = Q_H3,
      modality = f$modality, term = f$rhs, label = f$exposure,
      estimate = f$est, se = f$se, ci.lower = f$ci.lower, ci.upper = f$ci.upper,
      p = f$pvalue, p_adjusted = f$p_adj, decision = f$decision,
      stringsAsFactors = FALSE)
  }
}


# ---------------------------------------------------------------------
# 5. Save
# ---------------------------------------------------------------------
# Matrix is loaded elsewhere in this project and makes rbind an S4 generic
# whose dispatch indexes list(...)[[1]], so an empty list is an error rather
# than NULL. Guard it.
if (!length(rows)) {
  cat("\nNothing was corrected. Run 04, 05_h2_single and 06_h3_single first.\n")
} else {
  out <- do.call(rbind, rows)
  p <- file.path(RES, sprintf("multiplicity_%s.csv", TAG))
  write.csv(out, p, row.names = FALSE)
  cat("\nwrote", p, "\n")

  cat("\n=== every confirmatory decision in one place ===\n")
  print(data.frame(family = out$family, modality = out$modality,
                   label = out$label, p = signif(out$p, 3),
                   adjusted = signif(out$p_adjusted, 3),
                   decision = out$decision),
        row.names = FALSE)
}

cat("\nThe adjusted column is a Bonferroni-adjusted p for H1 and a BH q for\n")
cat("H2 and for H3. They are not the same quantity and should be labelled\n")
cat("separately in the thesis tables.\n")

close_log()
