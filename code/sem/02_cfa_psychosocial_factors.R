# ====================================================================
# Measurement model 2 of 2.
# Establishes two correlated latent exposures, psychological well-being and social
# connectedness. Perceived stress has a single indicator and therefore enters the
# structural models as an observed variable instead.
# These factors are two of the six exposures in H2 and H3.
# ====================================================================

# =====================================================================
# CFA3 - Psychosocial measurement model
#
#   PWB  Psychological Well-Being  (positive affect, life satisfaction,
#                                   meaning and purpose)
#   SOC  Social Connectedness      (emotional support, instrumental
#                                   support, friendship, loneliness rev.)
#
# Perceived stress is a single indicator and cannot support a latent
# factor. It enters the structural models as an observed variable and is
# therefore not part of this measurement model.
#
# DATA PIPELINE
#   This script reads data/clean/analysis_phenotype.csv and nothing else,
#   for the same reason as the cognitive model: every merge, recode and
#   exclusion lives in code/brain_age/03_build_phenotype.py, so the
#   measurement model and the hypothesis tests are guaranteed to describe
#   the same participants. Loneliness is already reverse-scored there, as
#   loneliness_rev, and is not reversed again here.
#
# ESTIMATION
#   MLR with FIML, matching every other model in the protocol.
# =====================================================================

library(lavaan)
library(dplyr)

ROOT  <- "C:/Users/inaye/Desktop/MASTERTHESIS"
PHENO <- file.path(ROOT, "data", "clean", "analysis_phenotype.csv")
OUT   <- file.path(ROOT, "data", "clean")
RES   <- file.path(ROOT, "results", "sem")
dir.create(RES, recursive = TRUE, showWarnings = FALSE)

# HOW TO RUN
#   Rscript code/sem/02_cfa_psychosocial_factors.R   (or source())
#
# Not line by line from the console: the admissibility stop() in section 3
# is a guard, and submitting the file as a selection lets execution continue
# past it.
if (interactive())
  message("NOTE: run with Rscript or source(), not line by line.")

# Everything printed below is also written to a log file. Table 2 comes out
# of this script, and a run that is not recorded has to be repeated to be
# checked. The earlier version of this script wrote no log at all, so the
# psychosocial numbers reported in the analysis plan had no saved output
# behind them.
while (sink.number() > 0) sink()
try(sink(type = "message"), silent = TRUE)

con <- file(file.path(RES, "cfa_psychosocial_log.txt"), open = "wt")
sink(con, split = TRUE); sink(con, type = "message")

# on.exit only fires inside a function frame, so it covers Rscript and
# source() but not line-by-line submission. close_log() is idempotent and is
# also called explicitly at the end of the script.
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
# 0. Provenance
# ---------------------------------------------------------------------
# The input hash matters here specifically: the measurement models were
# first run a month before the phenotype file was last rebuilt, so
# recording which version of it produced Table 2 is the only way to show
# that the measurement model and the hypothesis tests describe the same
# participants.
this_script <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(normalizePath(f[1], winslash = "/", mustWork = FALSE))
  fr <- sys.frames()
  for (i in rev(seq_along(fr)))
    if (!is.null(fr[[i]]$ofile))
      return(normalizePath(fr[[i]]$ofile, winslash = "/", mustWork = FALSE))
  normalizePath(file.path(ROOT, "code", "sem", "02_cfa_psychosocial_factors.R"),
                winslash = "/", mustWork = FALSE)
})

PROV <- list(
  run_time       = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  input_file     = PHENO,
  input_md5      = if (file.exists(PHENO))
                     unname(tools::md5sum(PHENO)) else NA_character_,
  r_version      = R.version.string,
  lavaan_version = as.character(utils::packageVersion("lavaan")),
  script_path    = this_script,
  script_md5     = if (file.exists(this_script))
                     unname(tools::md5sum(this_script)) else NA_character_)

cat("=== provenance ===\n")
for (k in names(PROV)) cat(sprintf("  %-15s %s\n", k, PROV[[k]]))
writeLines(vapply(names(PROV),
                  function(k) sprintf("%-15s %s", k, PROV[[k]]), character(1)),
           file.path(RES, "cfa_psychosocial_provenance.txt"))
cat("\n")

# na.strings must include "". An empty field in a character column comes back
# from read.csv as the empty string rather than NA, so it survives every
# is.na() test and is then treated as a real value. The phenotype file has
# empty fields in education, scanner, apoe and mr_qc_code. Numeric columns
# are unaffected; read.csv already returns an empty numeric field as NA.
d <- read.csv(PHENO, stringsAsFactors = FALSE, na.strings = c("", "NA"))
cat("phenotype rows read:", nrow(d), "\n")

# The phenotype file is written by pandas, which encodes booleans as the
# strings "True" and "False". read.csv keeps them as character, so sum(),
# subsetting and filtering all fail or behave unexpectedly. Convert any column
# whose non-missing values are exactly those two strings.
bool_cols <- names(d)[vapply(d, function(x)
  is.character(x) && length(x[!is.na(x)]) > 0 &&
    all(x[!is.na(x)] %in% c("True", "False")), logical(1))]
d[bool_cols] <- lapply(d[bool_cols], function(x) x == "True")
cat("converted to logical:", paste(bool_cols, collapse = ", "), "\n")


# ---------------------------------------------------------------------
# 1. Indicators
# ---------------------------------------------------------------------
# All seven point the same way: higher = greater well-being or greater
# connectedness. loneliness_rev is (observed maximum - loneliness), which
# preserves the original scale range and so keeps the descriptive table
# interpretable; simply negating the score would fit identically but would
# report as negative means.

ind <- c("positive_affect",       # PWB
         "life_satisfaction",     # PWB
         "meaning_purpose",       # PWB
         "emotional_support",     # SOC
         "instrumental_support",  # SOC
         "friendship",            # SOC
         "loneliness_rev")        # SOC, already reversed upstream

stopifnot(all(ind %in% names(d)))

cat("\nnon-missing per indicator:\n"); print(colSums(!is.na(d[ind])))
cat("at least one indicator:", sum(rowSums(!is.na(d[ind])) > 0), "\n")
cat("complete on all seven :", sum(complete.cases(d[ind])), "\n")


# ---------------------------------------------------------------------
# 2. Primary model: two correlated factors
# ---------------------------------------------------------------------
# std.lv = TRUE fixes both factor variances at 1, so every loading is freely
# estimated and coefficients are read per standard deviation of the factor.
# Indicator residuals are left uncorrelated: no modification indices are used
# to add data-driven correlated residuals or cross-loadings, as pre-specified.
#
# The two factors are allowed to correlate, because they are theoretically
# related. Whether they are empirically distinguishable is tested in step 5.

mod2 <- '
  PWB =~ positive_affect + life_satisfaction + meaning_purpose
  SOC =~ emotional_support + instrumental_support + friendship + loneliness_rev
  PWB ~~ SOC
'

fit2 <- cfa(mod2, data = d,
            estimator = "MLR",
            missing   = "fiml",
            std.lv    = TRUE)


# ---------------------------------------------------------------------
# 3. Admissibility - check before reading any estimate
# ---------------------------------------------------------------------
# The protocol specifies that an inadmissible solution is reported as a
# model-fitting failure, with the released summary scores used only as a
# sensitivity analysis. These checks are what would trigger that.

converged <- lavInspect(fit2, "converged")
theta     <- diag(lavInspect(fit2, "theta"))
r_factors <- lavInspect(fit2, "cov.lv")["PWB", "SOC"]  # "psi" not accepted by lavaan 0.7

# The other Heywood case, tested here for the same reason as in 01: a
# standardised loading above one is inadmissible and is not caught by the
# residual variance or the factor correlation.
std_adm <- standardizedSolution(fit2)
load    <- std_adm$est.std[std_adm$op == "=~"]

cat("\n--- admissibility ---\n")
cat("converged                :", converged, "\n")
cat("min residual variance    :", round(min(theta), 4), "\n")
cat("factor correlation       :", round(r_factors, 3), "\n")
cat("max |standardised loading|:", round(max(abs(load)), 4), "\n")

if (!converged || min(theta) < 0 || abs(r_factors) >= 1 || max(abs(load)) > 1) {
  stop("Model is inadmissible. Report as a model-fitting failure and use the ",
       "released summary scores as a sensitivity analysis only.")
}


# ---------------------------------------------------------------------
# 4. Fit and estimates
# ---------------------------------------------------------------------
summary(fit2, fit.measures = TRUE, standardized = TRUE)

cat("\n--- fit indices for Table 2 ---\n")
print(round(fitMeasures(fit2, c(
  "chisq.scaled", "df", "pvalue.scaled",
  "cfi.robust", "tli.robust",
  "rmsea.robust", "rmsea.ci.lower.robust", "rmsea.ci.upper.robust",
  "srmr")), 3))

std <- standardizedSolution(fit2)

cat("\n--- standardised loadings ---\n")
print(std %>% filter(op == "=~") %>%
        select(factor = lhs, indicator = rhs, loading = est.std, se, pvalue))

cat("\n--- factor correlation ---\n")
print(std %>% filter(op == "~~", lhs != rhs) %>%
        select(lhs, rhs, r = est.std, ci.lower, ci.upper))

# Model-based composite reliability (omega) per factor. Computed explicitly so
# the formula is visible: squared sum of loadings over that plus the residuals.
omega <- function(f) {
  lam <- std$est.std[std$op == "=~" & std$lhs == f]
  sum(lam)^2 / (sum(lam)^2 + sum(1 - lam^2))
}
cat("\n--- reliability for Table 2 ---\n")
cat("omega, PWB               :", round(omega("PWB"), 3), "\n")
cat("omega, SOC               :", round(omega("SOC"), 3), "\n")


# ---------------------------------------------------------------------
# 5. Discriminant validity: do the two factors collapse into one?
# ---------------------------------------------------------------------
# A single factor is nested within the two-factor model (it is the special
# case where their correlation is 1), so the two can be compared by a
# likelihood-ratio test. With MLR, lavTestLRT applies the scaled correction
# automatically; the naive difference of scaled chi-squares would be wrong.

mod1 <- '
  G =~ positive_affect + life_satisfaction + meaning_purpose +
       emotional_support + instrumental_support + friendship + loneliness_rev
'
fit1 <- cfa(mod1, data = d, estimator = "MLR", missing = "fiml", std.lv = TRUE)

cat("\n--- one-factor vs two-factor (scaled LRT) ---\n")
print(lavTestLRT(fit1, fit2))


# ---------------------------------------------------------------------
# 6. Reference comparisons
# ---------------------------------------------------------------------
# (a) Social factor without loneliness. Loneliness is the only reversed
#     indicator on that factor, so it is worth knowing whether the factor
#     depends on it. Reported for completeness, not used.
mod2b <- '
  PWB =~ positive_affect + life_satisfaction + meaning_purpose
  SOC =~ emotional_support + instrumental_support + friendship
  PWB ~~ SOC
'
fit2b <- cfa(mod2b, data = d, estimator = "MLR", missing = "fiml", std.lv = TRUE)
cat("\n--- (a) social factor without loneliness ---\n")
print(round(fitMeasures(fit2b, c("cfi.robust", "tli.robust",
                                 "rmsea.robust", "srmr")), 3))

# (b) Listwise version. Omitting the missing argument makes lavaan delete
#     incomplete cases, so this shows whether FIML and complete-case analysis
#     agree. If they diverge, the missing-at-random assumption deserves a
#     closer look.
fit2_cc <- cfa(mod2, data = d, estimator = "MLR", std.lv = TRUE)
cat("\n--- (b) listwise deletion, n =", lavInspect(fit2_cc, "nobs"), "---\n")
print(round(fitMeasures(fit2_cc, c("cfi.robust", "tli.robust",
                                   "rmsea.robust", "srmr")), 3))


# ---------------------------------------------------------------------
# 7. Factor scores - descriptive use only
# ---------------------------------------------------------------------
# In the hypothesis tests both factors stay latent inside a single-step SEM,
# so that their measurement error is modelled rather than passed silently
# into the regression. These scores are for descriptive tables and plots.
#
# lavPredict returns one row per case lavaan actually used, which is not
# every row of d. lavInspect(..., "case.idx") gives those row numbers.
# Matching them any other way risks assigning scores to the wrong people
# without any error being raised.

fs  <- lavPredict(fit2, method = "regression")
idx <- lavInspect(fit2, "case.idx")
stopifnot(length(idx) == nrow(fs))

d$PWB_score <- NA_real_
d$SOC_score <- NA_real_
d$PWB_score[idx] <- fs[, "PWB"]
d$SOC_score[idx] <- fs[, "SOC"]

cat("\nfactor scores assigned to", sum(!is.na(d$PWB_score)), "participants\n")

write.csv(d[, c("participant_id", "PWB_score", "SOC_score")],
          file.path(OUT, "psychosocial_factor_scores_V1.csv"), row.names = FALSE)

cat("wrote", file.path(OUT, "psychosocial_factor_scores_V1.csv"), "\n")
cat("wrote", file.path(RES, "cfa_psychosocial_log.txt"), "\n")

close_log()
