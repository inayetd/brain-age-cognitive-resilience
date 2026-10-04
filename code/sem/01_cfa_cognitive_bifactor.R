# ====================================================================
# Measurement model 1 of 2. Run before any hypothesis test.
# Establishes latent general cognitive performance from seven indicators. The
# bifactor-(S-1) form anchors the general factor to the two processing-speed tasks,
# so that the executive and memory factors describe what those domains share beyond
# speed rather than an uninterpretable residual.
# This general factor is the outcome in H1 and H3.
# ====================================================================

# =====================================================================
# CFA1 / CFA2 - Latent cognitive performance
#
# Bifactor-(S-1) model over seven aging-sensitive indicators. The general
# factor is the measure of latent cognitive performance used in every
# downstream model.
#
# DATA PIPELINE
#   This script reads data/clean/analysis_phenotype.csv and nothing else.
#   That file is produced by code/brain_age/03_build_phenotype.py, which
#   performs every merge, recode and exclusion in one place: the age floor,
#   the neurological exclusions, the reverse-scoring of loneliness, and the
#   sign of Trail Making B. Reading raw release files here would risk
#   estimating the measurement model on a different sample from the one the
#   hypothesis tests use, which is exactly the inconsistency this layout
#   prevents.
#
# ESTIMATION
#   MLR (robust ML) with FIML for missing indicator values, matching the
#   estimator used for every other model in the protocol.
# =====================================================================

library(lavaan)
library(dplyr)

ROOT  <- "C:/Users/inaye/Desktop/MASTERTHESIS"
PHENO <- file.path(ROOT, "data", "clean", "analysis_phenotype.csv")
OUT   <- file.path(ROOT, "data", "clean")

RES <- file.path(ROOT, "results", "sem")
dir.create(RES, recursive = TRUE, showWarnings = FALSE)

# HOW TO RUN
#   Rscript code/sem/01_cfa_cognitive_bifactor.R   (or source())
#
# Not line by line from the console: the admissibility stop() in section 3
# is a guard, and submitting the file as a selection lets execution continue
# past it.
if (interactive())
  message("NOTE: run with Rscript or source(), not line by line.")

# Everything printed below is also written to a log file. The numbers that go
# into Table 1 are easy to lose in a scroll buffer, and a run that is not
# recorded has to be repeated to be checked.
while (sink.number() > 0) sink()
try(sink(type = "message"), silent = TRUE)

con <- file(file.path(RES, "cfa_cognitive_log.txt"), open = "wt")
sink(con, split = TRUE)      # split = TRUE keeps the console output as well
sink(con, type = "message")

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
# The input hash matters here specifically: this model was first run a month
# before the phenotype file was last rebuilt, so recording which version of
# it produced Table 1 is the only way to show that the measurement model and
# the hypothesis tests describe the same participants.
this_script <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(normalizePath(f[1], winslash = "/", mustWork = FALSE))
  fr <- sys.frames()
  for (i in rev(seq_along(fr)))
    if (!is.null(fr[[i]]$ofile))
      return(normalizePath(fr[[i]]$ofile, winslash = "/", mustWork = FALSE))
  normalizePath(file.path(ROOT, "code", "sem", "01_cfa_cognitive_bifactor.R"),
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
           file.path(RES, "cfa_cognitive_provenance.txt"))
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
# All seven are already oriented so that higher = better performance.
# trails_b_rev is the negative of log completion time, computed upstream;
# it is not recomputed here, so the sign convention lives in exactly one
# place in the whole project.

ind <- c("flanker",          # Flanker inhibitory control      (Toolbox SS)
         "card_sort",        # Dimensional change card sort    (Toolbox SS)
         "list_sorting",     # List sorting working memory     (Toolbox SS)
         "pattern_comp",     # Pattern comparison speed        (Toolbox SS)
         "picture_seq",      # Picture sequence memory         (Toolbox SS)
         "ravlt",            # RAVLT immediate recall          (raw sum)
         "trails_b_rev")     # Trail Making B, -log(seconds)

stopifnot(all(ind %in% names(d)))

cat("\nnon-missing per indicator (all", nrow(d), "rows):\n")
print(colSums(!is.na(d[ind])))
cat("at least one indicator:", sum(rowSums(!is.na(d[ind])) > 0), "\n")
cat("at least five (sem_eligible):", sum(d$sem_eligible), "\n")
cat("complete on all seven:", sum(complete.cases(d[ind])), "\n")

# Restrict to the analysis sample BEFORE fitting.
#
# sem_eligible marks the participants with at least five of the seven
# indicators, and it is the filter 04_h1_base.R applies before testing H1.
# Fitting the measurement model on a wider set would estimate the loadings
# on people the hypothesis test excludes, which is the inconsistency this
# script's header says the layout exists to prevent. The restriction is
# applied here rather than left to FIML because FIML would happily include
# a participant with one indicator, and a factor score resting on one
# indicator is not a measurement of the factor.
#
# FIML still does its work inside this sample, on the participants who
# have five or six indicators rather than all seven.
n_before <- nrow(d)
d <- d[d$sem_eligible, , drop = FALSE]
cat(sprintf("\nanalysis sample: %d of %d rows (sem_eligible, >= 5 indicators)\n",
            nrow(d), n_before))
if (!nrow(d)) stop("no participants are sem_eligible")

# Standardising the indicators does not change model fit, standardised
# loadings, or the factor correlations - only the unstandardised estimates.
# It is done because the protocol specifies it, and because it puts the
# seven indicators, which arrive on very different scales, on a common one.
# Computed within the analysis sample, so the scores are centred on the
# people the model actually describes.
d[ind] <- scale(d[ind])


# ---------------------------------------------------------------------
# 2. The bifactor-(S-1) model
# ---------------------------------------------------------------------
# What makes this S-1 rather than a symmetric bifactor: the two
# processing-speed indicators (pattern_comp, trails_b_rev) load on the
# general factor ONLY. They are the reference domain. The general factor
# is therefore interpreted as processing-speed-anchored general ability,
# and the two specific factors describe what executive and memory tasks
# share over and above it.
#
# a* constrains the two executive loadings to equality. A specific factor
# with only two indicators is not identified on its own; the constraint
# supplies the missing information.
#
# The three fixed zero covariances make all factors orthogonal, which is
# what allows the general factor to be read as a single dimension rather
# than as one correlated factor among several.

mod_cog <- '
  g    =~ flanker + card_sort + list_sorting + pattern_comp +
          picture_seq + ravlt + trails_b_rev
  Exec =~ a*flanker + a*card_sort
  Mem  =~ picture_seq + ravlt + list_sorting

  g    ~~ 0*Exec
  g    ~~ 0*Mem
  Exec ~~ 0*Mem
'

fit_cog <- cfa(mod_cog, data = d,
               estimator = "MLR",     # robust to the non-normality present here
               missing   = "fiml",    # uses everyone with >= 1 indicator
               std.lv    = TRUE)      # factor variances fixed at 1


# ---------------------------------------------------------------------
# 3. Admissibility - check before reading any estimate
# ---------------------------------------------------------------------
# The protocol specifies a fallback if this model is inadmissible: drop the
# memory and executive specific factors, in that order. These checks are
# what trigger that rule, so they run before anything is reported.

converged <- lavInspect(fit_cog, "converged")
theta     <- diag(lavInspect(fit_cog, "theta"))

# "psi" is not an accepted what= argument in lavaan 0.7; "cov.lv" returns the
# same latent covariance matrix. Same substitution as in 02.
#
# With std.lv = TRUE the factor variances are fixed at 1, so that diagonal
# carries no information and cannot fail. The check that can fail in a
# bifactor is on the other two quantities: a negative residual variance, and
# a standardised loading above one. Both are Heywood cases and both are what
# the fallback rule exists for, so both are tested here rather than the
# vacuous one.
psi     <- diag(lavInspect(fit_cog, "cov.lv"))
std_adm <- standardizedSolution(fit_cog)
load    <- std_adm$est.std[std_adm$op == "=~"]

cat("\n--- admissibility ---\n")
cat("converged                :", converged, "\n")
cat("min residual variance    :", round(min(theta), 4), "\n")
cat("factor variances         :", paste(round(psi, 3), collapse = ", "),
    "(fixed at 1 by std.lv)\n")
cat("max |standardised loading|:", round(max(abs(load)), 4), "\n")

if (!converged || min(theta) < 0 || max(abs(load)) > 1) {
  stop("Model is inadmissible. Apply the pre-specified fallback: remove the ",
       "memory specific factor, then the executive one, and report the ",
       "reduced structure as a departure.")
}


# ---------------------------------------------------------------------
# 4. Fit and estimates
# ---------------------------------------------------------------------
summary(fit_cog, fit.measures = TRUE, standardized = TRUE)

cat("\n--- fit indices for Table 1 ---\n")
print(round(fitMeasures(fit_cog, c(
  "chisq.scaled", "df", "pvalue.scaled",
  "cfi.robust", "tli.robust",
  "rmsea.robust", "rmsea.ci.lower.robust", "rmsea.ci.upper.robust",
  "srmr")), 3))

std <- standardizedSolution(fit_cog)
cat("\n--- standardised loadings ---\n")
print(std %>% filter(op == "=~") %>%
        select(factor = lhs, indicator = rhs, loading = est.std, se, pvalue))


# ---------------------------------------------------------------------
# 5. Strength of the general factor
# ---------------------------------------------------------------------
# ECV and omega hierarchical are computed explicitly rather than through a
# package function, so the formulas are visible and version-independent.
#
# ECV: share of common variance carried by the general factor.
# omegaH: share of total variance in a unit-weighted composite attributable
#         to the general factor, treating the specific factors as error.
#
# Both were developed for SYMMETRIC bifactor models. In an S-1 specification
# they describe the dominance of the reference-anchored general factor and do
# NOT establish domain-neutral unidimensionality. The protocol interprets
# them accordingly.

# Pull the standardised loadings as named vectors, indexed by indicator, so
# that every indicator is matched to its own general and specific loading
# rather than relying on the order in which lavaan happens to return them.
load_of <- function(f) {
  rows <- std$op == "=~" & std$lhs == f
  setNames(std$est.std[rows], std$rhs[rows])
}
lg <- load_of("g")            # one per indicator
le <- load_of("Exec")         # two indicators
lm <- load_of("Mem")          # three indicators

# Specific loading per indicator; zero for the two reference-domain tests,
# which by design load on g only.
spec        <- setNames(rep(0, length(lg)), names(lg))
spec[names(le)] <- le
spec[names(lm)] <- lm

# With standardised indicators, each indicator's total variance is 1, so its
# residual is whatever the two factors do not explain.
resid <- 1 - lg^2 - spec^2

# ECV: general-factor share of the variance the factors explain in common.
ecv <- sum(lg^2) / (sum(lg^2) + sum(spec^2))

# omega hierarchical: general-factor share of the variance of a unit-weighted
# composite of all seven indicators, with the specific factors counted as error.
omega_h <- sum(lg)^2 /
  (sum(lg)^2 + sum(le)^2 + sum(lm)^2 + sum(resid))

cat("\n--- general factor strength for Table 1 ---\n")
cat("ECV                      :", round(ecv, 3), "\n")
cat("omega hierarchical       :", round(omega_h, 3), "\n")
cat("executive loadings       :", paste(round(le, 3), collapse = ", "), "\n")
cat("memory loadings          :", paste(round(lm, 3), collapse = ", "), "\n")


# ---------------------------------------------------------------------
# 6. Factor scores - descriptive use only
# ---------------------------------------------------------------------
# In every hypothesis test the factor stays latent inside a single-step SEM.
# These scores exist only for descriptive figures and tables.
#
# Alignment matters: lavPredict returns one row per case that lavaan actually
# used, which is not necessarily every row of d. Asking lavaan which rows
# those were is the only safe way to merge them back. Guessing the rows would
# misassign scores silently.

fs  <- lavPredict(fit_cog, method = "regression")
idx <- lavInspect(fit_cog, "case.idx")
stopifnot(length(idx) == nrow(fs))

d$cog_g_score <- NA_real_
d$cog_g_score[idx] <- fs[, "g"]

cat("\nfactor scores assigned to", sum(!is.na(d$cog_g_score)), "participants\n")

write.csv(d[, c("participant_id", "cog_g_score")],
          file.path(OUT, "cognitive_factor_scores_V1.csv"), row.names = FALSE)


# ---------------------------------------------------------------------
# 7. Pre-specified fallback models, fitted for reference
# ---------------------------------------------------------------------
# Not used unless the primary model is inadmissible. Fitted here so their
# behaviour is known in advance rather than discovered under pressure.

mod_no_mem <- '
  g    =~ flanker + card_sort + list_sorting + pattern_comp +
          picture_seq + ravlt + trails_b_rev
  Exec =~ a*flanker + a*card_sort
  g    ~~ 0*Exec
'
mod_g_only <- '
  g =~ flanker + card_sort + list_sorting + pattern_comp +
       picture_seq + ravlt + trails_b_rev
'

for (nm in c("mod_no_mem", "mod_g_only")) {
  f <- cfa(get(nm), data = d, estimator = "MLR", missing = "fiml", std.lv = TRUE)
  cat("\n---", nm, "---\n")
  print(round(fitMeasures(f, c("cfi.robust", "tli.robust",
                               "rmsea.robust", "srmr")), 3))
}


# ---------------------------------------------------------------------
# 8. A note for the write-up
# ---------------------------------------------------------------------
# The two executive loadings are constrained equal by a* in the model
# syntax above, because a specific factor with only two indicators is not
# identified on its own. With the indicators standardised, this makes the
# two STANDARDISED executive loadings identical as well, which is why they
# print as the same number rather than as two slightly different ones.
#
# The analysis plan describes the executive factor as indicated by Flanker
# and Card Sort but does not state the equality constraint, and the two
# executive loadings it reports (.28 and .26) differ from each other, so
# they cannot have come from this specification. The constraint therefore
# has to be stated explicitly in the thesis Methods as a clarification of
# the registered model, not left for a reader to infer from two identical
# numbers.

cat("\n--- for the Methods section ---\n")
cat("The two executive loadings are constrained equal (a*) because a\n")
cat("two-indicator specific factor is not identified without it. State this\n")
cat("in the thesis; the analysis plan does not, and its reported .28/.26\n")
cat("came from a run that did not impose it.\n")

close_log()

