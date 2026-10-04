# ====================================================================
# Missing data. Run after the measurement models, before the hypothesis tests.
# Household income is 25% missing and is the variable carrying the main finding, so
# how it is handled matters. Chained equations impute the observed predictors; the
# latent indicators are deliberately restored to missing afterwards and handled by
# full-information maximum likelihood inside lavaan instead.
# Produces the 50 completed datasets that every hypothesis test pools over.
# ====================================================================

# =====================================================================
# MI - multiple imputation of the observed external predictors
#
# Produces the completed datasets that every structural model reads.
# Nothing is modelled here.
#
# WHAT IS IMPUTED AND WHAT IS NOT
#   Two mechanisms are used, for different parts of the model.
#
#   Chained equations, kept:
#     log household income, sleep quality, perceived stress, education,
#     physical activity category. These are observed predictors on mixed
#     measurement scales, which is what chained equations handle well.
#
#   Chained equations, then discarded:
#     the seven cognitive indicators, the seven psychosocial indicators,
#     and the modality-specific brain age gap. These are imputed only so
#     that they can inform the imputation of the variables above; their
#     imputed values are then overwritten with NA, restoring the original
#     pattern. Downstream they are handled by full-information maximum
#     likelihood inside the structural equation model, which is the
#     natural treatment for indicators of a latent variable.
#
#   Never imputed:
#     age, sex, and the acquisition summaries, which are complete.
#
#   The consequence is that the imputation and analysis models are not
#   fully congenial. This is recorded in the analysis plan as a limitation.
#
# THE BRAIN AGE GAP IS AN AUXILIARY VARIABLE HERE, NOT AN OUTCOME
#   The gap is a predictor in every downstream model, but inside this
#   script it is treated exactly like the cognitive and psychosocial
#   indicators: imputed so that it can inform the imputation of income,
#   education, sleep, stress and activity, then restored to NA in all
#   M datasets before anything is saved.
#
#   It is never imputed for analysis. A participant with no usable
#   imaging has no brain age gap; there is no information anywhere in
#   the data from which one could be constructed, and a value invented
#   by predictive mean matching for the central predictor of every
#   hypothesis would be a fabrication carried straight into the
#   estimates. The participants who arrive here without a gap are the
#   six who fail the imaging criteria, and they must leave the
#   structural models on that ground, not be filled in.
#
#   The restoration is verified per dataset in section 5 and the script
#   stops if any gap cell survived.
#
# CONTRACT WITH THE DOWNSTREAM SEM SCRIPTS
#   Each saved object is a list with three fields:
#     $completed  list of M data frames. THIS is what the SEM scripts read.
#     $mids       the raw mice object, for convergence diagnostics ONLY.
#     $provenance seed, package and R versions, script checksum.
#
#   The SEM scripts must iterate over $completed and pool with Rubin's
#   rules, and must never call complete(x$mids, i). The mids object still
#   holds the imputed cognitive, psychosocial and brain age gap values,
#   because mice has no notion that they were meant to be discarded.
#   Reading from it would silently feed fabricated indicator values and a
#   fabricated brain age gap into the analysis. Section 6 prints the
#   difference between the two so that it is visible in the log.
#
#     functional models  ->  readRDS("mice_<TAG>_functional.rds")$completed
#     structural models  ->  readRDS("mice_<TAG>_structural.rds")$completed
#
# SEPARATE RUNS PER MODALITY
#   The functional run includes the functional gap as a predictor and
#   excludes the structural one, and the reverse for the structural run.
#   Carrying the other modality's gap into the imputation would let one
#   modality inform the other's analysis sample, which is the asymmetry
#   the H4 comparison is meant to avoid.
#
# INPUT   data/clean/analysis_table_<TAG>.csv   (from 06_merge_bags.py)
# OUTPUT  data/clean/mice_<TAG>_functional.rds
#         data/clean/mice_<TAG>_structural.rds
#         (convergence diagnostics are produced by 03b_convergence.R)
#         results/sem/mice_provenance_<TAG>.txt
#         results/sem/mice_log_<TAG>.txt
# =====================================================================

library(mice)

TAG  <- Sys.getenv("FUNCTIONAL_TAG", "pca")   # which functional run to use
ROOT <- "C:/Users/inaye/Desktop/MASTERTHESIS"
IN   <- file.path(ROOT, "data", "clean", sprintf("analysis_table_%s.csv", TAG))
OUTD <- file.path(ROOT, "data", "clean")
RES  <- file.path(ROOT, "results", "sem")
dir.create(RES, recursive = TRUE, showWarnings = FALSE)

M_IMP  <- 50    # imputed datasets, as pre-specified
N_ITER <- 50    # was 20
SEED   <- 2026

PSRF_MAX <- 1.05   # convergence threshold for the potential scale reduction factor

# The analysis plan specifies 20 iterations per chain. A run with any other
# value is a convergence sensitivity analysis, not the pre-specified one, so
# it writes under its own names and cannot overwrite the registered output.
RUN <- if (N_ITER == 20L) TAG else sprintf("%s_it%d", TAG, N_ITER)
if (RUN != TAG)
  message("N_ITER = ", N_ITER, " departs from the pre-specified 20; ",
          "output will be written under the label '", RUN, "'.")

# HOW TO RUN THIS SCRIPT
#   Rscript code/sem/03_impute.R      (or source() from the console)
#
# Do not select the file and submit it to the console line by line. RStudio
# sends the whole selection as separate top-level expressions, so an error
# aborts only the expression that raised it and execution continues with the
# next line. Every stop() below - the education level check and the
# post-imputation sanity checks - is a guard that is supposed to halt the run,
# and submitting the file this way silently disables all of them.
if (interactive())
  message("NOTE: run with Rscript or source(). Submitting this file to the ",
          "console line by line lets execution continue past stop() and ",
          "disables every guard in it.")

# Any sink left open by an aborted earlier run would swallow this run's output.
while (sink.number() > 0) sink()
try(sink(type = "message"), silent = TRUE)

con <- file(file.path(RES, sprintf("mice_log_%s.txt", RUN)), open = "wt")
sink(con, split = TRUE); sink(con, type = "message")

# on.exit only fires inside a function frame, so it covers Rscript and source()
# but not line-by-line submission. close_log() is therefore idempotent and is
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
# Recorded before anything runs and saved inside every .rds, so that a
# result can always be traced back to the exact script and package
# versions that produced it.

this_script <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(normalizePath(f[1], winslash = "/", mustWork = FALSE))
  fr <- sys.frames()
  for (i in rev(seq_along(fr)))
    if (!is.null(fr[[i]]$ofile))
      return(normalizePath(fr[[i]]$ofile, winslash = "/", mustWork = FALSE))
  # Neither Rscript nor source(): fall back to the canonical location.
  normalizePath(file.path(ROOT, "code", "sem", "03_impute.R"),
                winslash = "/", mustWork = FALSE)
})

script_md5 <- if (file.exists(this_script))
  unname(tools::md5sum(this_script)) else NA_character_

git_commit <- tryCatch({
  g <- suppressWarnings(system2("git", c("-C", shQuote(ROOT), "rev-parse", "HEAD"),
                                stdout = TRUE, stderr = FALSE))
  if (length(g) && nzchar(g[1])) g[1] else "not a git repository"
}, error = function(e) "git unavailable")

git_dirty <- tryCatch({
  s <- suppressWarnings(system2("git", c("-C", shQuote(ROOT), "status", "--porcelain"),
                                stdout = TRUE, stderr = FALSE))
  if (length(s)) "yes (uncommitted changes present)" else "no"
}, error = function(e) NA_character_)

PROV <- list(
  run_time      = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  tag           = TAG,
  input_file    = IN,
  input_md5     = if (file.exists(IN)) unname(tools::md5sum(IN)) else NA_character_,
  seed          = SEED,
  m_imputations = M_IMP,
  n_iterations  = N_ITER,
  psrf_max      = PSRF_MAX,
  r_version     = R.version.string,
  platform      = R.version$platform,
  mice_version  = as.character(utils::packageVersion("mice")),
  script_path   = this_script,
  script_md5    = script_md5,
  git_commit    = git_commit,
  git_dirty     = git_dirty)

cat("=== provenance ===\n")
for (k in names(PROV)) cat(sprintf("  %-14s %s\n", k, PROV[[k]]))
cat("\n")

writeLines(
  vapply(names(PROV), function(k) sprintf("%-14s %s", k, PROV[[k]]), character(1)),
  file.path(RES, sprintf("mice_provenance_%s.txt", RUN)))


# na.strings must include "". An empty field in a character column comes back
# from read.csv as the empty string, not NA, so it survives every is.na() test
# and is then treated as a real category. In this table education, scanner,
# apoe and mr_qc_code all contain empty fields; without this the 26 blank
# education values reach the level check and report themselves as an
# unrecognised category. Numeric columns are unaffected, read.csv already
# returns an empty numeric field as NA.
d <- read.csv(IN, stringsAsFactors = FALSE, na.strings = c("", "NA"))
cat("rows read:", nrow(d), "\n")

# pandas writes booleans as the strings "True"/"False"; read.csv keeps them as
# character, so any logical column has to be converted before it can be used.
bool_cols <- names(d)[vapply(d, function(x)
  is.character(x) && length(x[!is.na(x)]) > 0 &&
    all(x[!is.na(x)] %in% c("True", "False")), logical(1))]
d[bool_cols] <- lapply(d[bool_cols], function(x) x == "True")
cat("converted to logical:", paste(bool_cols, collapse = ", "), "\n")


# ---------------------------------------------------------------------
# 1. Recode the three variables that do not arrive in modelling form
# ---------------------------------------------------------------------
# Education is released as text. The order below is the substantive order of
# the six levels and is fixed here rather than left to alphabetical sorting,
# which would scramble it: sorted alphabetically, "Graduate Degree" would be
# coded 2 and "Some Secondary Level Education" 6, a coding that correlates
# r = -0.67 with the correct one and raises no error. Mapped in the order
# below it reproduces the descriptives given in the protocol
# (M = 4.81, SD = 1.04, median = 5).
EDU_LEVELS <- c(
  "Primary Level Education",
  "Some Secondary Level Education",
  "Secondary Level Education or Equivalent",
  "Some College / Associates Degree / Trade or Vocational Training",
  "Four Years College or Undergraduate Degree",
  "Graduate Degree")

# factor() silently maps any label outside levels= to NA, so an unexpected
# label would turn into missing data instead of an error. Check first.
unknown <- setdiff(unique(d$education[!is.na(d$education)]), EDU_LEVELS)
if (length(unknown))
  stop("unrecognised education level(s): ", paste(unknown, collapse = " | "))

d$education_ord <- factor(d$education, levels = EDU_LEVELS, ordered = TRUE)
d$education_num <- as.integer(d$education_ord)

cat(sprintf("\neducation: M = %.2f, SD = %.2f, median = %d, n = %d\n",
            mean(d$education_num, na.rm = TRUE), sd(d$education_num, na.rm = TRUE),
            median(d$education_num, na.rm = TRUE), sum(!is.na(d$education_num))))
print(table(d$education_ord, useNA = "ifany"))

d$sex <- factor(d$sex)

# Physical activity is ordinal for the purpose of imputation but enters the
# structural models as indicator variables with the low category as reference,
# because the hypotheses specify a high-versus-low contrast rather than a
# linear trend. Both forms are carried forward.
d$ipaq_ord  <- factor(d$ipaq_category, levels = 1:3, ordered = TRUE)
d$ipaq_mod  <- as.integer(d$ipaq_category == 2)
d$ipaq_high <- as.integer(d$ipaq_category == 3)


# ---------------------------------------------------------------------
# 2. Variable roles
# ---------------------------------------------------------------------
COG <- c("flanker", "card_sort", "list_sorting", "pattern_comp",
         "picture_seq", "ravlt", "trails_b_rev")
PSY <- c("positive_affect", "life_satisfaction", "meaning_purpose",
         "emotional_support", "instrumental_support", "friendship",
         "loneliness_rev")
KEEP_IMPUTED <- c("income_log", "sleep_quality", "perceived_stress",
                  "education_ord", "ipaq_ord")
COMPLETE     <- c("age", "sex", "rsfmri_pct_complete", "rsfmri_runs")

missing_cols <- setdiff(c(COG, PSY, KEEP_IMPUTED, COMPLETE,
                          "participant_id", "pedigree", "sem_eligible"), names(d))
if (length(missing_cols))
  stop("columns absent from the analysis table: ", paste(missing_cols, collapse = ", "))

cat("\nmissingness before imputation:\n")
NA_BEFORE <- setNames(integer(0), character(0))
for (v in c(KEEP_IMPUTED, COG, PSY, "BAG_functional", "BAG_structural")) {
  k <- sum(is.na(d[[v]]))
  NA_BEFORE[v] <- k
  cat(sprintf("  %-22s %4d missing (%4.1f%%)\n", v, k, 100 * k / nrow(d)))
}


# ---------------------------------------------------------------------
# 3. One imputation per modality
# ---------------------------------------------------------------------
impute_for <- function(gap_var) {

  vars <- c(KEEP_IMPUTED, COG, PSY, COMPLETE, gap_var)
  sub  <- d[, vars]

  # The variables in COMPLETE are given the empty imputation method, which
  # means mice leaves them exactly as they are. That is correct only if they
  # really are complete. If one of them arrived with missing values, mice
  # would not fill it and would not complain, and every downstream model
  # would quietly drop those cases through listwise deletion inside lavaan,
  # with no trace in this log. Assert the premise instead of assuming it.
  if (anyNA(sub[COMPLETE])) {
    bad <- colSums(is.na(sub[COMPLETE]))
    stop("variables designated as complete contain missing values: ",
         paste(sprintf("%s (%d)", names(bad)[bad > 0], bad[bad > 0]),
               collapse = ", "),
         "\n  Either fix them upstream in 03_build_phenotype.py or move them ",
         "out of COMPLETE and give them an imputation method.")
  }

  # Which cells must be restored to missing afterwards. Recorded before mice
  # runs, because after it every cell is filled. The brain age gap is in this
  # list for the reason given in the header: it is an auxiliary variable here,
  # informing the imputation of the other predictors and nothing else.
  restore <- lapply(c(COG, PSY, gap_var), function(v) is.na(sub[[v]]))
  names(restore) <- c(COG, PSY, gap_var)

  meth <- make.method(sub)
  meth[c("income_log", "sleep_quality", "perceived_stress")] <- "pmm"
  meth[COG] <- "pmm"
  meth[PSY] <- "pmm"
  meth[gap_var] <- "pmm"          # auxiliary only, discarded in the loop below
  meth["education_ord"] <- "polr"     # ordinal
  meth["ipaq_ord"]      <- "polr"     # ordinal
  meth[COMPLETE] <- ""                # complete, never imputed

  pred <- make.predictorMatrix(sub)
  diag(pred) <- 0

  cat("\n--- imputing with", gap_var, "as an auxiliary predictor ---\n")
  cat("methods:\n"); print(meth[meth != ""])

  imp <- mice(sub, m = M_IMP, maxit = N_ITER, method = meth,
              predictorMatrix = pred, seed = SEED, printFlag = FALSE)

  if (length(imp$loggedEvents)) {
    cat("\nmice logged events (constant or collinear predictors):\n")
    print(head(imp$loggedEvents, 20))
  }

  # Discard the imputed values for the variables that are only there to inform
  # the others, and restore the original missing-data pattern in every dataset.
  completed <- lapply(seq_len(M_IMP), function(i) {
    x <- complete(imp, i)
    for (v in names(restore)) x[restore[[v]], v] <- NA
    # carry through the identifiers and the numeric/indicator forms the models need
    x$participant_id <- d$participant_id
    x$pedigree       <- d$pedigree
    x$sem_eligible   <- d$sem_eligible
    x$education_num  <- as.integer(x$education_ord)
    x$ipaq_mod       <- as.integer(x$ipaq_ord == 2)
    x$ipaq_high      <- as.integer(x$ipaq_ord == 3)
    x
  })

  list(mids = imp, completed = completed, restore = names(restore),
       gap_var = gap_var,
       modality = if (gap_var == "BAG_functional") "functional" else "structural",
       provenance = PROV)
}

fun  <- impute_for("BAG_functional")
str_ <- impute_for("BAG_structural")


# ---------------------------------------------------------------------
# 4. Post-imputation sanity check, both modalities, every dataset
# ---------------------------------------------------------------------
# Runs BEFORE the results are saved and before any diagnostic that depends
# on an external package, because it is the check that decides whether the
# imputations may be used at all, and because it has no dependencies of its
# own that could fail.
#
# Three things must hold in each of the M datasets:
#   the kept predictors are complete;
#   the cognitive, psychosocial and brain age gap columns carry exactly the
#     missingness they had before imputation, not one cell fewer;
#   the derived numeric and indicator forms agree with the factors they
#     were built from.
# A failure here means the restoration step did not do what it claims, so
# the script stops rather than saving datasets that look complete.

sanity <- function(obj) {
  cat("\n--- sanity check:", obj$modality, "---\n")
  gap <- obj$gap_var

  for (i in seq_len(M_IMP)) {
    x <- obj$completed[[i]]

    left <- vapply(KEEP_IMPUTED, function(v) sum(is.na(x[[v]])), integer(1))
    if (any(left > 0))
      stop(sprintf("dataset %d (%s): kept predictor(s) still missing: %s",
                   i, obj$modality,
                   paste(sprintf("%s=%d", names(left)[left > 0], left[left > 0]),
                         collapse = ", ")))

    for (v in obj$restore) {
      got <- sum(is.na(x[[v]]))
      if (got != NA_BEFORE[[v]])
        stop(sprintf("dataset %d (%s): %s has %d missing after restoration, expected %d",
                     i, obj$modality, v, got, NA_BEFORE[[v]]))
    }

    if (!identical(x$education_num, as.integer(x$education_ord)))
      stop(sprintf("dataset %d (%s): education_num does not match education_ord",
                   i, obj$modality))
    if (!identical(x$ipaq_mod, as.integer(x$ipaq_ord == 2)) ||
        !identical(x$ipaq_high, as.integer(x$ipaq_ord == 3)))
      stop(sprintf("dataset %d (%s): ipaq indicators do not match ipaq_ord",
                   i, obj$modality))
  }

  cat(sprintf("  all %d datasets pass\n", M_IMP))
  x <- obj$completed[[1]]
  cat("  dataset 1, missing counts:\n")
  for (v in c(KEEP_IMPUTED, "ravlt", "life_satisfaction", gap))
    cat(sprintf("    %-20s %4d  (before imputation: %d)\n",
                v, sum(is.na(x[[v]])), NA_BEFORE[[v]]))
  invisible(TRUE)
}

sanity(fun)
sanity(str_)


# ---------------------------------------------------------------------
# 5. Save
# ---------------------------------------------------------------------
# Deliberately before the convergence diagnostics. The imputation itself
# takes around twenty minutes; a diagnostic that fails afterwards should
# cost a plot, not the whole run. The sanity checks above have already
# established that these objects are fit to save.

saveRDS(fun,  file.path(OUTD, sprintf("mice_%s_functional.rds", RUN)))
saveRDS(str_, file.path(OUTD, sprintf("mice_%s_structural.rds", RUN)))

cat("\nwrote", file.path(OUTD, sprintf("mice_%s_functional.rds", RUN)), "\n")
cat("wrote", file.path(OUTD, sprintf("mice_%s_structural.rds", RUN)), "\n")
cat("wrote", file.path(RES, sprintf("mice_provenance_%s.txt", RUN)), "\n")


# ---------------------------------------------------------------------
# 6. What the mids object would have given instead
# ---------------------------------------------------------------------
# The SEM scripts must read $completed. This block shows, in the log, what
# reading $mids directly would have done: the brain age gap and the
# indicators come back filled, because mice has no notion that those
# imputations were auxiliary. Printed rather than asserted, so the
# difference is on record every time the script runs.

cat("\n=== why the SEM scripts must read $completed, not $mids ===\n")
for (obj in list(fun, str_)) {
  raw  <- complete(obj$mids, 1)
  done <- obj$completed[[1]]
  cat(sprintf("  %-11s %-18s mids: %4d missing   completed: %4d missing\n",
              obj$modality, obj$gap_var,
              sum(is.na(raw[[obj$gap_var]])), sum(is.na(done[[obj$gap_var]]))))
  cat(sprintf("  %-11s %-18s mids: %4d missing   completed: %4d missing\n",
              "", "ravlt",
              sum(is.na(raw$ravlt)), sum(is.na(done$ravlt))))
}
cat("  The mids columns are complete because mice imputed them. Those values\n")
cat("  are auxiliary and must not reach an estimate.\n")

# ---------------------------------------------------------------------
# 7. Convergence is diagnosed by a separate script
# ---------------------------------------------------------------------
# The trace plots and the mixing diagnostic live in 03b_convergence.R,
# which reads the two .rds files written above. Keeping them out of this
# script means a diagnostic that fails cannot cost the imputation, and
# that the plots can be regenerated in seconds without repeating a run
# that takes twenty minutes. Run it next:
#
#     Rscript code/sem/03b_convergence.R
#
# Do not treat these imputations as usable until it reports that the
# chains have mixed.

cat("\n=== next step ===\n")
cat("  Convergence diagnostics: Rscript code/sem/03b_convergence.R\n")
cat("  The imputations are saved but are not to be used until it passes.\n")

cat("\nSEM scripts read the $completed field of these objects. See the contract\n")
cat("in the header.\n")

close_log()
