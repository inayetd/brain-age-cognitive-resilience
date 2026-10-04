# ====================================================================
# H2, SENSITIVITY. The same six exposures entered simultaneously.
# Answers a different question from 05: each coefficient is now the association
# unique to that exposure. Reported to show that the primary finding does not
# depend on the one-at-a-time specification.
# ====================================================================

# =====================================================================
# H2, SENSITIVITY - Lifestyle and psychosocial correlates of the brain
#                   age gap, all six exposures entered together
#
# The six exposures are entered SIMULTANEOUSLY, so each coefficient is
# the association of that exposure with the gap after the other five and
# the covariates. Fitted separately for the functional and the structural
# gap.
#
# THIS IS NOT THE PRIMARY ANALYSIS
#   The primary analysis is 05_h2_single.R, which fits one exposure per
#   model. This script answers a different question and is reported
#   beside it, not instead of it:
#
#     05_h2_single   is X associated with the gap?
#                    covariate-adjusted total association
#
#     05b_h2_joint   is X associated with the gap over and above the
#                    other five exposures?
#                    the part of X not shared with them
#
#   The six exposures are substantially intercorrelated, so a coefficient
#   that shrinks between the two is the expected result and is itself the
#   finding: the association is largely shared across psychosocial
#   domains rather than specific to one of them. Shrinkage here is not
#   evidence that the primary association is absent, and a sign flip
#   under mutual adjustment is suppression, not a contradiction.
#
#   Nothing in this script decides support for H2. The Benjamini-Hochberg
#   family is the twelve rows from 05_h2_single. The estimates here are
#   descriptive and go in the supplementary table.
#
# HOW TO RUN
#   Rscript code/sem/05b_h2_joint.R      (or source())
#
# Do not submit this file to the console line by line. Every stop() below
# is a guard meant to halt the run, and line-by-line submission lets
# execution continue past it.
#
# THE SIX EXPOSURES
#   (i)   log household income          income_z    higher -> smaller gap
#   (ii)  sleep quality                 sleep_z     higher -> smaller gap
#   (iii) physical activity, high vs low ipaq_high  higher -> smaller gap
#   (iv)  psychological well-being      PWB latent  higher -> smaller gap
#   (v)   social connectedness          SOC latent  higher -> smaller gap
#   (vi)  perceived stress              stress_z    higher -> LARGER gap
#
#   ipaq_mod is in the model as the middle category of the same contrast
#   set but is not one of the six: the hypothesis is high versus low, so
#   ipaq_high carries the comparison and ipaq_mod is there so that the low
#   category remains the reference.
#
# INFERENCE
#   Two-tailed alpha = .05 with 95% confidence intervals, reported
#   unadjusted. No multiplicity correction is applied to this script's
#   output, because it is not a test family.
#
#   THE GAP STAYS IN YEARS, exactly as in 05_h2_single, so the two sets of
#   coefficients are on the same scale and can be tabulated side by side.
#
# SAMPLE AND VARIABLES
#   Identical to 05_h2_single: everyone with a gap for the modality in
#   question, same standardisation, same covariates. The two scripts must
#   describe the same participants or the comparison between them is
#   between two datasets rather than between two specifications.
#
# MEASUREMENT MODEL
#   The correlated two-factor psychosocial model from
#   02_cfa_psychosocial_factors.R, repeated verbatim. Well-being and
#   connectedness enter as latent variables so that their measurement
#   error is modelled rather than passed into the regression; the other
#   four exposures are observed.
#
# WHAT THIS SCRIPT READS
#   data/clean/mice_<TAG>_functional.rds   $completed
#   data/clean/mice_<TAG>_structural.rds   $completed
#   The $completed field, never $mids. See the contract in 03_impute.R.
#
# OUTPUT  results/sem/h2_joint_log_<TAG>.txt
#         results/sem/h2_joint_estimates_<TAG>.csv
#         results/sem/h2_joint_variance_<TAG>.csv
#         results/sem/h2_joint_variance_dropone_<TAG>.csv
#         results/sem/h2_joint_fit_<TAG>.csv
#         results/sem/h2_joint_provenance_<TAG>.txt
# =====================================================================

library(lavaan)
library(lavaan.mi)

TAG  <- Sys.getenv("FUNCTIONAL_TAG", "pca")
ROOT <- "C:/Users/inaye/Desktop/MASTERTHESIS"
IND  <- file.path(ROOT, "data", "clean")
RES  <- file.path(ROOT, "results", "sem")
dir.create(RES, recursive = TRUE, showWarnings = FALSE)

ALPHA <- 0.05
CONF  <- 1 - ALPHA          # 0.95

# Drop-one refits give each predictor's unique contribution to the gap.
# Six refits per modality on top of the two main models. Set FALSE to skip
# if the run is too slow; the block comparison in section 8 does not depend
# on it.
DROP_ONE <- TRUE

SMOKE_TEST <- TRUE

if (interactive())
  message("NOTE: run with Rscript or source(), not line by line.")

while (sink.number() > 0) sink()
try(sink(type = "message"), silent = TRUE)

con <- file(file.path(RES, sprintf("h2_joint_log_%s.txt", TAG)), open = "wt")
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
  normalizePath(file.path(ROOT, "code", "sem", "05b_h2_joint.R"),
                winslash = "/", mustWork = FALSE)
})

PROV <- list(
  run_time          = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  tag               = TAG,
  role              = "SENSITIVITY (mutually adjusted); 05_h2_single is the primary",
  alpha             = ALPHA,
  conf_level        = CONF,
  outcome_scale     = "years (gap left unstandardised, as pre-specified)",
  drop_one          = DROP_ONE,
  r_version         = R.version.string,
  lavaan_version    = as.character(utils::packageVersion("lavaan")),
  lavaan_mi_version = as.character(utils::packageVersion("lavaan.mi")),
  script_path       = this_script,
  script_md5        = if (file.exists(this_script))
                        unname(tools::md5sum(this_script)) else NA_character_)

cat("=== provenance ===\n")
for (k in names(PROV)) cat(sprintf("  %-18s %s\n", k, PROV[[k]]))
writeLines(vapply(names(PROV),
                  function(k) sprintf("%-18s %s", k, PROV[[k]]), character(1)),
           file.path(RES, sprintf("h2_joint_provenance_%s.txt", TAG)))


# ---------------------------------------------------------------------
# 2. Load the imputations and verify the contract
# ---------------------------------------------------------------------
PSY_IND <- c("positive_affect", "life_satisfaction", "meaning_purpose",
             "emotional_support", "instrumental_support", "friendship",
             "loneliness_rev")

load_imp <- function(modality) {
  f <- file.path(IND, sprintf("mice_%s_%s.rds", TAG, modality))
  if (!file.exists(f))
    stop("imputations not found: ", f, "\n  Run 03_impute.R first.")
  obj <- readRDS(f)

  for (fld in c("completed", "gap_var", "modality"))
    if (is.null(obj[[fld]]))
      stop(f, ": field '", fld, "' is absent. This file was not written by ",
           "the current 03_impute.R.")
  if (obj$modality != modality)
    stop(f, ": contains modality '", obj$modality, "', expected '", modality, "'")

  # The gap must still be missing where it was missing before imputation.
  # It is the outcome here, so an imputed value would be a fabricated
  # outcome, which is worse than a fabricated predictor.
  x <- obj$completed[[1]]
  if (!anyNA(x[[obj$gap_var]]))
    stop(f, ": ", obj$gap_var, " has no missing values in dataset 1. The ",
         "auxiliary imputations were not restored to NA. Do not use this file.")
  if (all(vapply(PSY_IND, function(v) sum(is.na(x[[v]])), integer(1)) == 0))
    stop(f, ": no psychosocial indicator has missing values in dataset 1. ",
         "This looks like $mids output rather than $completed.")

  cat(sprintf("\n%s: %d datasets, %d rows, %s has %d missing (correct)\n",
              basename(f), length(obj$completed), nrow(x), obj$gap_var,
              sum(is.na(x[[obj$gap_var]]))))
  obj
}

imp <- list(functional = load_imp("functional"),
            structural = load_imp("structural"))


# ---------------------------------------------------------------------
# 3. Analysis variables
# ---------------------------------------------------------------------
z <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) stop("cannot standardise a constant variable")
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

prepare <- function(obj) {
  gap <- obj$gap_var

  mask1 <- NULL
  out <- lapply(seq_along(obj$completed), function(i) {
    d <- obj$completed[[i]]
    # H2 sample: a usable gap, and nothing else. The cognitive-data
    # requirement belongs to the models that contain the cognitive factor.
    keep <- !is.na(d[[gap]])
    if (is.null(mask1)) mask1 <<- keep
    else if (!identical(keep, mask1))
      stop("dataset ", i, " (", obj$modality,
           "): analysis sample differs from dataset 1")
    d <- d[keep, , drop = FALSE]

    if (!all(as.character(d$sex) %in% c("F", "M")))
      stop("unexpected value in sex")

    data.frame(
      participant_id = d$participant_id,
      pedigree       = d$pedigree,
      # The psychosocial indicators are standardised for the same reason as
      # the cognitive ones in 04: it changes no fit index and no standardised
      # estimate, and puts scales with very different ranges onto one.
      setNames(as.data.frame(lapply(d[PSY_IND], z)), PSY_IND),
      # The outcome stays in years. This is the pre-specified exception.
      BAG       = as.numeric(d[[gap]]),
      income_z  = z(d$income_log),
      sleep_z   = z(d$sleep_quality),
      stress_z  = z(d$perceived_stress),
      ipaq_mod  = as.integer(d$ipaq_mod),
      ipaq_high = as.integer(d$ipaq_high),
      age_z     = z(d$age),
      edu_z     = z(as.numeric(d$education_num)),
      sex_male  = as.integer(as.character(d$sex) == "M"),
      stringsAsFactors = FALSE)
  })

  x <- out[[1]]
  cat(sprintf("\n%s H2 sample: n = %d\n", obj$modality, nrow(x)))
  cat(sprintf("  gap in years: M = %.2f, SD = %.2f, range [%.2f, %.2f]\n",
              mean(x$BAG), sd(x$BAG), min(x$BAG), max(x$BAG)))
  cat("  missing per psychosocial indicator (dataset 1, handled by FIML):\n")
  print(colSums(is.na(x[PSY_IND])))
  cat(sprintf("  activity: %d low (reference), %d moderate, %d high\n",
              sum(x$ipaq_mod == 0 & x$ipaq_high == 0),
              sum(x$ipaq_mod == 1), sum(x$ipaq_high == 1)))
  out
}

dat <- lapply(imp, prepare)

cat("\nNote on sample sizes: the H2 samples above are larger than the\n")
cat("cognitive-resilience samples used in 04 and 06, because H2 does not\n")
cat("involve the cognitive factor. This is intended, not an inconsistency.\n")


# ---------------------------------------------------------------------
# 4. Model
# ---------------------------------------------------------------------
# Measurement part: the correlated two-factor psychosocial model from
# 02_cfa_psychosocial_factors.R. Perceived stress is a single indicator
# and cannot support a latent factor, so it enters as an observed
# predictor; this is why five of the six predictors are observed and only
# two are latent.
#
# Structural part: the gap on all six predictors at once, plus the
# covariates. Simultaneous entry is what the analysis plan specifies, so
# each coefficient is the unique association of that predictor with the
# gap after the other five and the covariates. Well-being and
# connectedness correlate substantially at the latent level, so their
# simultaneous coefficients describe variance not shared with each other,
# and the supplementary single-predictor models in 07 exist to show the
# unadjusted associations alongside.

MEASUREMENT <- '
  PWB =~ positive_affect + life_satisfaction + meaning_purpose
  SOC =~ emotional_support + instrumental_support + friendship + loneliness_rev
  PWB ~~ SOC
'

SIX      <- c("income_z", "sleep_z", "ipaq_high", "PWB", "SOC", "stress_z")
COVARS   <- c("age_z", "sex_male", "edu_z")
IN_MODEL <- c("income_z", "sleep_z", "ipaq_mod", "ipaq_high",
              "PWB", "SOC", "stress_z")

model_h2 <- function(predictors) {
  rhs <- paste(c(predictors, COVARS), collapse = " + ")
  paste0(MEASUREMENT, "\n  BAG ~ ", rhs, "\n")
}

MODEL_FULL   <- model_h2(IN_MODEL)   # mutually adjusted
MODEL_COVAR  <- model_h2(character(0))  # covariates only, for the block test

cat("\n=== mutually adjusted model ===\n")
cat(MODEL_FULL)


# ---------------------------------------------------------------------
# 5. Fit
# ---------------------------------------------------------------------
fit_mi <- function(model, dlist, label, quiet = FALSE) {
  if (!quiet) cat("\n\n############ ", label, " ############\n")
  t0 <- Sys.time()
  fit <- lavaan.mi::sem.mi(model = model, data = dlist,
                           estimator = "MLR", missing = "fiml", std.lv = TRUE)
  if (!quiet)
    cat(sprintf("fitted %d imputations in %.1f s\n", length(dlist),
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  fit
}

if (SMOKE_TEST) {
  cat("\n=== smoke test: 2 imputations ===\n")
  tf <- fit_mi(MODEL_FULL, dat$structural[1:2], "smoke", quiet = TRUE)
  s  <- lavaan.mi::standardizedSolution.mi(tf)
  miss <- setdiff(SIX, s$rhs[s$op == "~" & s$lhs == "BAG"])
  if (length(miss))
    stop("smoke test: predictors absent from the fitted model: ",
         paste(miss, collapse = ", "))
  cat("  all six predictors present; proceeding to the full run\n")
}

fits <- list()
for (m in c("functional", "structural")) {
  fits[[m]] <- list(
    full  = fit_mi(MODEL_FULL,  dat[[m]], sprintf("H2 joint %s, mutually adjusted", m)),
    covar = fit_mi(MODEL_COVAR, dat[[m]], sprintf("H2 %s, covariates only", m)))
}


# ---------------------------------------------------------------------
# 6. Pooled estimates
# ---------------------------------------------------------------------
get_est <- function(fit, standardized) {
  f <- if (standardized) lavaan.mi::standardizedSolution.mi else
                         lavaan.mi::parameterEstimates.mi
  out <- tryCatch(f(fit, level = CONF), error = function(e) NULL)
  if (is.null(out)) {
    out <- f(fit)
    ecol <- if (standardized) "est.std" else "est"
    if (all(c(ecol, "se") %in% names(out))) {
      k <- stats::qnorm(1 - ALPHA / 2)
      out$ci.lower <- out[[ecol]] - k * out$se
      out$ci.upper <- out[[ecol]] + k * out$se
      cat("  note: level= not supported; interval rebuilt as est +/- ",
          round(k, 3), " * SE\n", sep = "")
    }
  }
  out
}

# R-squared for a variable from its standardised residual variance. For the
# gap this is the share of gap variance the model accounts for, which is the
# number the "brain age gap explains little" criticism is about.
r2_of <- function(fit, v) {
  s <- tryCatch(lavaan.mi::standardizedSolution.mi(fit), error = function(e) NULL)
  if (is.null(s)) return(NA_real_)
  row <- s$op == "~~" & s$lhs == v & s$rhs == v
  if (!any(row)) return(NA_real_)
  1 - s$est.std[row][1]
}

DIRECTION <- c(income_z = "-", sleep_z = "-", ipaq_high = "-",
               PWB = "-", SOC = "-", stress_z = "+")

rows <- list(); fitrows <- list(); varrows <- list()

# Fit indices to request, both flavours. See the note at the fitMeasures call.
FIT_WANT <- c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr",
              "cfi.robust", "tli.robust", "rmsea.robust")
robust_failed <- character(0)

for (m in names(fits)) {
  fit <- fits[[m]]$full
  cat("\n\n=== H2", m, ": mutually adjusted model ===\n")
  print(tryCatch(summary(fit, standardized = TRUE),
                 error = function(e) paste("summary() failed:", conditionMessage(e))))

  raw <- get_est(fit, FALSE)
  std <- get_est(fit, TRUE)

  # The gap is the outcome and stays in years, so the UNSTANDARDISED
  # coefficient is the one to read: years of gap per SD of predictor.
  b <- raw[raw$op == "~" & raw$lhs == "BAG", , drop = FALSE]
  s <- std[std$op == "~" & std$lhs == "BAG", , drop = FALSE]

  cat("\n--- gap in years per SD of predictor,", CONF * 100, "% CI ---\n")
  tab <- data.frame(
    predictor = b$rhs,
    years     = round(b$est, 3),
    se        = round(b$se, 3),
    ci.lower  = round(b$ci.lower, 3),
    ci.upper  = round(b$ci.upper, 3),
    p         = signif(b$pvalue, 3),
    beta_std  = round(s$est.std[match(b$rhs, s$rhs)], 3),
    is_exposure = ifelse(b$rhs %in% SIX, "yes", "covariate/contrast"),
    predicted = DIRECTION[b$rhs],
    row.names = NULL)
  tab$predicted[is.na(tab$predicted)] <- ""
  tab$as_predicted <- ifelse(tab$is_exposure != "yes", "",
                      ifelse(sign(tab$years) == ifelse(tab$predicted == "-", -1, 1),
                             "yes", "NO"))
  print(tab, row.names = FALSE)

  cat("\n  p values above are UNADJUSTED and stay that way. This script is\n")
  cat("  the sensitivity analysis; the H2 test family is the twelve rows\n")
  cat("  from 05_h2_single.R, and 08 corrects those. Read the estimates\n")
  cat("  here against the primary ones: shrinkage means shared variance.\n")

  b$modality <- m; b$type <- "unstandardized_years"
  s$modality <- m; s$type <- "standardized"
  b$is_exposure <- b$rhs %in% SIX
  s$is_exposure <- s$rhs %in% SIX
  rows[[paste(m, "raw")]] <- b
  rows[[paste(m, "std")]] <- s

  # Both flavours by name. Asking for "cfi" with pool.robust = TRUE is not
  # enough: when the robust computation fails lavaan WARNS and silently
  # returns the naive index, and a tryCatch on errors never sees it. Naming
  # cfi.robust makes a failure show up as NA rather than as a wrong number
  # under the right label.
  fm <- withCallingHandlers(
    tryCatch(
      lavaan::fitMeasures(fit, FIT_WANT,
                          pool.method = "D2", pool.robust = TRUE),
      error = function(e) tryCatch(lavaan::fitMeasures(fit, FIT_WANT),
                                   error = function(e2) NULL)),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("robust", msg, ignore.case = TRUE))
        robust_failed <<- union(robust_failed, trimws(msg))
      invokeRestart("muffleWarning")
    })
  if (!is.null(fm)) {
    cat("\n--- pooled fit ---\n"); print(round(fm, 3))
    fitrows[[m]] <- data.frame(modality = m, index = names(fm),
                               value = as.numeric(fm))
  }
}


# ---------------------------------------------------------------------
# 7. How much of the gap is accounted for
# ---------------------------------------------------------------------
# The direct answer to the standing criticism of brain-age research, put
# the other way round: given that the gap explains little cognitive
# variance, how much of the GAP itself do lifestyle and psychosocial
# factors explain?
#
# Three quantities per modality:
#   R2 with covariates only          age, sex, education
#   R2 with the six added            the mutually adjusted model
#   the difference                   what the six contribute as a block
#
# The block difference is the number to report. A single predictor's p
# value says whether an association is distinguishable from zero; the
# block R2 says whether the set of them accounts for an appreciable share
# of the gap. Those are different questions and a small p with a tiny
# block R2 is a real and reportable result, not a contradiction.

cat("\n\n=== variance in the gap accounted for ===\n")
cat(sprintf("%-12s %14s %14s %14s\n", "modality", "covariates", "full", "added by six"))
for (m in names(fits)) {
  r2c <- r2_of(fits[[m]]$covar, "BAG")
  r2f <- r2_of(fits[[m]]$full,  "BAG")
  cat(sprintf("%-12s %13.4f %14.4f %14.4f\n", m, r2c, r2f, r2f - r2c))
  varrows[[m]] <- data.frame(modality = m, r2_covariates = r2c,
                             r2_full = r2f, r2_added_by_six = r2f - r2c)
}
cat("\n  Read as proportions of gap variance. Report these next to the\n")
cat("  coefficients: they are what a reader who knows the brain-age\n")
cat("  literature will ask for, and giving the number first is stronger\n")
cat("  than being asked for it.\n")


# ---------------------------------------------------------------------
# 8. Unique contribution of each predictor
# ---------------------------------------------------------------------
# Drop-one refits. The unique contribution of a predictor is how much R2
# falls when it alone is removed from the full model. This is preferred
# here over squaring the standardised coefficient, because two of the six
# predictors are latent and have no observed column to regress on the
# others, and because the well-being and connectedness factors correlate
# substantially, which squaring would ignore.

if (DROP_ONE) {
  cat("\n\n=== unique contribution of each predictor (drop-one delta R2) ===\n")
  drop_rows <- list()
  for (m in names(fits)) {
    r2f <- r2_of(fits[[m]]$full, "BAG")
    cat(sprintf("\n%s, full model R2 = %.4f\n", m, r2f))
    for (p in SIX) {
      reduced <- setdiff(IN_MODEL, p)
      f <- fit_mi(model_h2(reduced), dat[[m]],
                  sprintf("%s minus %s", m, p), quiet = TRUE)
      r2r <- r2_of(f, "BAG")
      cat(sprintf("  without %-10s R2 = %.4f   unique = %.4f\n",
                  p, r2r, r2f - r2r))
      drop_rows[[paste(m, p)]] <-
        data.frame(modality = m, predictor = p, r2_full = r2f,
                   r2_without = r2r, unique_r2 = r2f - r2r)
    }
  }
  if (length(drop_rows))
    varrows[["drop_one"]] <- do.call(rbind, drop_rows)
} else {
  cat("\n\nDROP_ONE is FALSE; per-predictor contributions were not computed.\n")
}


# ---------------------------------------------------------------------
# 9. Save
# ---------------------------------------------------------------------
# rbind over a collected list must be guarded: Matrix, loaded elsewhere in
# the project, makes rbind an S4 generic whose dispatch indexes
# list(...)[[1]], so an empty list is an error rather than NULL. The two
# accessors also return different columns, so everything is forced onto one
# schema first.
CANON <- c("modality", "type", "is_exposure", "lhs", "op", "rhs",
           "est", "se", "z", "pvalue", "ci.lower", "ci.upper")

harmonise <- function(x) {
  if (!"est" %in% names(x) && "est.std" %in% names(x)) x$est <- x$est.std
  if (!"z"   %in% names(x) && "t"       %in% names(x)) x$z   <- x$t
  for (cc in setdiff(CANON, names(x))) x[[cc]] <- NA
  x[, CANON, drop = FALSE]
}

if (!length(rows)) stop("no estimates were collected; nothing to write")
est <- do.call(rbind, lapply(rows, harmonise))
write.csv(est, file.path(RES, sprintf("h2_joint_estimates_%s.csv", TAG)),
          row.names = FALSE)
cat("\nwrote", file.path(RES, sprintf("h2_joint_estimates_%s.csv", TAG)), "\n")

if (length(varrows)) {
  blocks <- varrows[names(varrows) != "drop_one"]
  if (length(blocks))
    write.csv(do.call(rbind, blocks),
              file.path(RES, sprintf("h2_joint_variance_%s.csv", TAG)), row.names = FALSE)
  if (!is.null(varrows$drop_one))
    write.csv(varrows$drop_one,
              file.path(RES, sprintf("h2_joint_variance_dropone_%s.csv", TAG)),
              row.names = FALSE)
  cat("wrote", file.path(RES, sprintf("h2_joint_variance_%s.csv", TAG)), "\n")
}

if (length(fitrows))
  write.csv(do.call(rbind, fitrows),
            file.path(RES, sprintf("h2_joint_fit_%s.csv", TAG)), row.names = FALSE)

if (length(robust_failed)) {
  cat("\n\n=== fit indices: robust versions unavailable ===\n")
  for (msg in robust_failed) cat("  ", msg, "\n", sep = "")
  cat("  The *.robust rows in h2_joint_fit are NA and the plain cfi/tli/rmsea\n")
  cat("  are the NAIVE indices, not the scaled ones. Label them as such.\n")
}

cat("\nH2 sensitivity complete. These are mutually adjusted estimates and\n")
cat("no multiplicity correction is applied to them. Tabulate them beside\n")
cat("h2_single_estimates_", TAG, ".csv; the difference between the two is\n", sep = "")
cat("the share of each association that is common to the other exposures.\n")

close_log()
