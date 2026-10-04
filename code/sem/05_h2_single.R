# ====================================================================
# H2, PRIMARY. What is associated with the brain-age gap itself?
# Six exposures, two modalities, one model each: twelve models. Separate models
# because the question is whether each exposure is associated with the gap, not
# what is left of it once the others are partialled out.
# The gap stays in years here, so coefficients read as years of brain age.
# ====================================================================

# =====================================================================
# H2, PRIMARY - Lifestyle and psychosocial correlates of the brain age gap
#                one exposure at a time
#
# Each of the six exposures is regressed on the modality-specific brain
# age gap in its OWN model, adjusted for chronological age, sex and
# education and for nothing else. Twelve models in total, six per
# modality.
#
#     BAG ~ X + age + sex + education
#
# WHY THIS IS THE PRIMARY SPECIFICATION
#   The question H2 asks is whether each exposure is associated with the
#   gap, and how large that association is. A separate model answers
#   exactly that: the covariate-adjusted total association, uncontaminated
#   by the choice of which other exposures happen to be in the model.
#
#   The six exposures are substantially intercorrelated - well-being,
#   social connectedness and low stress travel together - so entering
#   them simultaneously estimates something else: the part of each that
#   is not shared with the others. That is a real quantity and it is
#   worth reporting, but it is a second question, and a coefficient that
#   shrinks under mutual adjustment does not mean the association is
#   absent. It means the association is shared.
#
#   The mutually adjusted model is therefore fitted in 05b as a
#   sensitivity analysis, and the two sets of estimates are reported side
#   by side. Attenuation between them is the finding, not a problem.
#
# HOW TO RUN
#   Rscript code/sem/05_h2_single.R      (or source())
#
# Do not submit this file to the console line by line. Every stop() below
# is a guard meant to halt the run, and line-by-line submission lets
# execution continue past it.
#
# THE SIX EXPOSURES, ONE MODEL EACH
#   (i)   log household income          income_z    higher -> smaller gap
#   (ii)  sleep quality                 sleep_z     higher -> smaller gap
#   (iii) physical activity, high vs low ipaq_high  higher -> smaller gap
#   (iv)  psychological well-being      PWB latent  higher -> smaller gap
#   (v)   social connectedness          SOC latent  higher -> smaller gap
#   (vi)  perceived stress              stress_z    higher -> LARGER gap
#
#   Activity is the one exposure that needs two terms. It is a
#   three-level variable and the hypothesis is high versus low, so
#   ipaq_mod and ipaq_high enter together to keep low as the reference
#   and ipaq_high alone carries the test. Two terms, one exposure, one
#   test.
#
#   Well-being and connectedness stay latent. Their measurement model is
#   carried into the regression so that their indicators' measurement
#   error is modelled rather than averaged away. The measurement part is
#   included only in the models that need it; a model for income does not
#   estimate PWB.
#
# THE TEST FAMILY
#   Six exposures x two modalities = twelve tests. These twelve are the
#   H2 family. This script reports unadjusted p values and writes them
#   out; 08 applies Benjamini-Hochberg across all twelve and that is what
#   decides support for H2. Both modalities are kept in one family
#   because they are two measurements of the same hypothesis, not two
#   hypotheses.
#
# INFERENCE
#   Two-tailed alpha = .05 with 95% confidence intervals, before the BH
#   step in 08.
#
#   THE GAP STAYS IN YEARS. Continuous exposures are standardised but the
#   gap is left on its original scale wherever it is the outcome, so a
#   coefficient reads as years of brain-age gap per standard deviation of
#   the exposure. The standardised beta is reported beside it.
#
# EFFECT SIZE
#   For every exposure, the increment over a covariates-only model:
#
#       delta R2 = R2(covariates + X) - R2(covariates)
#
#   This is the number that answers "significant, but how much?", and it
#   is the one a reader who knows the brain-age literature will ask for.
#   Giving it unprompted is stronger than being asked.
#
# SAMPLE
#   Everyone with a gap for the modality in question. This is NOT the
#   cognitive-resilience sample used in 04 and 06: H2 does not involve
#   the cognitive factor, so the cognitive-data requirement does not
#   apply and excluding those participants would discard usable
#   information. The samples differ in size on purpose.
#
# WHAT THIS SCRIPT READS
#   data/clean/mice_<TAG>_functional.rds   $completed
#   data/clean/mice_<TAG>_structural.rds   $completed
#   The $completed field, never $mids. See the contract in 03_impute.R.
#
# OUTPUT  results/sem/h2_single_log_<TAG>.txt
#         results/sem/h2_single_estimates_<TAG>.csv   (feeds BH in 08)
#         results/sem/h2_single_variance_<TAG>.csv
#         results/sem/h2_single_fit_<TAG>.csv
#         results/sem/h2_single_provenance_<TAG>.txt
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

SMOKE_TEST <- TRUE

if (interactive())
  message("NOTE: run with Rscript or source(), not line by line.")

while (sink.number() > 0) sink()
try(sink(type = "message"), silent = TRUE)

con <- file(file.path(RES, sprintf("h2_single_log_%s.txt", TAG)), open = "wt")
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
  normalizePath(file.path(ROOT, "code", "sem", "05_h2_single.R"),
                winslash = "/", mustWork = FALSE)
})

PROV <- list(
  run_time          = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  tag               = TAG,
  role              = "PRIMARY (one exposure per model); 05b is the sensitivity",
  alpha             = ALPHA,
  conf_level        = CONF,
  outcome_scale     = "years (gap left unstandardised)",
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
           file.path(RES, sprintf("h2_single_provenance_%s.txt", TAG)))


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
# Identical to 05b, deliberately. The two scripts must describe the same
# participants on the same scales, or the comparison between them is
# between two datasets rather than between two specifications.

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
      setNames(as.data.frame(lapply(d[PSY_IND], z)), PSY_IND),
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
# 4. The twelve models
# ---------------------------------------------------------------------
# One entry per exposure. `terms` is what goes on the right-hand side,
# `test` is the single coefficient that carries the hypothesis, and
# `latent` says whether the measurement model has to come along.
#
# Activity is the only exposure whose terms and test differ: both
# contrasts are needed to define the reference category, but only the
# high-versus-low contrast is the test.

MEAS_PWB <- "  PWB =~ positive_affect + life_satisfaction + meaning_purpose\n"
MEAS_SOC <- paste0("  SOC =~ emotional_support + instrumental_support + ",
                   "friendship + loneliness_rev\n")

COVARS <- c("age_z", "sex_male", "edu_z")

EXPOSURES <- list(
  income   = list(terms = "income_z",               test = "income_z",
                  meas = "",         direction = "-",
                  label = "log household income"),
  sleep    = list(terms = "sleep_z",                test = "sleep_z",
                  meas = "",         direction = "-",
                  label = "sleep quality"),
  activity = list(terms = c("ipaq_mod", "ipaq_high"), test = "ipaq_high",
                  meas = "",         direction = "-",
                  label = "physical activity, high vs low"),
  pwb      = list(terms = "PWB",                    test = "PWB",
                  meas = MEAS_PWB,   direction = "-",
                  label = "psychological well-being (latent)"),
  soc      = list(terms = "SOC",                    test = "SOC",
                  meas = MEAS_SOC,   direction = "-",
                  label = "social connectedness (latent)"),
  stress   = list(terms = "stress_z",               test = "stress_z",
                  meas = "",         direction = "+",
                  label = "perceived stress"))

build <- function(e) {
  rhs <- paste(c(e$terms, COVARS), collapse = " + ")
  paste0(e$meas, "  BAG ~ ", rhs, "\n")
}

# The baseline every increment is measured against. Covariates only, no
# exposure, no measurement model.
MODEL_COVAR <- paste0("  BAG ~ ", paste(COVARS, collapse = " + "), "\n")

cat("\n\n=== the twelve primary models ===\n")
for (nm in names(EXPOSURES)) {
  cat("\n--", EXPOSURES[[nm]]$label, "--\n")
  cat(build(EXPOSURES[[nm]]))
}
cat("\n-- baseline for every increment --\n")
cat(MODEL_COVAR)


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

# A mistyped variable name is silently dropped by lavaan rather than
# raising an error, so each model is checked on two imputations before
# the full run commits to fifty.
if (SMOKE_TEST) {
  cat("\n=== smoke test: every model on 2 imputations ===\n")
  for (nm in names(EXPOSURES)) {
    e  <- EXPOSURES[[nm]]
    tf <- fit_mi(build(e), dat$structural[1:2], nm, quiet = TRUE)
    s  <- lavaan.mi::standardizedSolution.mi(tf)
    present <- s$rhs[s$op == "~" & s$lhs == "BAG"]
    miss <- setdiff(c(e$terms, COVARS), present)
    if (length(miss))
      stop("smoke test, ", nm, ": absent from the fitted model: ",
           paste(miss, collapse = ", "))
    cat(sprintf("  %-9s ok\n", nm))
  }
  tf <- fit_mi(MODEL_COVAR, dat$structural[1:2], "covariates", quiet = TRUE)
  cat("  baseline ok; proceeding to the full run\n")
}

fits <- list()
for (m in c("functional", "structural")) {
  fits[[m]] <- list(covar = fit_mi(MODEL_COVAR, dat[[m]],
                                   sprintf("H2 %s, covariates only", m)))
  for (nm in names(EXPOSURES))
    fits[[m]][[nm]] <- fit_mi(build(EXPOSURES[[nm]]), dat[[m]],
                              sprintf("H2 %s, %s", m, EXPOSURES[[nm]]$label))
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

# R-squared for a variable from its standardised residual variance.
r2_of <- function(fit, v) {
  s <- tryCatch(lavaan.mi::standardizedSolution.mi(fit), error = function(e) NULL)
  if (is.null(s)) return(NA_real_)
  row <- s$op == "~~" & s$lhs == v & s$rhs == v
  if (!any(row)) return(NA_real_)
  1 - s$est.std[row][1]
}

rows <- list(); fitrows <- list(); varrows <- list()

# Fit indices to request, both flavours. See the note at the fitMeasures call.
FIT_WANT <- c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr",
              "cfi.robust", "tli.robust", "rmsea.robust")
robust_failed <- character(0)

for (m in names(fits)) {
  r2c <- r2_of(fits[[m]]$covar, "BAG")
  cat(sprintf("\n\n=== H2 %s: twelve-family rows, covariates-only R2 = %.4f ===\n",
              m, r2c))

  tabs <- list()
  for (nm in names(EXPOSURES)) {
    e   <- EXPOSURES[[nm]]
    fit <- fits[[m]][[nm]]

    raw <- get_est(fit, FALSE)
    std <- get_est(fit, TRUE)
    b <- raw[raw$op == "~" & raw$lhs == "BAG", , drop = FALSE]
    s <- std[std$op == "~" & std$lhs == "BAG", , drop = FALSE]

    r2f <- r2_of(fit, "BAG")
    hit <- b$rhs == e$test
    if (!any(hit))
      stop(m, " / ", nm, ": the test coefficient ", e$test,
           " is not in the fitted model")

    tabs[[nm]] <- data.frame(
      exposure     = e$label,
      term         = e$test,
      years        = round(b$est[hit], 3),
      se           = round(b$se[hit], 3),
      ci.lower     = round(b$ci.lower[hit], 3),
      ci.upper     = round(b$ci.upper[hit], 3),
      p            = signif(b$pvalue[hit], 3),
      beta_std     = round(s$est.std[match(e$test, s$rhs)], 3),
      r2_model     = round(r2f, 4),
      delta_r2     = round(r2f - r2c, 4),
      predicted    = e$direction,
      as_predicted = ifelse(sign(b$est[hit]) ==
                            ifelse(e$direction == "-", -1, 1), "yes", "NO"),
      row.names = NULL)

    # Every coefficient from every model is written out, test or not, so
    # that the covariate estimates are on record too.
    b$modality <- m; b$type <- "unstandardized_years"; b$exposure <- nm
    s$modality <- m; s$type <- "standardized";         s$exposure <- nm
    b$is_test <- b$rhs == e$test
    s$is_test <- s$rhs == e$test
    rows[[paste(m, nm, "raw")]] <- b
    rows[[paste(m, nm, "std")]] <- s

    varrows[[paste(m, nm)]] <- data.frame(
      modality = m, exposure = nm, term = e$test,
      r2_covariates = r2c, r2_model = r2f, delta_r2 = r2f - r2c)

    # Both flavours by name. Asking for "cfi" with pool.robust = TRUE is not
    # enough: when the robust computation fails lavaan WARNS and silently
    # returns the naive index, and a tryCatch on errors never sees it. Then
    # this script reports naive indices while 01 and 02 report robust ones,
    # with nothing saying which is which. Naming cfi.robust makes a failure
    # show up as NA rather than as a wrong number under the right label.
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
    # Long format, as in 04 and 05b. Transposing the named vector into wide
    # columns recycled the values across columns and produced an unreadable
    # table; index/value cannot do that.
    if (!is.null(fm))
      fitrows[[paste(m, nm)]] <-
        data.frame(modality = m, exposure = nm,
                   index = names(fm), value = as.numeric(fm),
                   stringsAsFactors = FALSE)
  }

  cat("\n--- gap in years per SD of exposure,", CONF * 100, "% CI ---\n")
  print(do.call(rbind, tabs), row.names = FALSE)

  cat("\n  Six models, six rows. Each row comes from its own model and is\n")
  cat("  adjusted for age, sex and education only. delta_r2 is the share of\n")
  cat("  gap variance the exposure adds over those covariates.\n")
  cat("  p values are UNADJUSTED; 08 applies Benjamini-Hochberg across the\n")
  cat("  twelve rows of this family and that decides support for H2.\n")
}

if (length(robust_failed)) {
  cat("\n\n=== fit indices: robust versions unavailable ===\n")
  for (msg in robust_failed) cat("  ", msg, "\n", sep = "")
  cat("  The *.robust columns in h2_single_fit are NA and the plain cfi/tli/\n")
  cat("  rmsea are the NAIVE indices, not the scaled ones. Label them as such;\n")
  cat("  01 and 02 report robust indices and the two do not belong in the\n")
  cat("  same column of a table without a note.\n")
}


# ---------------------------------------------------------------------
# 7. Save
# ---------------------------------------------------------------------
# rbind over a collected list must be guarded: Matrix, loaded elsewhere in
# the project, makes rbind an S4 generic whose dispatch indexes
# list(...)[[1]], so an empty list is an error rather than NULL. The two
# accessors also return different columns, so everything is forced onto
# one schema first.
CANON <- c("modality", "exposure", "type", "is_test", "lhs", "op", "rhs",
           "est", "se", "z", "pvalue", "ci.lower", "ci.upper")

harmonise <- function(x) {
  if (!"est" %in% names(x) && "est.std" %in% names(x)) x$est <- x$est.std
  if (!"z"   %in% names(x) && "t"       %in% names(x)) x$z   <- x$t
  for (cc in setdiff(CANON, names(x))) x[[cc]] <- NA
  x[, CANON, drop = FALSE]
}

if (!length(rows)) stop("no estimates were collected; nothing to write")
est <- do.call(rbind, lapply(rows, harmonise))
write.csv(est, file.path(RES, sprintf("h2_single_estimates_%s.csv", TAG)),
          row.names = FALSE)
cat("\nwrote", file.path(RES, sprintf("h2_single_estimates_%s.csv", TAG)), "\n")

if (length(varrows)) {
  write.csv(do.call(rbind, varrows),
            file.path(RES, sprintf("h2_single_variance_%s.csv", TAG)),
            row.names = FALSE)
  cat("wrote", file.path(RES, sprintf("h2_single_variance_%s.csv", TAG)), "\n")
}

if (length(fitrows)) {
  write.csv(do.call(rbind, fitrows),
            file.path(RES, sprintf("h2_single_fit_%s.csv", TAG)),
            row.names = FALSE)
  cat("wrote", file.path(RES, sprintf("h2_single_fit_%s.csv", TAG)), "\n")
}

cat("\nH2 primary complete. The twelve family rows are those with\n")
cat("is_test = TRUE in h2_single_estimates_", TAG, ".csv. They carry\n", sep = "")
cat("unadjusted p values; 08 applies Benjamini-Hochberg across the family.\n")
cat("Run 05b_h2_joint.R next for the mutually adjusted sensitivity\n")
cat("analysis, and report the two side by side.\n")

close_log()
