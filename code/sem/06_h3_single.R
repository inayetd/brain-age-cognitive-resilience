# ====================================================================
# H3, PRIMARY. What is associated with cognitive resilience?
# The same six exposures, now predicting cognition with the gap already in the
# model. By the Frisch-Waugh-Lovell theorem the exposure coefficient is exactly its
# association with the part of cognition the gap and the demographics do not
# explain - which is what this thesis calls resilience. No separate residual score
# is computed or used, because a two-step residual would understate the standard
# errors.
# The gap is standardised here, unlike in H2.
# ====================================================================

# =====================================================================
# H3, PRIMARY - Lifestyle and psychosocial correlates of cognitive
#                resilience, one exposure at a time
#
# Cognitive resilience is latent cognitive performance conditional on the
# brain age gap, chronological age, sex and education. Each of the six
# exposures is tested in its OWN model, twelve models in total:
#
#     g ~ BAG + X + age + sex + education
#
# WHY THE EXPOSURE COEFFICIENT IS THE RESILIENCE ASSOCIATION
#   There is no separate residual score in this model and none is needed.
#   By the Frisch-Waugh-Lovell theorem the coefficient on X in the model
#   above is identical to the coefficient from regressing the residual of
#   g on the residual of X, both residualised on BAG, age, sex and
#   education. So the X coefficient already IS the association between X
#   and the part of cognition that the brain age gap and the demographics
#   do not account for. That is the definition of resilience used here.
#
#   The two-step alternative - save a resilience score, then regress it on
#   raw X - is not the same quantity and is worse in two ways. A score
#   treated as data is treated as error-free, so the second-step standard
#   errors do not carry the uncertainty of the first step; and regressing
#   the residual of g on an unresidualised X does not return the partial
#   coefficient. With income and education correlated as strongly as they
#   are here, that difference is not academic.
#
#   A descriptive resilience score IS produced, in section 8, for figures
#   and for describing who the high-resilience participants are. It is
#   never used in an inferential model.
#
# WHY ONE EXPOSURE PER MODEL
#   The question H3 asks is whether each exposure is associated with
#   resilience, and how large that association is. A separate model answers
#   that directly, free of the choice of which other exposures happen to
#   be in the model. The six are substantially intercorrelated, so
#   entering them simultaneously estimates a different thing: the part of
#   each that is not shared with the others. That is the sensitivity
#   analysis in 06b, and the two are reported side by side.
#
# HOW TO RUN
#   Rscript code/sem/06_h3_single.R      (or source())
#
# Do not submit this file to the console line by line. Every stop() below
# is a guard meant to halt the run, and line-by-line submission lets
# execution continue past it.
#
# THE SIX EXPOSURES, ONE MODEL EACH
#   Directions are the MIRROR of H2, because the outcome is cognition
#   rather than the gap:
#
#   (i)   log household income           higher -> higher resilience  (+)
#   (ii)  sleep quality                  higher -> higher resilience  (+)
#   (iii) physical activity, high vs low higher -> higher resilience  (+)
#   (iv)  psychological well-being       higher -> higher resilience  (+)
#   (v)   social connectedness           higher -> higher resilience  (+)
#   (vi)  perceived stress               higher -> LOWER  resilience  (-)
#
#   Activity needs two terms: it is a three-level variable and the
#   hypothesis is high versus low, so ipaq_mod and ipaq_high enter
#   together to keep low as the reference, and ipaq_high alone carries the
#   test. Two terms, one exposure, one test.
#
#   Well-being and connectedness stay latent, and their measurement model
#   is carried only into the models that need it.
#
# THE TEST FAMILY
#   Six exposures x two modalities = twelve tests. 08_multiplicity.R
#   applies Benjamini-Hochberg across all twelve and that decides support
#   for H3. This script reports unadjusted p values and writes them out.
#
# THE MEASUREMENT MODEL AND THE SPECIFIC FACTORS
#   The cognitive measurement model is the bifactor-(S-1) from 01,
#   reproduced verbatim. REGRESS_SPECIFICS matches 04: when TRUE the two
#   specific factors are regressed on the same right-hand side as g, so
#   that the g coefficient is not absorbing variance that belongs to Exec
#   or Mem. It is set to the same value as in 04 so that H1 and H3 are
#   nested specifications of one model rather than two unrelated ones.
#
# SCALE
#   The outcome is latent and standardised, and BAG is standardised, so
#   the reported coefficients are standardised throughout. This differs
#   from H2, where the gap is the outcome and stays in years. The two sets
#   of coefficients are therefore NOT comparable in size, and the thesis
#   should not table them together without saying so.
#
# EFFECT SIZE
#   For every exposure, the increment over the H1 model:
#
#       delta R2 = R2(BAG + covariates + X) - R2(BAG + covariates)
#
#   This is the share of cognitive variance the exposure adds beyond brain
#   status and demographics, which is the quantity the resilience claim
#   rests on.
#
# SAMPLE
#   sem_eligible and a usable gap for the modality: the same sample as H1,
#   n = 1342 in both modalities. This is SMALLER than the H2 sample, which
#   has no cognitive requirement. The difference is intended.
#
# WHAT THIS SCRIPT READS
#   data/clean/mice_<TAG>_functional.rds   $completed
#   data/clean/mice_<TAG>_structural.rds   $completed
#   The $completed field, never $mids. See the contract in 03_impute.R.
#
# OUTPUT  results/sem/h3_single_log_<TAG>.txt
#         results/sem/h3_estimates_<TAG>.csv          (feeds BH in 08)
#         results/sem/h3_single_variance_<TAG>.csv
#         results/sem/h3_single_fit_<TAG>.csv
#         results/sem/h3_resilience_scores_<TAG>.csv  (DESCRIPTIVE ONLY)
#         results/sem/h3_single_provenance_<TAG>.txt
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

# Must match 04_h1_base.R. H1 and H3 are nested specifications of the same
# model and a mismatch here would make them incomparable.
REGRESS_SPECIFICS <- TRUE

# Section 8. Costs one extra lavaan fit per imputation per modality and
# produces nothing any inferential model reads.
RESILIENCE_SCORES <- TRUE

SMOKE_TEST <- TRUE

if (interactive())
  message("NOTE: run with Rscript or source(), not line by line.")

while (sink.number() > 0) sink()
try(sink(type = "message"), silent = TRUE)

con <- file(file.path(RES, sprintf("h3_single_log_%s.txt", TAG)), open = "wt")
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
  normalizePath(file.path(ROOT, "code", "sem", "06_h3_single.R"),
                winslash = "/", mustWork = FALSE)
})

PROV <- list(
  run_time          = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  tag               = TAG,
  role              = "PRIMARY (one exposure per model); 06b is the sensitivity",
  alpha             = ALPHA,
  conf_level        = CONF,
  regress_specifics = REGRESS_SPECIFICS,
  outcome_scale     = "standardised latent cognition",
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
           file.path(RES, sprintf("h3_single_provenance_%s.txt", TAG)))


# ---------------------------------------------------------------------
# 2. Load the imputations and verify the contract
# ---------------------------------------------------------------------
IND_COLS <- c("flanker", "card_sort", "list_sorting", "pattern_comp",
              "picture_seq", "ravlt", "trails_b_rev")
PSY_IND  <- c("positive_affect", "life_satisfaction", "meaning_purpose",
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

  x <- obj$completed[[1]]
  if (!anyNA(x[[obj$gap_var]]))
    stop(f, ": ", obj$gap_var, " has no missing values in dataset 1. The ",
         "auxiliary imputations were not restored to NA. Do not use this file.")
  if (all(vapply(IND_COLS, function(v) sum(is.na(x[[v]])), integer(1)) == 0))
    stop(f, ": no cognitive indicator has missing values in dataset 1. ",
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
    # The H1 sample, exactly: cognitive eligibility and a usable gap.
    keep <- d$sem_eligible & !is.na(d[[gap]])
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
      setNames(as.data.frame(lapply(d[IND_COLS], z)), IND_COLS),
      setNames(as.data.frame(lapply(d[PSY_IND],  z)), PSY_IND),
      BAG       = z(d[[gap]]),
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
  cat(sprintf("\n%s H3 sample: n = %d\n", obj$modality, nrow(x)))
  cat("  missing per cognitive indicator (dataset 1, handled by FIML):\n")
  print(colSums(is.na(x[IND_COLS])))
  cat("  missing per psychosocial indicator (dataset 1, handled by FIML):\n")
  print(colSums(is.na(x[PSY_IND])))
  cat(sprintf("  activity: %d low (reference), %d moderate, %d high\n",
              sum(x$ipaq_mod == 0 & x$ipaq_high == 0),
              sum(x$ipaq_mod == 1), sum(x$ipaq_high == 1)))
  out
}

dat <- lapply(imp, prepare)

cat("\nNote on sample sizes: the H3 samples match H1 and are SMALLER than\n")
cat("the H2 samples, because H3 contains the cognitive factor and H2 does\n")
cat("not. This is intended, not an inconsistency.\n")


# ---------------------------------------------------------------------
# 4. The twelve models
# ---------------------------------------------------------------------
MEASUREMENT <- '
  g    =~ flanker + card_sort + list_sorting + pattern_comp +
          picture_seq + ravlt + trails_b_rev
  Exec =~ a*flanker + a*card_sort
  Mem  =~ picture_seq + ravlt + list_sorting

  g    ~~ 0*Exec
  g    ~~ 0*Mem
  Exec ~~ 0*Mem
'
MEAS_PWB <- "  PWB =~ positive_affect + life_satisfaction + meaning_purpose\n"
MEAS_SOC <- paste0("  SOC =~ emotional_support + instrumental_support + ",
                   "friendship + loneliness_rev\n")

COVARS <- c("BAG", "age_z", "sex_male", "edu_z")

EXPOSURES <- list(
  income   = list(terms = "income_z",                test = "income_z",
                  meas = "", direction = "+",
                  label = "log household income"),
  sleep    = list(terms = "sleep_z",                 test = "sleep_z",
                  meas = "", direction = "+",
                  label = "sleep quality"),
  activity = list(terms = c("ipaq_mod", "ipaq_high"), test = "ipaq_high",
                  meas = "", direction = "+",
                  label = "physical activity, high vs low"),
  pwb      = list(terms = "PWB",                     test = "PWB",
                  meas = MEAS_PWB, direction = "+",
                  label = "psychological well-being (latent)"),
  soc      = list(terms = "SOC",                     test = "SOC",
                  meas = MEAS_SOC, direction = "+",
                  label = "social connectedness (latent)"),
  stress   = list(terms = "stress_z",                test = "stress_z",
                  meas = "", direction = "-",
                  label = "perceived stress"))

# The specific factors are regressed on the same right-hand side as g, so
# that the g coefficient is not quietly absorbing Exec or Mem variance.
build <- function(terms, meas = "") {
  rhs <- paste(c(terms, COVARS), collapse = " + ")
  s <- paste0(MEASUREMENT, meas, "\n  g ~ ", rhs, "\n")
  if (REGRESS_SPECIFICS)
    s <- paste0(s, "  Exec ~ ", rhs, "\n  Mem  ~ ", rhs, "\n")
  s
}

# The baseline every increment is measured against: the H1 model. What each
# exposure adds over this is the resilience signal.
MODEL_BASE <- build(character(0))

cat("\n\n=== the twelve primary models ===\n")
for (nm in names(EXPOSURES)) {
  cat("\n--", EXPOSURES[[nm]]$label, "--\n")
  cat(build(EXPOSURES[[nm]]$terms, EXPOSURES[[nm]]$meas))
}
cat("\n-- baseline (the H1 model) --\n")
cat(MODEL_BASE)


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

# lavaan drops a mistyped variable silently rather than raising an error, so
# every model is checked on two imputations before fifty are committed.
if (SMOKE_TEST) {
  cat("\n=== smoke test: every model on 2 imputations ===\n")
  for (nm in names(EXPOSURES)) {
    e  <- EXPOSURES[[nm]]
    tf <- fit_mi(build(e$terms, e$meas), dat$structural[1:2], nm, quiet = TRUE)
    s  <- lavaan.mi::standardizedSolution.mi(tf)
    present <- s$rhs[s$op == "~" & s$lhs == "g"]
    miss <- setdiff(c(e$terms, COVARS), present)
    if (length(miss))
      stop("smoke test, ", nm, ": absent from the fitted model: ",
           paste(miss, collapse = ", "))
    cat(sprintf("  %-9s ok\n", nm))
  }
  tf <- fit_mi(MODEL_BASE, dat$structural[1:2], "baseline", quiet = TRUE)
  cat("  baseline ok; proceeding to the full run\n")
}

fits <- list()
for (m in c("functional", "structural")) {
  fits[[m]] <- list(base = fit_mi(MODEL_BASE, dat[[m]],
                                  sprintf("H3 %s, baseline (H1 model)", m)))
  for (nm in names(EXPOSURES))
    fits[[m]][[nm]] <- fit_mi(build(EXPOSURES[[nm]]$terms, EXPOSURES[[nm]]$meas),
                              dat[[m]],
                              sprintf("H3 %s, %s", m, EXPOSURES[[nm]]$label))
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

r2_of <- function(fit, v) {
  s <- tryCatch(lavaan.mi::standardizedSolution.mi(fit), error = function(e) NULL)
  if (is.null(s)) return(NA_real_)
  row <- s$op == "~~" & s$lhs == v & s$rhs == v
  if (!any(row)) return(NA_real_)
  1 - s$est.std[row][1]
}

rows <- list(); fitrows <- list(); varrows <- list()

# Both flavours by name. Asking for "cfi" with pool.robust = TRUE is not
# enough: when the robust computation fails lavaan WARNS and silently
# returns the naive index, and tryCatch on errors never sees it.
FIT_WANT <- c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr",
              "cfi.robust", "tli.robust", "rmsea.robust")
robust_failed <- character(0)

for (m in names(fits)) {
  r2b <- r2_of(fits[[m]]$base, "g")
  cat(sprintf("\n\n=== H3 %s: twelve-family rows, baseline R2(g) = %.4f ===\n",
              m, r2b))

  tabs <- list()
  for (nm in names(EXPOSURES)) {
    e   <- EXPOSURES[[nm]]
    fit <- fits[[m]][[nm]]

    std <- get_est(fit, TRUE)
    raw <- get_est(fit, FALSE)
    s <- std[std$op == "~" & std$lhs == "g", , drop = FALSE]
    b <- raw[raw$op == "~" & raw$lhs == "g", , drop = FALSE]

    r2f <- r2_of(fit, "g")
    hit <- s$rhs == e$test
    if (!any(hit))
      stop(m, " / ", nm, ": the test coefficient ", e$test,
           " is not in the fitted model")

    tabs[[nm]] <- data.frame(
      exposure     = e$label,
      term         = e$test,
      beta         = round(s$est.std[hit], 3),
      se           = round(s$se[hit], 3),
      ci.lower     = round(s$ci.lower[hit], 3),
      ci.upper     = round(s$ci.upper[hit], 3),
      p            = signif(s$pvalue[hit], 3),
      r2_model     = round(r2f, 4),
      delta_r2     = round(r2f - r2b, 4),
      predicted    = e$direction,
      as_predicted = ifelse(sign(s$est.std[hit]) ==
                            ifelse(e$direction == "-", -1, 1), "yes", "NO"),
      row.names = NULL)

    # Everything is written out, test rows and covariates alike, so that the
    # BAG coefficient in each H3 model is on record next to H1's.
    s$modality <- m; s$type <- "standardized";         s$exposure <- nm
    b$modality <- m; b$type <- "unstandardized";       b$exposure <- nm
    s$is_test <- s$rhs == e$test
    b$is_test <- b$rhs == e$test
    rows[[paste(m, nm, "std")]] <- s
    rows[[paste(m, nm, "raw")]] <- b

    varrows[[paste(m, nm)]] <- data.frame(
      modality = m, exposure = nm, term = e$test,
      r2_baseline = r2b, r2_model = r2f, delta_r2 = r2f - r2b)

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

  cat("\n--- standardised beta on latent cognition,", CONF * 100, "% CI ---\n")
  print(do.call(rbind, tabs), row.names = FALSE)

  cat("\n  Each row comes from its own model, adjusted for BAG, age, sex and\n")
  cat("  education. The coefficient is the association with the part of\n")
  cat("  cognition those four do not account for, which is resilience.\n")
  cat("  delta_r2 is what the exposure adds over the H1 model.\n")
  cat("  p values are UNADJUSTED; 08 applies Benjamini-Hochberg across the\n")
  cat("  twelve rows of this family and that decides support for H3.\n")

  # The BAG coefficient should stay close to its H1 value across the six
  # models. A large shift would mean an exposure is carrying brain-status
  # variance, which is worth seeing rather than discovering later.
  cat("\n--- BAG coefficient across the six models (should be stable) ---\n")
  bag <- do.call(rbind, lapply(names(EXPOSURES), function(nm) {
    s <- get_est(fits[[m]][[nm]], TRUE)
    k <- s$op == "~" & s$lhs == "g" & s$rhs == "BAG"
    data.frame(model = nm, beta_BAG = round(s$est.std[k][1], 3))
  }))
  sb <- get_est(fits[[m]]$base, TRUE)
  kb <- sb$op == "~" & sb$lhs == "g" & sb$rhs == "BAG"
  bag <- rbind(data.frame(model = "baseline (H1)",
                          beta_BAG = round(sb$est.std[kb][1], 3)), bag)
  print(bag, row.names = FALSE)
}

if (length(robust_failed)) {
  cat("\n\n=== fit indices: robust versions unavailable ===\n")
  for (msg in robust_failed) cat("  ", msg, "\n", sep = "")
  cat("  The *.robust columns are NA and the plain cfi/tli/rmsea are the\n")
  cat("  NAIVE indices, not the scaled ones. Label them as such.\n")
}


# ---------------------------------------------------------------------
# 7. Save the estimates
# ---------------------------------------------------------------------
# Matrix, loaded elsewhere in this project, makes rbind an S4 generic whose
# dispatch indexes list(...)[[1]], so an empty list errors rather than
# returning NULL. The two accessors also return different columns, so
# everything is forced onto one schema first.
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
write.csv(est, file.path(RES, sprintf("h3_estimates_%s.csv", TAG)),
          row.names = FALSE)
cat("\nwrote", file.path(RES, sprintf("h3_estimates_%s.csv", TAG)), "\n")

if (length(varrows)) {
  write.csv(do.call(rbind, varrows),
            file.path(RES, sprintf("h3_single_variance_%s.csv", TAG)),
            row.names = FALSE)
  cat("wrote", file.path(RES, sprintf("h3_single_variance_%s.csv", TAG)), "\n")
}
if (length(fitrows)) {
  write.csv(do.call(rbind, fitrows),
            file.path(RES, sprintf("h3_single_fit_%s.csv", TAG)),
            row.names = FALSE)
  cat("wrote", file.path(RES, sprintf("h3_single_fit_%s.csv", TAG)), "\n")
}


# ---------------------------------------------------------------------
# 8. A descriptive resilience score
# ---------------------------------------------------------------------
# DESCRIPTIVE ONLY. Nothing above uses this, and nothing downstream should.
# It exists so that resilience can be plotted, and so that the question
# "who are the high-resilience participants?" can be answered.
#
# Why it must not enter an inferential model: a factor score is not the
# factor. It carries measurement error and factor-score indeterminacy, and
# a model that takes it as data treats it as error-free, which makes the
# standard errors too small. The confirmatory tests above avoid this by
# keeping cognition latent.
#
# Construction, per imputed dataset:
#   fit the baseline model with plain lavaan,
#   take the g factor score,
#   subtract the part predicted by BAG, age, sex and education.
# The remainder is cognition higher or lower than brain status and
# demographics would imply, in SD units of g.
#
# Averaged over the 50 datasets. The SD across datasets is kept beside the
# mean, for the same reason the brain age gap carries a within-participant
# SD: a score whose value depends on which imputation you drew is a score
# to be careful with.

if (RESILIENCE_SCORES) {
  cat("\n\n=== descriptive resilience scores ===\n")
  score_rows <- list()

  for (m in names(dat)) {
    dl <- dat[[m]]
    acc <- NULL
    for (i in seq_along(dl)) {
      f <- tryCatch(
        lavaan::sem(MODEL_BASE, data = dl[[i]], estimator = "MLR",
                    missing = "fiml", std.lv = TRUE),
        error = function(e) NULL)
      if (is.null(f) || !lavInspect(f, "converged")) next

      fs  <- lavPredict(f, method = "regression")
      idx <- lavInspect(f, "case.idx")
      gsc <- rep(NA_real_, nrow(dl[[i]]))
      gsc[idx] <- fs[, "g"]

      pe <- lavaan::parameterEstimates(f)
      bb <- pe$est[pe$op == "~" & pe$lhs == "g"]
      names(bb) <- pe$rhs[pe$op == "~" & pe$lhs == "g"]
      fitted_g <- rowSums(sapply(names(bb),
                                 function(v) bb[[v]] * dl[[i]][[v]]))
      r <- gsc - fitted_g

      acc <- if (is.null(acc)) matrix(r, ncol = 1) else cbind(acc, r)
    }

    if (is.null(acc)) {
      cat(sprintf("  %s: no imputation produced a usable score\n", m))
      next
    }
    cat(sprintf("  %s: scores from %d of %d imputations\n",
                m, ncol(acc), length(dl)))

    score_rows[[m]] <- data.frame(
      participant_id = dl[[1]]$participant_id,
      modality       = m,
      resilience     = rowMeans(acc, na.rm = TRUE),
      resilience_sd  = apply(acc, 1, stats::sd, na.rm = TRUE),
      n_imputations  = ncol(acc),
      stringsAsFactors = FALSE)

    v <- score_rows[[m]]$resilience
    cat(sprintf("    M = %.3f, SD = %.3f, range [%.2f, %.2f]\n",
                mean(v, na.rm = TRUE), sd(v, na.rm = TRUE),
                min(v, na.rm = TRUE), max(v, na.rm = TRUE)))
    cat(sprintf("    mean within-participant SD across imputations: %.3f\n",
                mean(score_rows[[m]]$resilience_sd, na.rm = TRUE)))
  }

  if (length(score_rows)) {
    p <- file.path(RES, sprintf("h3_resilience_scores_%s.csv", TAG))
    write.csv(do.call(rbind, score_rows), p, row.names = FALSE)
    cat("\nwrote", p, "\n")
    cat("DESCRIPTIVE ONLY. Do not regress anything on this column; the\n")
    cat("confirmatory answer to that question is in h3_estimates.\n")
  }
}

cat("\nH3 primary complete. The twelve family rows are those with\n")
cat("is_test = TRUE and type = standardized in h3_estimates_", TAG, ".csv.\n",
    sep = "")
cat("They carry unadjusted p values; 08_multiplicity.R applies\n")
cat("Benjamini-Hochberg across the family. Run 06b_h3_joint.R for the\n")
cat("mutually adjusted sensitivity analysis.\n")

close_log()
