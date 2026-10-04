# ====================================================================
# H
# 3
# ,
#  
# S
# E
# N
# S
# I
# T
# I
# V
# I
# T
# Y
# .
#  
# A
# l
# l
#  
# s
# i
# x
#  
# e
# x
# p
# o
# s
# u
# r
# e
# s
#  
# p
# l
# u
# s
#  
# t
# h
# e
#  
# g
# a
# p
# ,
#  
# e
# n
# t
# e
# r
# e
# d
#  
# t
# o
# g
# e
# t
# h
# e
# r
# .
# ====================================================================

# =====================================================================
# H3, SENSITIVITY - Lifestyle and psychosocial correlates of cognitive
#                   resilience, all six exposures entered together
#
#     g ~ BAG + income + sleep + activity + PWB + SOC + stress
#         + age + sex + education
#
# Each coefficient is the association of that exposure with cognitive
# resilience after the other five AND after BAG, age, sex and education.
#
# THIS IS NOT THE PRIMARY ANALYSIS
#   The primary analysis is 06b_h3_joint.R, which fits one exposure per
#   model. This script answers a different question and is reported beside
#   it, not instead of it:
#
#     06_h3_single   is X associated with resilience?
#                    the covariate-adjusted total association
#
#     06b_h3_joint   is X associated with resilience over and above the
#                    other five exposures?
#                    the part of X not shared with them
#
#   The six are substantially intercorrelated, so a coefficient that
#   shrinks between the two is the expected result and is itself the
#   finding: the association is largely shared across psychosocial
#   domains. Shrinkage is not evidence that the primary association is
#   absent, and a sign flip under mutual adjustment is suppression.
#
#   Nothing here decides support for H3. The Benjamini-Hochberg family is
#   the twelve rows from 06_h3_single. These estimates are descriptive and
#   belong in the supplementary table.
#
# HOW TO RUN
#   Rscript code/sem/06b_h3_joint.R      (or source())
#
# DIRECTIONS, as in the primary (mirror of H2, because the outcome is
# cognition rather than the gap):
#   income +, sleep +, activity +, well-being +, connectedness +,
#   perceived stress -
#
#   ipaq_mod is in the model to keep low as the reference category but is
#   not one of the six.
#
# INFERENCE
#   Two-tailed alpha = .05 with 95% intervals, reported unadjusted. No
#   multiplicity correction is applied, because this is not a test family.
#
# SAMPLE, SCALE AND MEASUREMENT
#   Identical to 06_h3_single: same participants, same standardisation,
#   same bifactor measurement model, same REGRESS_SPECIFICS. The two
#   scripts must describe the same data or the comparison between them is
#   between two datasets rather than between two specifications.
#
# WHAT THIS SCRIPT READS
#   data/clean/mice_<TAG>_functional.rds   $completed
#   data/clean/mice_<TAG>_structural.rds   $completed
#
# OUTPUT  results/sem/h3_joint_log_<TAG>.txt
#         results/sem/h3_joint_estimates_<TAG>.csv
#         results/sem/h3_joint_variance_<TAG>.csv
#         results/sem/h3_joint_fit_<TAG>.txt is not written; see the csv
#         results/sem/h3_joint_provenance_<TAG>.txt
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
RESILIENCE_SCORES <- FALSE   # produced by the primary script

SMOKE_TEST <- TRUE

if (interactive())
  message("NOTE: run with Rscript or source(), not line by line.")

while (sink.number() > 0) sink()
try(sink(type = "message"), silent = TRUE)

con <- file(file.path(RES, sprintf("h3_joint_log_%s.txt", TAG)), open = "wt")
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
  normalizePath(file.path(ROOT, "code", "sem", "06b_h3_joint.R"),
                winslash = "/", mustWork = FALSE)
})

PROV <- list(
  run_time          = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  tag               = TAG,
  role              = "SENSITIVITY (mutually adjusted); 06_h3_single is the primary",
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
           file.path(RES, sprintf("h3_joint_provenance_%s.txt", TAG)))


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

# The baseline: the H1 model. What the six add as a BLOCK over this is the
# number to report from this script.
MODEL_BASE <- build(character(0))

# All six at once. ipaq_mod rides along to keep low as the reference.
ALL_TERMS <- c("income_z", "sleep_z", "ipaq_mod", "ipaq_high",
               "PWB", "SOC", "stress_z")
SIX       <- c("income_z", "sleep_z", "ipaq_high", "PWB", "SOC", "stress_z")
MODEL_FULL <- build(ALL_TERMS, paste0(MEAS_PWB, MEAS_SOC, "  PWB ~~ SOC\n"))

DIRECTION <- c(income_z = "+", sleep_z = "+", ipaq_high = "+",
               PWB = "+", SOC = "+", stress_z = "-")

cat("\n\n=== mutually adjusted model ===\n")
cat(MODEL_FULL)
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
  cat("\n=== smoke test: 2 imputations ===\n")
  tf <- fit_mi(MODEL_FULL, dat$structural[1:2], "smoke", quiet = TRUE)
  st <- lavaan.mi::standardizedSolution.mi(tf)
  present <- st$rhs[st$op == "~" & st$lhs == "g"]
  miss <- setdiff(c(ALL_TERMS, COVARS), present)
  if (length(miss))
    stop("smoke test: absent from the fitted model: ",
         paste(miss, collapse = ", "))
  cat("  all six exposures and the covariates are present\n")
  tf <- fit_mi(MODEL_BASE, dat$structural[1:2], "baseline", quiet = TRUE)
  cat("  baseline ok; proceeding to the full run\n")
}

fits <- list()
for (m in c("functional", "structural")) {
  fits[[m]] <- list(
    base = fit_mi(MODEL_BASE, dat[[m]],
                  sprintf("H3 joint %s, baseline (H1 model)", m)),
    full = fit_mi(MODEL_FULL, dat[[m]],
                  sprintf("H3 joint %s, mutually adjusted", m)))
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
  fit <- fits[[m]]$full
  r2b <- r2_of(fits[[m]]$base, "g")
  r2f <- r2_of(fit, "g")

  cat(sprintf("\n\n=== H3 joint %s: mutually adjusted ===\n", m))
  cat(sprintf("  R2(g), baseline = %.4f, with the six = %.4f, added = %.4f\n",
              r2b, r2f, r2f - r2b))

  std <- get_est(fit, TRUE)
  s_g <- std[std$op == "~" & std$lhs == "g", , drop = FALSE]

  tab <- data.frame(
    term      = s_g$rhs,
    beta      = round(s_g$est.std, 3),
    se        = round(s_g$se, 3),
    ci.lower  = round(s_g$ci.lower, 3),
    ci.upper  = round(s_g$ci.upper, 3),
    p         = signif(s_g$pvalue, 3),
    role      = ifelse(s_g$rhs %in% SIX, "exposure",
                ifelse(s_g$rhs == "ipaq_mod", "contrast", "covariate")),
    predicted = ifelse(s_g$rhs %in% names(DIRECTION),
                       DIRECTION[s_g$rhs], ""),
    row.names = NULL)
  tab$predicted[is.na(tab$predicted)] <- ""
  tab$as_predicted <- ifelse(tab$role != "exposure", "",
                      ifelse(sign(tab$beta) ==
                             ifelse(tab$predicted == "-", -1, 1), "yes", "NO"))

  cat("\n--- standardised beta on latent cognition,", CONF * 100, "% CI ---\n")
  print(tab, row.names = FALSE)

  cat("\n  p values above are UNADJUSTED and stay that way. This script is\n")
  cat("  the sensitivity analysis; the H3 test family is the twelve rows\n")
  cat("  from 06_h3_single.R, and 08 corrects those. Read these estimates\n")
  cat("  against the primary ones: shrinkage means shared variance.\n")

  # The BAG coefficient with all six in the model, against H1's value. A
  # large shift would mean the exposures are carrying brain-status variance.
  kb <- s_g$rhs == "BAG"
  sb <- get_est(fits[[m]]$base, TRUE)
  kk <- sb$op == "~" & sb$lhs == "g" & sb$rhs == "BAG"
  cat(sprintf("\n  BAG beta: %.3f at baseline, %.3f with the six exposures\n",
              sb$est.std[kk][1], s_g$est.std[kb][1]))

  std$modality <- m; std$type <- "standardized"
  std$is_exposure <- std$rhs %in% SIX
  raw <- get_est(fit, FALSE)
  raw$modality <- m; raw$type <- "unstandardized"
  raw$is_exposure <- raw$rhs %in% SIX
  rows[[paste(m, "std")]] <- std
  rows[[paste(m, "raw")]] <- raw

  varrows[[m]] <- data.frame(modality = m, r2_baseline = r2b,
                             r2_full = r2f, r2_added_by_six = r2f - r2b)

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
  if (!is.null(fm))
    fitrows[[m]] <- data.frame(modality = m, index = names(fm),
                               value = as.numeric(fm),
                               stringsAsFactors = FALSE)
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
write.csv(est, file.path(RES, sprintf("h3_joint_estimates_%s.csv", TAG)),
          row.names = FALSE)
cat("\nwrote", file.path(RES, sprintf("h3_joint_estimates_%s.csv", TAG)), "\n")

if (length(varrows)) {
  write.csv(do.call(rbind, varrows),
            file.path(RES, sprintf("h3_joint_variance_%s.csv", TAG)),
            row.names = FALSE)
  cat("wrote", file.path(RES, sprintf("h3_joint_variance_%s.csv", TAG)), "\n")
}
if (length(fitrows)) {
  write.csv(do.call(rbind, fitrows),
            file.path(RES, sprintf("h3_joint_fit_%s.csv", TAG)),
            row.names = FALSE)
  cat("wrote", file.path(RES, sprintf("h3_joint_fit_%s.csv", TAG)), "\n")
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

cat("\nH3 sensitivity complete. These are mutually adjusted estimates and\n")
cat("no multiplicity correction is applied to them. Tabulate them beside\n")
cat("h3_estimates_", TAG, ".csv; the difference between the two is the\n", sep = "")
cat("share of each association that is common to the other exposures.\n")

close_log()
