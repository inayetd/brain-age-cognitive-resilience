# ====================================================================
# H1. Does the brain-age gap relate to cognition at all?
# Latent general cognition regressed on the gap plus age, sex and education, fitted
# separately per modality. The executive and memory factors are regressed on the
# same predictors so that the specificity of the association is observed rather
# than assumed.
# If H1 were null, H3 would have nothing to be a residual of.
# ====================================================================

# =====================================================================
# H1 - Brain age gap and latent cognitive performance
#
# The base structural model. Latent cognitive performance, the general
# factor of the bifactor-(S-1) measurement model, is regressed on the
# modality-specific brain age gap and the pre-specified covariates
# (chronological age, sex, education). Fitted separately for the
# functional and the structural gap.
#
# HOW TO RUN
#   Rscript code/sem/04_h1_base.R      (or source() from the console)
#
# Do not submit this file to the console line by line. Every stop() below
# is a guard meant to halt the run, and line-by-line submission lets
# execution continue past it.
#
# INFERENCE CRITERIA, AS PRE-SPECIFIED
#   The two modality-specific H1 tests each use alpha = .025, a Bonferroni
#   correction across the functional and structural modalities, and are
#   therefore reported with 97.5% confidence intervals. This differs from
#   the two-tailed alpha = .05 and 95% intervals used elsewhere in the
#   protocol, so the level is set once in ALPHA below and used everywhere.
#
#   Standardised coefficients are reported. Continuous predictors are
#   standardised before modelling; sex enters as an indicator variable
#   with female as the reference category. The gap is a predictor here,
#   so it is standardised too. It stays in years only in H2, where it is
#   the outcome.
#
# WHAT THIS SCRIPT READS
#   data/clean/mice_<TAG>_functional.rds   $completed
#   data/clean/mice_<TAG>_structural.rds   $completed
#
#   The $completed field, never $mids. The mids object still holds the
#   imputed cognitive indicator and brain age gap values, which were
#   auxiliary and were restored to NA on the way into $completed. The
#   contract is stated in the header of 03_impute.R and is verified in
#   section 2 below rather than trusted.
#
#   Indicator values remain missing on purpose and are handled by
#   full-information maximum likelihood inside each fit. Participants
#   without a gap for the modality in question are excluded from that
#   modality's model; a gap cannot be imputed.
#
# MEASUREMENT MODEL
#   Identical syntax to 01_cfa_cognitive_bifactor.R. It is repeated here
#   rather than sourced because a structural model that silently drifts
#   from the measurement model it claims to extend is the failure this
#   duplication is meant to make visible: if the two ever disagree, the
#   diff is in one place.
#
# OUTPUT  results/sem/h1_log_<TAG>.txt
#         results/sem/h1_estimates_<TAG>.csv
#         results/sem/h1_fit_<TAG>.csv
#         results/sem/h1_provenance_<TAG>.txt
# =====================================================================

library(lavaan)
library(lavaan.mi)

TAG  <- Sys.getenv("FUNCTIONAL_TAG", "pca")
ROOT <- "C:/Users/inaye/Desktop/MASTERTHESIS"
IND  <- file.path(ROOT, "data", "clean")
RES  <- file.path(ROOT, "results", "sem")
dir.create(RES, recursive = TRUE, showWarnings = FALSE)

# Bonferroni across the two modalities, as pre-specified for H1 only.
ALPHA <- 0.025
CONF  <- 1 - ALPHA          # 0.975

# Whether the two specific factors are also regressed on the covariates.
# Fixed at TRUE, on the a priori grounds given in section 5, and fixed
# before any H1 estimate was produced. The specification with the specific
# factors left exogenous is also fitted, but it is a labelled supplementary
# analysis: the primary specification is NOT revised on the basis of what
# the two of them return. Changing it after seeing the coefficients would
# be a researcher degree of freedom, which is the thing the analysis plan
# exists to remove.
REGRESS_SPECIFICS <- TRUE

# Fit two imputations first and stop on failure. The full run is 50
# datasets times two specifications times two modalities; an identification
# or syntax problem should surface in seconds rather than after minutes.
SMOKE_TEST <- TRUE

if (interactive())
  message("NOTE: run with Rscript or source(), not line by line.")

while (sink.number() > 0) sink()
try(sink(type = "message"), silent = TRUE)

con <- file(file.path(RES, sprintf("h1_log_%s.txt", TAG)), open = "wt")
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
  normalizePath(file.path(ROOT, "code", "sem", "04_h1_base.R"),
                winslash = "/", mustWork = FALSE)
})

PROV <- list(
  run_time        = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  tag             = TAG,
  alpha           = ALPHA,
  conf_level      = CONF,
  regress_specifics = REGRESS_SPECIFICS,
  r_version       = R.version.string,
  lavaan_version  = as.character(utils::packageVersion("lavaan")),
  lavaan_mi_version = as.character(utils::packageVersion("lavaan.mi")),
  script_path     = this_script,
  script_md5      = if (file.exists(this_script))
                      unname(tools::md5sum(this_script)) else NA_character_)

cat("=== provenance ===\n")
for (k in names(PROV)) cat(sprintf("  %-18s %s\n", k, PROV[[k]]))

writeLines(vapply(names(PROV),
                  function(k) sprintf("%-18s %s", k, PROV[[k]]), character(1)),
           file.path(RES, sprintf("h1_provenance_%s.txt", TAG)))

# semTools is named in the analysis plan for multiple-imputation model
# comparison. Its MI functions have since been deprecated and moved to
# lavaan.mi, which is what is used here. Same procedures, different
# package; recorded so the deviation is not silent.
cat("\nnote: multiple-imputation SEM uses lavaan.mi; semTools::runMI is deprecated\n")


# ---------------------------------------------------------------------
# 2. Load the imputations and verify the contract
# ---------------------------------------------------------------------
IND_COLS <- c("flanker", "card_sort", "list_sorting", "pattern_comp",
              "picture_seq", "ravlt", "trails_b_rev")

load_imp <- function(modality) {
  f <- file.path(IND, sprintf("mice_%s_%s.rds", TAG, modality))
  if (!file.exists(f))
    stop("imputations not found: ", f,
         "\n  Run 03_impute.R first.")
  obj <- readRDS(f)

  for (fld in c("completed", "mids", "gap_var", "modality", "restore"))
    if (is.null(obj[[fld]]))
      stop(f, ": field '", fld, "' is absent. This file was not written by ",
           "the current 03_impute.R.")
  if (obj$modality != modality)
    stop(f, ": contains modality '", obj$modality, "', expected '", modality, "'")

  # The contract: the gap and the indicators must still be missing where
  # they were missing before imputation. If a previous step filled them,
  # every estimate below would be computed from fabricated values, and
  # nothing in the output would look wrong.
  x <- obj$completed[[1]]
  if (!anyNA(x[[obj$gap_var]]))
    stop(f, ": ", obj$gap_var, " has no missing values in dataset 1. The ",
         "auxiliary imputations were not restored to NA. Do not use this file.")
  n_na_ind <- vapply(IND_COLS, function(v) sum(is.na(x[[v]])), integer(1))
  if (all(n_na_ind == 0))
    stop(f, ": no cognitive indicator has missing values in dataset 1. ",
         "This looks like $mids output rather than $completed.")

  cat(sprintf("\n%s: %d imputed datasets, %d rows each\n",
              basename(f), length(obj$completed), nrow(x)))
  cat(sprintf("  gap column %s, %d missing (correctly not imputed)\n",
              obj$gap_var, sum(is.na(x[[obj$gap_var]]))))
  obj
}

imp <- list(functional = load_imp("functional"),
            structural = load_imp("structural"))


# ---------------------------------------------------------------------
# 3. Analysis variables
# ---------------------------------------------------------------------
# One list of M data frames per modality, each holding only the columns
# the model needs, on the analysis sample for that modality.
#
# Standardisation is done within each imputed dataset. The alternative,
# using the observed-data moments for every dataset, would ignore that
# the imputed values shift the mean and spread of income and education
# slightly from one dataset to the next; standardising within keeps each
# dataset internally consistent, which is what the pooling assumes.

z <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) stop("cannot standardise a constant variable")
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

prepare <- function(obj) {
  gap <- obj$gap_var

  # The analysis sample must be identical across the M datasets, otherwise
  # the pooled estimates would refer to M different samples. sem_eligible
  # and the gap missingness are both carried through unchanged, so this is
  # a check rather than a computation.
  mask1 <- NULL
  out <- lapply(seq_along(obj$completed), function(i) {
    d <- obj$completed[[i]]
    keep <- d$sem_eligible & !is.na(d[[gap]])
    if (is.null(mask1)) mask1 <<- keep
    else if (!identical(keep, mask1))
      stop("dataset ", i, " (", obj$modality,
           "): analysis sample differs from dataset 1")
    d <- d[keep, , drop = FALSE]

    # lavaan does not accept factor predictors. Sex enters as a single
    # indicator variable with female as the reference category.
    if (!all(as.character(d$sex) %in% c("F", "M")))
      stop("unexpected value in sex: ",
           paste(setdiff(unique(as.character(d$sex)), c("F", "M")),
                 collapse = ", "))

    data.frame(
      participant_id = d$participant_id,
      pedigree       = d$pedigree,
      # Indicators standardised, as in the measurement model. This changes
      # no fit index and no standardised estimate; it puts seven variables
      # that arrive on very different scales onto one.
      setNames(as.data.frame(lapply(d[IND_COLS], z)), IND_COLS),
      BAG      = z(d[[gap]]),
      age_z    = z(d$age),
      edu_z    = z(as.numeric(d$education_num)),
      sex_male = as.integer(as.character(d$sex) == "M"),
      stringsAsFactors = FALSE)
  })

  cat(sprintf("\n%s analysis sample: n = %d\n", obj$modality, nrow(out[[1]])))
  cat("  missing per indicator (dataset 1, handled by FIML):\n")
  print(colSums(is.na(out[[1]][IND_COLS])))
  cat("  complete predictors:",
      sum(stats::complete.cases(out[[1]][c("BAG", "age_z", "edu_z", "sex_male")])),
      "of", nrow(out[[1]]), "\n")
  out
}

dat <- lapply(imp, prepare)


# ---------------------------------------------------------------------
# 4. Descriptives of the modelled variables
# ---------------------------------------------------------------------
# On the unstandardised scale, so the numbers are readable. Taken from the
# first imputed dataset and labelled as such: for income and education the
# value depends on the imputation, and reporting one dataset as if it were
# the data would misstate the precision.

cat("\n=== unstandardised descriptives, imputed dataset 1 ===\n")
for (m in names(imp)) {
  d <- imp[[m]]$completed[[1]]
  keep <- d$sem_eligible & !is.na(d[[imp[[m]]$gap_var]])
  d <- d[keep, ]
  cat(sprintf("\n%s (n = %d)\n", m, nrow(d)))
  cat(sprintf("  %-16s M = %7.2f  SD = %6.2f  [%.2f, %.2f]\n",
              imp[[m]]$gap_var, mean(d[[imp[[m]]$gap_var]]),
              sd(d[[imp[[m]]$gap_var]]),
              min(d[[imp[[m]]$gap_var]]), max(d[[imp[[m]]$gap_var]])))
  cat(sprintf("  %-16s M = %7.2f  SD = %6.2f\n", "age",
              mean(d$age), sd(d$age)))
  cat(sprintf("  %-16s M = %7.2f  SD = %6.2f\n", "education_num",
              mean(as.numeric(d$education_num)), sd(as.numeric(d$education_num))))
  cat(sprintf("  %-16s %d female, %d male\n", "sex",
              sum(d$sex == "F"), sum(d$sex == "M")))
  cat(sprintf("  gap correlation with age: r = %+.3f\n",
              cor(d[[imp[[m]]$gap_var]], d$age)))
}


# ---------------------------------------------------------------------
# 5. Model
# ---------------------------------------------------------------------
# The measurement part is the bifactor-(S-1) specification from
# 01_cfa_cognitive_bifactor.R: pattern_comp and trails_b_rev, the two
# processing-speed indicators, load on the general factor only and form
# the reference domain; a* constrains the two executive loadings to
# equality, which identifies a two-indicator specific factor; the three
# fixed zero covariances make the factors orthogonal.
#
# The structural part regresses the general factor on the gap and the
# covariates. b_bag is the H1 coefficient.
#
# WHETHER THE SPECIFIC FACTORS ARE ALSO REGRESSED
#   The analysis plan states that latent cognitive performance is regressed
#   on the gap and the covariates, and says nothing about the two specific
#   factors. That silence still has to be resolved, because leaving Exec
#   and Mem exogenous is not a neutral default. With observed exogenous
#   covariates lavaan conditions on them, and an exogenous latent factor is
#   then taken to be independent of them. The model would be asserting that
#   age, sex, education and the brain age gap are unrelated to
#   executive-specific and memory-specific variance.
#
#   That assertion is implausible on prior grounds alone. Differential
#   decline across cognitive domains in ageing is among the better
#   established findings in the literature, so age in particular is
#   expected to relate to domain-specific variance and not only to the
#   general factor. An assumption known in advance to be false, and whose
#   violation can bias the coefficient H1 is about, is not one to impose.
#
#   The primary specification therefore regresses all three factors on the
#   covariates. This leaves the b_bag equation meaning exactly what the
#   analysis plan says it means, while dropping the false zero constraints.
#
#   The specification that leaves Exec and Mem exogenous is fitted too and
#   reported as a supplementary analysis, so a reader can see what the
#   constraint costs. It is not a candidate for promotion to primary. The
#   choice was made here, from the argument above, before any H1 estimate
#   existed, and section 8 reports the comparison without acting on it.

MEASUREMENT <- '
  g    =~ flanker + card_sort + list_sorting + pattern_comp +
          picture_seq + ravlt + trails_b_rev
  Exec =~ a*flanker + a*card_sort
  Mem  =~ picture_seq + ravlt + list_sorting

  g    ~~ 0*Exec
  g    ~~ 0*Mem
  Exec ~~ 0*Mem
'

COVARIATES <- "BAG + age_z + sex_male + edu_z"

model_h1 <- function(with_specifics) {
  s <- paste0(MEASUREMENT, "\n  g ~ b_bag*BAG + age_z + sex_male + edu_z\n")
  if (with_specifics)
    s <- paste0(s,
                "  Exec ~ ", COVARIATES, "\n",
                "  Mem  ~ ", COVARIATES, "\n")
  s
}

cat("\n=== primary model syntax (REGRESS_SPECIFICS =",
    REGRESS_SPECIFICS, ") ===\n")
cat(model_h1(REGRESS_SPECIFICS))


# ---------------------------------------------------------------------
# 6. Fit
# ---------------------------------------------------------------------
# MLR with FIML for the missing indicator values, matching the estimator
# used for the measurement models. std.lv fixes the latent scale; note
# that for an endogenous factor this fixes the RESIDUAL variance to one,
# not the total variance, which is why the standardised solution rather
# than the raw loadings is what gets reported.
#
# sem.mi fits the model to each of the M datasets and pools with Rubin's
# rules. data= must be a list of data frames.

fit_one <- function(dlist, with_specifics, label) {
  cat("\n\n############ ", label, " ############\n")
  t0 <- Sys.time()
  fit <- lavaan.mi::sem.mi(
    model     = model_h1(with_specifics),
    data      = dlist,
    estimator = "MLR",
    missing   = "fiml",
    std.lv    = TRUE)
  cat(sprintf("fitted %d imputations in %.1f s\n",
              length(dlist), as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  # How many imputations produced an admissible solution. A pooled estimate
  # built from a minority of converged fits is not interpretable, so this is
  # reported before any coefficient.
  nconv <- tryCatch({
    cv <- slot(fit, "convergence")
    sum(vapply(cv, function(x) isTRUE(x$converged), logical(1)))
  }, error = function(e) NA_integer_)
  if (!is.na(nconv)) {
    cat(sprintf("converged: %d of %d imputations\n", nconv, length(dlist)))
    if (nconv < length(dlist))
      warning(label, ": only ", nconv, " of ", length(dlist),
              " imputations converged.", call. = FALSE)
  } else {
    cat("convergence count not available from this lavaan.mi version; ",
        "see the summary below\n", sep = "")
  }
  fit
}

# Smoke test before the full run. Two imputations of the structural data,
# both specifications, and the three accessors the rest of the script needs.
# A model that is not identified, a syntax slip, or a lavaan.mi signature
# that differs from the one assumed here will all fail here in seconds
# instead of after several minutes of fitting.
if (SMOKE_TEST) {
  cat("\n\n=== smoke test: 2 imputations, both specifications ===\n")
  for (ws in c(TRUE, FALSE)) {
    tf <- lavaan.mi::sem.mi(model = model_h1(ws), data = dat$structural[1:2],
                            estimator = "MLR", missing = "fiml", std.lv = TRUE)
    for (acc in c("parameterEstimates.mi", "standardizedSolution.mi")) {
      r <- get(acc, envir = asNamespace("lavaan.mi"))(tf)
      if (!nrow(r))
        stop("smoke test: ", acc, " returned no rows (specifics = ", ws, ")")
    }
    b <- lavaan.mi::standardizedSolution.mi(tf)
    b <- b[b$op == "~" & b$lhs == "g" & b$rhs == "BAG", ]
    if (nrow(b) != 1)
      stop("smoke test: g ~ BAG not found (specifics = ", ws, ")")
    cat(sprintf("  specifics regressed = %-5s  ok, b_bag = %+.4f\n", ws, b$est.std))
  }
  cat("  smoke test passed; proceeding to the full run\n")
}

# primary is the pre-specified specification. supplementary is the other
# one, fitted for reporting only; nothing downstream reads it as the H1
# result. The names are used in the output files, so that no table can be
# mistaken for the confirmatory one.
fits <- list()
for (m in c("functional", "structural")) {
  fits[[m]] <- list(
    primary       = fit_one(dat[[m]], REGRESS_SPECIFICS,
                            sprintf("H1 %s, PRIMARY (pre-specified)", m)),
    supplementary = fit_one(dat[[m]], !REGRESS_SPECIFICS,
                            sprintf("H1 %s, supplementary (not confirmatory)", m)))
}


# ---------------------------------------------------------------------
# 7. Pooled estimates
# ---------------------------------------------------------------------
# Confidence intervals at CONF = 0.975 for H1, as pre-specified. If the
# installed lavaan.mi does not accept level=, the interval is rebuilt from
# the pooled estimate and standard error rather than silently reported at
# the wrong level.

get_est <- function(fit, standardized) {
  f <- if (standardized) lavaan.mi::standardizedSolution.mi else
                         lavaan.mi::parameterEstimates.mi
  out <- tryCatch(f(fit, level = CONF),
                  error = function(e) NULL)
  if (is.null(out)) {
    out <- f(fit)
    ecol <- if (standardized) "est.std" else "est"
    if (all(c(ecol, "se") %in% names(out))) {
      k <- stats::qnorm(1 - ALPHA / 2)
      out$ci.lower <- out[[ecol]] - k * out$se
      out$ci.upper <- out[[ecol]] + k * out$se
      cat("  note: level= not supported; ", CONF * 100,
          "% interval rebuilt as est +/- ", round(k, 3),
          " * SE (normal approximation)\n", sep = "")
    }
  }
  out
}

rows <- list()
fitrows <- list()

for (m in names(fits)) {
  for (spec in c("primary", "supplementary")) {
    fit <- fits[[m]][[spec]]
    lbl <- sprintf("%s / %s", m, spec)

    cat("\n\n=== ", lbl, " ===\n")
    print(tryCatch(summary(fit, standardized = TRUE),
                   error = function(e) paste("summary() failed:", conditionMessage(e))))

    std <- get_est(fit, TRUE)
    raw <- get_est(fit, FALSE)

    reg <- std[std$op == "~", , drop = FALSE]
    cat("\n--- standardised structural coefficients, ",
        CONF * 100, "% CI ---\n", sep = "")
    print(reg, row.names = FALSE)

    b <- reg[reg$lhs == "g" & reg$rhs == "BAG", , drop = FALSE]
    if (nrow(b) != 1)
      stop(lbl, ": could not locate the g ~ BAG coefficient")

    cat(sprintf("\nH1 (%s): standardised beta = %+.4f, SE = %.4f, %s%% CI [%+.4f, %+.4f], p = %.4g\n",
                lbl, b$est.std, b$se, format(CONF * 100),
                b$ci.lower, b$ci.upper, b$pvalue))
    cat(sprintf("  decision rule: alpha = %.3f (Bonferroni across modalities)\n", ALPHA))
    cat(sprintf("  direction predicted: negative. observed: %s\n",
                if (b$est.std < 0) "negative" else "positive"))
    cat(sprintf("  p < alpha: %s\n", b$pvalue < ALPHA))

    reg$modality <- m; reg$spec <- spec; reg$type <- "standardized"
    rows[[paste(m, spec, "std")]] <- reg
    if (!is.null(raw)) {
      r2 <- raw[raw$op == "~", , drop = FALSE]
      if (nrow(r2)) {
        r2$modality <- m; r2$spec <- spec; r2$type <- "unstandardized"
        rows[[paste(m, spec, "raw")]] <- r2
      }
    }

    # Pooled fit. The chi-square for a model fitted to imputations has to be
    # pooled as well; D2 with the robust statistic is requested explicitly so
    # that the reported value is not whichever default the installed version
    # happens to use.
    # Both flavours are requested by name. Asking for "cfi" with
    # pool.robust = TRUE is not enough: when the robust computation fails
    # lavaan emits a WARNING and silently returns the naive index, and a
    # tryCatch on errors does not see it. The run then reports naive indices
    # while 01 and 02 report robust ones, with nothing in the output saying
    # which is which. Requesting cfi.robust explicitly makes the failure
    # visible as an NA instead of a wrong number wearing the right name.
    WANT <- c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr",
              "cfi.robust", "tli.robust", "rmsea.robust")

    robust_failed <- character(0)
    fm <- withCallingHandlers(
      tryCatch(
        lavaan::fitMeasures(fit, WANT, pool.method = "D2", pool.robust = TRUE),
        error = function(e)
          tryCatch(lavaan::fitMeasures(fit, WANT), error = function(e2) NULL)),
      warning = function(w) {
        m <- conditionMessage(w)
        if (grepl("robust", m, ignore.case = TRUE))
          robust_failed <<- union(robust_failed, trimws(m))
        invokeRestart("muffleWarning")
      })

    if (!is.null(fm)) {
      cat("\n--- pooled fit ---\n"); print(round(fm, 3))
      if (length(robust_failed)) {
        cat("\n  NOTE: lavaan could not compute the robust indices here:\n")
        for (msg in robust_failed) cat("    ", msg, "\n", sep = "")
        cat("  The *.robust columns are therefore NA and the plain cfi/tli/\n")
        cat("  rmsea above are the NAIVE indices, not the scaled ones. Report\n")
        cat("  them as such; 01 and 02 report robust indices and the two are\n")
        cat("  not interchangeable in a table.\n")
      }
      fitrows[[lbl]] <- data.frame(modality = m, spec = spec,
                                   index = names(fm), value = as.numeric(fm),
                                   robust_available = !length(robust_failed))
    } else {
      cat("\n--- pooled fit unavailable from this lavaan.mi version ---\n")
    }
  }
}


# ---------------------------------------------------------------------
# 8. What the exogenous-specifics constraint costs
# ---------------------------------------------------------------------
# Reported, not acted on. The primary specification was fixed in advance
# from the argument in section 5, and this comparison does not revise it
# whatever it shows. Its purpose is to tell a reader how much the zero
# constraints on the specific factors would have mattered, which is a
# question about the model, not a basis for selecting one.
#
# Interpretation for the thesis: if the two b_bag estimates agree closely,
# the constraint was harmless and the choice of specification is immaterial
# to the H1 conclusion. If they differ, the supplementary model shows what
# imposing an assumption known to be false would have done, which is the
# reason the primary specification does not impose it.

cat("\n\n=== primary vs supplementary specification, b_bag (reported only) ===\n")
cmp <- do.call(rbind, lapply(names(fits), function(m) {
  g <- lapply(c("primary", "supplementary"), function(spec) {
    e <- get_est(fits[[m]][[spec]], TRUE)
    e <- e[e$op == "~" & e$lhs == "g" & e$rhs == "BAG", ]
    data.frame(modality = m, spec = spec, beta = e$est.std, se = e$se,
               p = e$pvalue)
  })
  do.call(rbind, g)
}))
cmp$specifics_regressed <- ifelse(cmp$spec == "primary",
                                  REGRESS_SPECIFICS, !REGRESS_SPECIFICS)
print(cmp, row.names = FALSE)

for (m in unique(cmp$modality)) {
  b <- cmp$beta[cmp$modality == m]
  cat(sprintf("  %s: absolute difference in beta = %.4f\n",
              m, abs(diff(b))))
}
cat("\n  The confirmatory H1 result is the primary row. This comparison is\n")
cat("  descriptive and does not change which row that is.\n")


# ---------------------------------------------------------------------
# 9. Save
# ---------------------------------------------------------------------
# Matrix, loaded as a dependency elsewhere in the project, makes rbind an S4
# generic whose dispatch indexes list(...)[[1]], so do.call(rbind, list()) is
# an error rather than NULL. Every rbind over a collected list is guarded.
if (!length(rows)) {
  stop("no estimates were collected; nothing to write")
}

# The two accessors return different columns: standardizedSolution.mi puts
# the estimate in est.std, parameterEstimates.mi in est, and neither returns
# exactly the same set of test-statistic columns. Taking whatever each frame
# happens to have and binding them fails, because rbind requires identical
# columns. Force every frame onto one schema first: the estimate always lands
# in est, and the type column records which scale it is on.
CANON <- c("modality", "spec", "type", "lhs", "op", "rhs", "label",
           "est", "se", "z", "pvalue", "ci.lower", "ci.upper")

harmonise <- function(x) {
  if (!"est" %in% names(x) && "est.std" %in% names(x)) x$est <- x$est.std
  # lavaan.mi names the test statistic z in some tables and t in others.
  if (!"z" %in% names(x) && "t" %in% names(x)) x$z <- x$t
  for (cc in setdiff(CANON, names(x))) x[[cc]] <- NA
  x[, CANON, drop = FALSE]
}

est <- do.call(rbind, lapply(rows, harmonise))
write.csv(est, file.path(RES, sprintf("h1_estimates_%s.csv", TAG)),
          row.names = FALSE)
cat("\nwrote", file.path(RES, sprintf("h1_estimates_%s.csv", TAG)), "\n")

if (length(fitrows)) {
  write.csv(do.call(rbind, fitrows),
            file.path(RES, sprintf("h1_fit_%s.csv", TAG)), row.names = FALSE)
  cat("wrote", file.path(RES, sprintf("h1_fit_%s.csv", TAG)), "\n")
}

cat("\nH1 complete. The confirmatory coefficients are the g ~ BAG rows of the\n")
cat("primary specification, at alpha =", ALPHA, "with", CONF * 100, "% intervals.\n")

close_log()
