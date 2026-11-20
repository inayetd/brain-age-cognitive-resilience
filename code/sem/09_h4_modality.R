# =====================================================================
# H4 - does cognitive resilience depend on which modality represents
#      brain status?  PRE-SPECIFIED EXPLORATORY
#
# WHAT THE PREREGISTRATION ASKS FOR
#   (i)   The Pearson correlation between functional and structural BAG
#         with a family-clustered bootstrap interval, Lin's concordance
#         correlation coefficient, and the correlation between the two
#         descriptive resilience scores.
#   (ii)  Latent cognitive performance regressed simultaneously on both
#         standardised gaps plus age, sex and education. Both coefficients,
#         the unique increment in explained variance from dropping each,
#         and -- in the words of the registered plan -- "a one-degree-of-
#         freedom Wald test of the null hypothesis that the two
#         coefficients are equal".
#
#         The Wald test is therefore the registered one. A scaled
#         likelihood-ratio test of the same hypothesis is reported beside
#         it as a supplementary check, because an earlier draft of the plan
#         specified that instead and the two should agree. They test the
#         same null but need not give identical finite-sample p values.
#   (iii) Six predictor-level comparisons of resilience coefficients across
#         modalities, by paired bootstrap.
#
#   This script does (i) and (ii). Part (iii) is a separate bootstrap.
#
# WHY THIS RUNS ON OBSERVED DATA AND NEEDS NO IMPUTATION
#   The H4 model contains both gaps, age, sex, education and the seven
#   cognitive indicators. In the intersection sample of 1,342:
#
#     both gaps              0 missing   -- H4 is defined on this sample
#     age, sex               0 missing
#     education             24 missing   1.8%
#     cognitive indicators   missing, but handled by FIML inside the
#                            measurement model, which is how every other
#                            script in this pipeline treats them. They are
#                            never imputed for analysis anywhere: 03 writes
#                            them back to NA on purpose.
#
#   So the only quantity an imputation could recover is education for 24
#   people. Full-information maximum likelihood is the appropriate
#   treatment for missingness on that scale in a single covariate, and it
#   is the treatment this pipeline already applies to every indicator
#   inside a measurement model. That is the reason, not the run time.
#
#   fixed.x = FALSE is what makes it work: missing = "fiml" alone leaves
#   exogenous predictors fixed, and lavaan then deletes any case with a
#   missing value among them. Setting it FALSE brings age, sex, education
#   and both gaps into the likelihood, so the 24 cases are retained.
#
#   The cost of that is an assumption worth stating. fixed.x = FALSE models
#   the joint distribution of the predictors, including their means,
#   variances and covariances, and education is a six-level ordinal score
#   while sex is binary. Neither is normal. With 1,342 observations and
#   missingness confined to 24 education values this is a mild assumption,
#   and MLR gives robust standard errors, but robustness to non-normality
#   is not the same thing as validity of the missing-data mechanism. The
#   fit on participants with observed education is reported beside it. That
#   comparison is reassuring; it cannot establish that the missingness is
#   ignorable.
#
# INFERENCE
#   H4 is exploratory. No support or rejection verdict is assigned to it,
#   the correlation and the equality test are each reported once without
#   multiplicity correction, and nothing here enters the confirmatory H1,
#   H2 or H3 families.
#
# WHAT LIMITS THE INTERPRETATION, FROM THE PLAN ITSELF
#   The plan states that a modality difference can only arise if the two
#   gaps carry partly non-overlapping information AND each explains
#   non-negligible variance in latent cognition. The first condition is
#   addressed in section 2, which reports the correlation between the two
#   gaps. The second is addressed in section 4, which computes each gap's
#   unique increment in the R-squared of latent g from the full and
#   reduced models fitted here. Read both before reading the test: if
#   neither gap explains much, there is little for a difference to be a
#   difference in.
#
#   Those increments are the quantity to quote. They are not the same as
#   a squared standardised coefficient from a single-gap model, which is
#   a different model with a different adjustment set.
#
#   A NOTE ON PRECISION, WHICH THE SCRIPT COMPUTES RATHER THAN ASSUMES
#   Entering correlated predictors together affects the precision of their
#   coefficients, and the precision of their difference follows from
#
#       Var(b_s - b_f) = Var(b_s) + Var(b_f) - 2 Cov(b_s, b_f)
#
#   so it turns entirely on the covariance between the two estimates. In
#   ordinary regression, positively correlated predictors usually give
#   negatively correlated coefficient estimates, which would make the
#   difference the less precisely estimated quantity -- but that is a
#   two-predictor result, and this model also carries covariates, a latent
#   measurement model, missing data and robust estimation. No amount of
#   inflation is assumed here.
#
#   Section 5 reads Var(b_s), Var(b_f) and Cov(b_s, b_f) off the fitted
#   model and prints the decomposition with the ratio of standard errors.
#   Both the coefficients and their contrast are informative: each
#   coefficient estimates one modality's association with cognition
#   conditional on the other modality and the covariates, and the Wald
#   test asks directly whether the two differ.
#
# HOW TO RUN
#   Rscript code/sem/09_h4_modality.R
#
# Do not submit this file to the console line by line. Every stop() below
# is a guard meant to halt the run.
#
# INPUT   data/clean/analysis_table_<TABLE_TAG>.csv
#         results/bag_sensitivity.csv
#         results/sem/h3_resilience_scores_<TAG>.csv   (optional)
# OUTPUT  results/sem/h4_estimates_<TAG>.csv
#         results/sem/h4_agreement_<TAG>.csv
#         results/sem/h4_log_<TAG>.txt
#         results/sem/h4_provenance_<TAG>.txt
# =====================================================================

suppressPackageStartupMessages(library(lavaan))

# Project root. The absolute path matches every other script in code/sem,
# so the script works from a session started anywhere -- including R's
# default working directory, which is where it fails otherwise. Set the
# PROJECT_ROOT environment variable to override it if the folder moves.
ROOT <- Sys.getenv("PROJECT_ROOT", "C:/Users/inaye/Desktop/MASTERTHESIS")
if (!dir.exists(file.path(ROOT, "data", "clean"))) {
  # Fall back to walking up from the working directory, for a machine on
  # which the path above does not exist.
  here <- normalizePath(".", mustWork = TRUE)
  for (up in c(".", "..", "../..", "../../..")) {
    cand <- normalizePath(file.path(here, up), mustWork = FALSE)
    if (dir.exists(file.path(cand, "data", "clean")) &&
        dir.exists(file.path(cand, "results"))) { ROOT <- cand; break }
  }
}
if (!dir.exists(file.path(ROOT, "data", "clean")))
  stop("cannot find the project. Edit ROOT at the top of this script, or\n",
       "  Sys.setenv(PROJECT_ROOT = \"C:/Users/inaye/Desktop/MASTERTHESIS\")\n",
       "Tried: ", ROOT, ", and upward from ", getwd())

TAG       <- Sys.getenv("FUNCTIONAL_TAG", "pca_it50")
TABLE_TAG <- Sys.getenv("TABLE_TAG", sub("_it[0-9]+$", "", TAG))
IND  <- file.path(ROOT, "data", "clean")
RES  <- file.path(ROOT, "results", "sem")
N_BOOT <- 2000
REGRESS_SPECIFICS <- TRUE   # as in 04, 06 and 06b

LOG <- file.path(RES, sprintf("h4_log_%s.txt", TAG))

# Release any sink left open by a previous run before opening a new one.
# Without this, re-running the script after an error stacks sinks: each
# run adds two more, and the console stays silent until as many sink()
# calls have been made as there were opens. That is the usual cause of
# "the output never came back".
release <- function() {
  # sink.number(type = "message") returns the CONNECTION NUMBER in use
  # for messages, not a count of diversions: it is 2 (stderr) when no
  # sink is in place, and message sinks have no stack. A while loop on
  # "> 0" therefore never terminates. Restore once, conditionally.
  if (sink.number(type = "message") != 2L) sink(type = "message")
  while (sink.number(type = "output") > 0L) sink()
  # The connection has to be closed as well, not just unsunk. On Windows
  # an open handle keeps the log file locked, so the next run cannot
  # overwrite it.
  if (exists("con", envir = globalenv())) {
    cc <- get("con", envir = globalenv())
    if (inherits(cc, "connection")) {
      ok <- tryCatch({ isOpen(cc) }, error = function(e) FALSE)
      if (isTRUE(ok)) try(close(cc), silent = TRUE)
    }
  }
  invisible(NULL)
}
release()

con <- file(LOG, open = "wt")
# split = TRUE on the output sink, so results appear on the console as
# well as in the log. The message sink has no split: R does not offer one,
# so error and warning text goes to the file only. That is why a failure
# can look like silence. The handler below therefore releases the sinks
# first and then reprints the actual error on the console, rather than
# only pointing at the log.
sink(con, split = TRUE); sink(con, type = "message")
options(error = function() {
  msg <- tryCatch(geterrmessage(), error = function(e) "")
  release()
  cat("\n*** THE SCRIPT STOPPED ***\n", msg,
      "\nFull output up to that point: ", LOG, "\n", sep = "")
})

PROV <- list(
  run_time   = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  analysis   = "H4, pre-specified exploratory",
  tag        = TAG, n_boot = N_BOOT,
  missing_data = "FIML on observed data; no imputation required",
  regress_specifics = REGRESS_SPECIFICS,
  r_version  = R.version.string,
  lavaan_version = as.character(utils::packageVersion("lavaan")))
cat("=== H4: modality dependence of cognitive resilience ===\n")
for (k in names(PROV)) cat(sprintf("  %-18s %s\n", k, PROV[[k]]))
writeLines(vapply(names(PROV), function(k) sprintf("%-18s %s", k, PROV[[k]]),
                  character(1)),
           file.path(RES, sprintf("h4_provenance_%s.txt", TAG)))


# ---------------------------------------------------------------------
# 1. Data
# ---------------------------------------------------------------------
IN <- file.path(IND, sprintf("analysis_table_%s.csv", TABLE_TAG))
if (!file.exists(IN)) stop("not found: ", IN)
# An empty cell arrives from read.csv as "" rather than NA and survives
# every is.na() test, so na.strings has to name it.
d <- read.csv(IN, stringsAsFactors = FALSE, na.strings = c("", "NA"))

# pandas writes booleans as "True"/"False"; d[d$sem_eligible, ] with a
# character vector performs NAME-based row indexing and silently returns a
# frame of all-NA rows.
if (is.character(d$sem_eligible)) d$sem_eligible <- d$sem_eligible == "True"
if (!is.logical(d$sem_eligible))
  stop("sem_eligible is ", class(d$sem_eligible), ", expected logical")

EDU_LEVELS <- c(
  "Primary Level Education",
  "Some Secondary Level Education",
  "Secondary Level Education or Equivalent",
  "Some College / Associates Degree / Trade or Vocational Training",
  "Four Years College or Undergraduate Degree",
  "Graduate Degree")
unknown <- setdiff(unique(d$education[!is.na(d$education)]), EDU_LEVELS)
if (length(unknown))
  stop("unrecognised education level(s): ", paste(unknown, collapse = " | "))
d$education_num <- as.integer(factor(d$education, levels = EDU_LEVELS,
                                     ordered = TRUE))
if (!all(stats::na.omit(unique(as.character(d$sex))) %in% c("F", "M")))
  stop("unexpected value in sex")

for (v in c("BAG_structural", "BAG_functional"))
  if (!v %in% names(d)) stop("absent from the analysis table: ", v)

IND_COLS <- c("flanker", "card_sort", "list_sorting", "pattern_comp",
              "picture_seq", "ravlt", "trails_b_rev")
z <- function(v) {
  s <- stats::sd(v, na.rm = TRUE)
  if (!is.finite(s) || s == 0) stop("cannot standardise a constant variable")
  as.numeric((v - mean(v, na.rm = TRUE)) / s)
}

# The intersection sample: cognitive eligibility and BOTH gaps.
keep <- d$sem_eligible & !is.na(d$BAG_structural) & !is.na(d$BAG_functional)
x <- d[keep, , drop = FALSE]
dat <- data.frame(
  BAG_s    = z(x$BAG_structural),
  BAG_f    = z(x$BAG_functional),
  age_z    = z(x$age),
  edu_z    = z(as.numeric(x$education_num)),
  sex_male = as.integer(as.character(x$sex) == "M"))
for (v in IND_COLS) dat[[v]] <- z(x[[v]])
N <- nrow(dat)

cat(sprintf("\n\nH4 intersection sample: n = %d\n", N))
cat("  missing in the H4 model variables:\n")
for (v in c("BAG_s", "BAG_f", "age_z", "sex_male", "edu_z"))
  cat(sprintf("    %-9s %3d\n", v, sum(is.na(dat[[v]]))))
cat("  cognitive indicators, handled by FIML:\n")
for (v in IND_COLS) cat(sprintf("    %-13s %3d\n", v, sum(is.na(dat[[v]]))))
cat(sprintf("\n  r(BAG_s, BAG_f) = %.3f\n",
            stats::cor(dat$BAG_s, dat$BAG_f)))


# ---------------------------------------------------------------------
# 2. (i) Agreement between the two gaps
# ---------------------------------------------------------------------
cat("\n\n=== (i) agreement between the two gaps ===\n")
bs <- read.csv(file.path(ROOT, "results", "bag_sensitivity.csv"),
               stringsAsFactors = FALSE)
bs <- merge(bs, d[, c("participant_id", "pedigree")], by = "participant_id",
            all.x = TRUE)
SC <- "BAG_structural_corrected"; FC <- "BAG_functional_pca_corrected"
for (v in c(SC, FC)) if (!v %in% names(bs)) stop("absent: ", v)
b <- bs[!is.na(bs[[SC]]) & !is.na(bs[[FC]]), ]
a <- b[[SC]]; f <- b[[FC]]

r   <- stats::cor(a, f)
rho <- stats::cor(a, f, method = "spearman")
# Lin's concordance correlation: agreement in value, not only in rank.
ccc <- 2 * stats::cov(a, f) /
  (stats::var(a) + stats::var(f) + (mean(a) - mean(f))^2)

# Family-clustered bootstrap rather than a Fisher z interval, which assumes
# independent observations; resampling families keeps relatives together,
# matching the fold structure of the brain age models.
fam <- as.character(b$pedigree)
if (anyNA(fam)) {
  warning("some pedigree values are missing; each becomes its own cluster")
  fam[is.na(fam)] <- paste0("solo_", which(is.na(fam)))
}
fams <- unique(fam)
idx  <- split(seq_along(fam), fam)
set.seed(2026)
boot <- vapply(seq_len(N_BOOT), function(i) {
  ii <- unlist(idx[sample(fams, length(fams), replace = TRUE)],
               use.names = FALSE)
  stats::cor(a[ii], f[ii])
}, numeric(1))
ci <- stats::quantile(boot, c(.025, .975), names = FALSE)

cat(sprintf("  n = %d in %d families\n", nrow(b), length(fams)))
cat(sprintf("  Pearson r       %.3f   95%% CI [%.3f, %.3f]\n", r, ci[1], ci[2]))
cat(sprintf("  Spearman rho    %.3f\n", rho))
cat(sprintf("  Lin's CCC       %.3f\n", ccc))
cat(sprintf("  shared variance %.1f%%\n", 100 * r^2))

# The correlation between the two descriptive resilience scores.
RS <- file.path(RES, sprintf("h3_resilience_scores_%s.csv", TAG))
res_r <- NA_real_
if (file.exists(RS)) {
  s <- read.csv(RS, stringsAsFactors = FALSE)
  # Merged rather than reshaped. reshape() names the new columns from the
  # values it finds in `modality`, so the result depends on what is in the
  # file; a merge names them here and fails loudly if a modality is absent.
  ss <- s[s$modality == "structural", c("participant_id", "resilience")]
  ff <- s[s$modality == "functional", c("participant_id", "resilience")]
  if (!nrow(ss) || !nrow(ff))
    stop("h3_resilience_scores has no rows for one of the modalities; ",
         "found: ", paste(unique(s$modality), collapse = ", "))
  names(ss)[2] <- "r_structural"; names(ff)[2] <- "r_functional"
  w <- merge(ss, ff, by = "participant_id")
  w <- w[stats::complete.cases(w), ]
  res_r <- stats::cor(w$r_structural, w$r_functional)
  cat(sprintf("\n  the two descriptive resilience scores correlate %.3f (n = %d)\n",
              res_r, nrow(w)))
  cat("  This is a descriptive comparison of the two resilience scores,\n")
  cat("  not a substitute for the coefficient comparisons in part (iii).\n")
  cat("  Both scores share the same cognitive outcome and covariates and\n")
  cat("  differ only in which gap was residualised out, so a high value\n")
  cat("  here indicates that the two operationalisations are close to one\n")
  cat("  another -- which is context for part (iii), not a result from it.\n")
} else {
  cat("\n  note: ", basename(RS), " not found; score correlation skipped\n",
      sep = "")
}

write.csv(data.frame(
  index = c("pearson_r", "ci_lower", "ci_upper", "spearman_rho", "lins_ccc",
            "shared_variance", "n", "n_families", "resilience_score_r"),
  value = c(r, ci[1], ci[2], rho, ccc, r^2, nrow(b), length(fams), res_r)),
  file.path(RES, sprintf("h4_agreement_%s.csv", TAG)), row.names = FALSE)


# ---------------------------------------------------------------------
# 3. Model syntax
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
COV <- "age_z + sex_male + edu_z"
# Braces are required. Without them the function body is the if-expression
# alone, the parser treats it as complete at the end of that line, and the
# else on the next line is orphaned: "unexpected 'else'".
spec <- function(terms) {
  if (REGRESS_SPECIFICS) {
    paste0("  Exec ~ ", terms, " + ", COV, "\n",
           "  Mem  ~ ", terms, " + ", COV, "\n")
  } else {
    ""
  }
}

# b_s and b_f are labelled so that their difference is a defined parameter
# with a standard error of its own.
MOD_FREE <- paste0(MEASUREMENT,
  "\n  g ~ b_s*BAG_s + b_f*BAG_f + ", COV, "\n",
  spec("BAG_s + BAG_f"), "  diff := b_s - b_f\n")
# The same model with one equality imposed: exactly one constraint, so the
# nested comparison has one degree of freedom.
MOD_EQ <- paste0(MEASUREMENT,
  "\n  g ~ b*BAG_s + b*BAG_f + ", COV, "\n", spec("BAG_s + BAG_f"))

# The reduced models for the increments in explained variance keep BOTH
# gaps as observed variables and fix only the relevant path to g at zero.
# Deleting a gap from the model entirely would change the set of observed
# variables, and under FIML that changes the information available for the
# missing education values and the missing cognitive indicators -- so the
# two models would no longer differ in one thing only.
MOD_S0 <- paste0(MEASUREMENT,
  "\n  g ~ 0*BAG_s + BAG_f + ", COV, "\n", spec("BAG_s + BAG_f"))
MOD_F0 <- paste0(MEASUREMENT,
  "\n  g ~ BAG_s + 0*BAG_f + ", COV, "\n", spec("BAG_s + BAG_f"))

# fixed.x = FALSE so that FIML covers the 24 missing education values
# rather than deleting those rows. missing = "fiml" alone does not: lavaan
# treats exogenous predictors as fixed by default and drops cases with any
# missing among them.
fit_one <- function(model, label) {
  cat("\n---- ", label, " ----\n")
  warns <- character(0)
  f <- withCallingHandlers(
    lavaan::sem(model, data = dat, estimator = "MLR", missing = "fiml",
                std.lv = TRUE, fixed.x = FALSE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (!lavaan::lavInspect(f, "converged")) stop("did not converge: ", label)
  cat(sprintf("  n used = %d\n", lavaan::lavInspect(f, "nobs")))

  # Convergence is not admissibility. A model can converge on a solution
  # with a negative variance or a non-positive-definite latent covariance
  # matrix, and lavaan reports that as a warning rather than an error, so
  # it is easy to read past. Both are checked here and reported, not
  # silently tolerated.
  bad <- character(0)
  th <- tryCatch(diag(lavaan::lavInspect(f, "theta")),
                 error = function(e) numeric(0))
  if (length(th) && any(th < 0))
    bad <- c(bad, paste0("negative residual variance: ",
                         paste(names(th)[th < 0], collapse = ", ")))
  pv <- tryCatch(diag(lavaan::lavInspect(f, "psi")),
                 error = function(e) numeric(0))
  if (length(pv) && any(pv < 0))
    bad <- c(bad, paste0("negative latent variance: ",
                         paste(names(pv)[pv < 0], collapse = ", ")))
  cl <- tryCatch(lavaan::lavInspect(f, "cov.lv"), error = function(e) NULL)
  if (!is.null(cl) && nrow(cl) > 1) {
    ev <- min(eigen(cl, symmetric = TRUE, only.values = TRUE)$values)
    if (ev <= 0)
      bad <- c(bad, sprintf(
        "latent covariance matrix not positive definite (smallest eigenvalue %.2e)",
        ev))
  }
  if (length(warns))
    for (w in unique(warns)) cat("  WARNING from lavaan: ", w, "\n", sep = "")
  if (length(bad)) {
    for (b in bad) cat("  INADMISSIBLE: ", b, "\n", sep = "")
    cat("  The estimates below should not be reported until this is\n")
    cat("  resolved. Convergence alone does not make a solution valid.\n")
  } else {
    cat("  admissible: no negative variances, latent covariance PD\n")
  }
  f
}

fit_free <- fit_one(MOD_FREE, "free")
fit_eq   <- fit_one(MOD_EQ,   "equality constrained")
fit_s0   <- fit_one(MOD_S0,   "structural path to g fixed at zero")
fit_f0   <- fit_one(MOD_F0,   "functional path to g fixed at zero")

FIT <- c("chisq.scaled", "df", "cfi.robust", "tli.robust", "rmsea.robust",
         "srmr")
cat("\n--- fit of the free model ---\n")
print(round(lavaan::fitMeasures(fit_free, FIT), 4))


# ---------------------------------------------------------------------
# 4. (ii) Coefficients, unique increments, equality test
# ---------------------------------------------------------------------
cat("\n\n=== (ii) both gaps in one model ===\n")
std <- lavaan::standardizedSolution(fit_free, ci = TRUE, level = .95)
rows <- list()
for (tm in c("BAG_s", "BAG_f", "age_z", "sex_male", "edu_z")) {
  r0 <- std[std$op == "~" & std$lhs == "g" & std$rhs == tm, , drop = FALSE]
  if (nrow(r0) != 1) stop("expected one row for g ~ ", tm)
  rows[[length(rows) + 1]] <- data.frame(
    quantity = "coefficient", term = tm, n = N,
    est = r0$est.std, se = r0$se, z = r0$z, df = NA_real_,
    pvalue = r0$pvalue,
    ci.lower = r0$ci.lower, ci.upper = r0$ci.upper,
    stringsAsFactors = FALSE)
}
cat("\n  standardised coefficients on latent g\n")
print(do.call(rbind, rows)[, c("term", "est", "se", "pvalue",
                               "ci.lower", "ci.upper")],
      row.names = FALSE, digits = 4)

r2 <- function(f) unname(lavaan::lavInspect(f, "rsquare")["g"])
R2_full <- r2(fit_free)
inc_s <- R2_full - r2(fit_s0)   # structural path to g set to zero
inc_f <- R2_full - r2(fit_f0)   # functional path to g set to zero
cat("\n  unique increment in the R-squared of latent g\n")
cat(sprintf("    full model                  R2 = %.4f\n", R2_full))
cat(sprintf("    structural gap, unique      %.4f\n", inc_s))
cat(sprintf("    functional gap, unique      %.4f\n", inc_f))
cat("    Each reduced model keeps both gaps as observed variables and\n")
cat("    fixes only that gap's path to g at zero, so the two models differ\n")
cat("    in one parameter and in nothing else. These are model-implied\n")
cat("    increments in a latent R-squared and should not be read as if\n")
cat("    they were increments in an ordinary regression on an observed\n")
cat("    outcome.\n")
for (nm in c("BAG_s", "BAG_f")) {
  rows[[length(rows) + 1]] <- data.frame(
    quantity = "unique_r2", term = nm, n = N,
    est = if (nm == "BAG_s") inc_s else inc_f, se = NA_real_,
    z = NA_real_, df = NA_real_, pvalue = NA_real_,
    ci.lower = NA_real_, ci.upper = NA_real_, stringsAsFactors = FALSE)
}

cat("\n\n=== the equality test ===\n")

# THE REGISTERED TEST: a one-degree-of-freedom Wald test of b_s = b_f,
# taken from the defined parameter.
#
# SCALE. The constraint tested is b_s == b_f, on the model's own
# coefficients for the two externally z-scored BAG predictors. That is not
# the same object as the fully standardised coefficients printed above: a
# standardised coefficient is also divided by its own predictor's
# model-implied SD, and under fixed.x = FALSE those SDs are estimated
# rather than fixed at one, so they need not be exactly equal to each
# other. The two contrasts will be close here because both predictors were
# z-scored before fitting, but do not describe the number below as the
# difference between the two standardised betas, and do not claim the test
# is numerically identical to testing their equality.
#
# std.lv = TRUE fixes the SCALE of each latent variable, which for an
# endogenous latent means its residual variance, not its total variance.
# Latent g does not have unit total variance here.
# The registered form: a one-degree-of-freedom Wald chi-square on the
# constraint itself. lavTestWald returns exactly that, so the output
# matches the wording of the plan rather than an equivalent restatement.
wald <- tryCatch(lavaan::lavTestWald(fit_free, constraints = "b_s == b_f"),
                 error = function(e) conditionMessage(e))
cat("\n  Wald test of b_s == b_f (registered, one degree of freedom)\n")
if (is.character(wald)) {
  cat("    unavailable: ", wald, "\n", sep = "")
} else {
  # Field names are read defensively: lavTestWald returns stat, df and
  # p.value in this version, but a renamed element would otherwise be
  # written to the results file as NA without comment.
  pick <- function(x, nm) if (nm %in% names(x)) x[[nm]] else NA_real_
  wstat <- pick(wald, "stat"); wdf <- pick(wald, "df")
  wp    <- pick(wald, "p.value")
  if (anyNA(c(wstat, wdf, wp)))
    stop("lavTestWald returned unexpected fields: ",
         paste(names(wald), collapse = ", "))
  cat(sprintf("    chi-square = %.4f   df = %g   p = %.4g\n",
              wstat, wdf, wp))
  if (!is.na(wdf) && wdf != 1)
    cat(sprintf("    WARNING: df = %g, expected 1\n", wdf))
  rows[[length(rows) + 1]] <- data.frame(
    quantity = "equality_wald_registered", term = "b_s == b_f", n = N,
    est = wstat, se = NA_real_, z = NA_real_, df = wdf,
    pvalue = wp, ci.lower = NA_real_, ci.upper = NA_real_,
    stringsAsFactors = FALSE)
}

# The same hypothesis as an estimate: the size of the difference, its
# standard error and its interval, which the chi-square above does not
# give. z here squares to the chi-square above.
pe <- lavaan::parameterEstimates(fit_free, ci = TRUE, level = .95)
dd <- pe[pe$op == ":=" & pe$lhs == "diff", , drop = FALSE]
if (nrow(dd) == 1) {
  cat(sprintf("    b_s - b_f = %+.4f   SE = %.4f   z = %+.3f   p = %.4g\n",
              dd$est, dd$se, dd$z, dd$pvalue))
  cat(sprintf("    95%% CI [%+.4f, %+.4f]\n", dd$ci.lower, dd$ci.upper))
  rows[[length(rows) + 1]] <- data.frame(
    quantity = "equality_difference", term = "b_s - b_f", n = N,
    est = dd$est, se = dd$se, z = dd$z, df = NA_real_,
    pvalue = dd$pvalue,
    ci.lower = dd$ci.lower, ci.upper = dd$ci.upper,
    stringsAsFactors = FALSE)
}

# Where the precision of the contrast comes from, computed rather than
# asserted:
#
#     Var(b_s - b_f) = Var(b_s) + Var(b_f) - 2 Cov(b_s, b_f)
#
# so whether the difference is more or less precisely estimated than
# either coefficient depends entirely on that covariance, which this model
# has to supply. Nothing is assumed about its sign or size.
#
# The indices come from the parameter table rather than from matching
# rownames: lavaan names a LABELLED parameter by its label, so the entries
# here are "b_s" and "b_f" and not "g~BAG_s". Matching on the latter finds
# nothing and would skip this block without saying so.
V  <- lavaan::vcov(fit_free)
pt <- lavaan::parTable(fit_free)
free_ix <- function(lbl) {
  i <- unique(pt$free[pt$label == lbl & pt$free > 0])
  if (length(i) != 1L)
    stop("expected exactly one free parameter labelled '", lbl,
         "'; found ", length(i))
  i
}
i_s <- free_ix("b_s")
i_f <- free_ix("b_f")
vs <- V[i_s, i_s]; vf <- V[i_f, i_f]; cv <- V[i_s, i_f]
se_d <- sqrt(vs + vf - 2 * cv)
cat("\n  precision of the contrast, from the fitted model\n")
cat(sprintf("    SE(b_s) = %.4f   SE(b_f) = %.4f\n", sqrt(vs), sqrt(vf)))
cat(sprintf("    Cov(b_s, b_f) = %+.6f   correlation of the estimates = %+.3f\n",
            cv, cv / sqrt(vs * vf)))
cat(sprintf("    SE(b_s - b_f) = %.4f\n", se_d))
cat(sprintf("    ratio to SE(b_s) = %.2f, to SE(b_f) = %.2f\n",
            se_d / sqrt(vs), se_d / sqrt(vf)))
cat("    A ratio above one means the contrast is the less precisely\n")
cat("    estimated quantity here; below one, the more. Both the contrast\n")
cat("    and the two coefficients are informative -- the coefficients\n")
cat("    estimate each modality's association conditional on the other,\n")
cat("    and the Wald test asks directly whether the two differ.\n")

# Supplementary: the scaled likelihood-ratio test of the same null. An
# earlier draft of the plan specified this instead, and the two should
# agree. Checked for nesting and for a one-degree-of-freedom difference
# before being reported, because a silent mismatch would be invisible.
lrt <- tryCatch(lavaan::lavTestLRT(fit_eq, fit_free,
                                   method = "satorra.bentler.2010"),
                error = function(e) conditionMessage(e))
cat("\n  scaled likelihood-ratio test (supplementary, same null)\n")
if (is.character(lrt)) {
  cat("    unavailable: ", lrt, "\n", sep = "")
} else {
  print(lrt)
  k <- nrow(lrt)
  dfd <- lrt[k, "Df diff"]
  if (!is.na(dfd) && dfd != 1)
    cat(sprintf("    WARNING: the models differ by %g degrees of freedom, ",
                dfd),
        "expected 1. Check that only the equality constraint differs.\n",
        sep = "")
  rows[[length(rows) + 1]] <- data.frame(
    quantity = "equality_scaled_LRT_supplementary", term = "b_s = b_f",
    n = N, est = lrt[k, "Chisq diff"], se = NA_real_, z = NA_real_,
    df = dfd, pvalue = lrt[k, "Pr(>Chisq)"],
    ci.lower = NA_real_, ci.upper = NA_real_, stringsAsFactors = FALSE)
}

# Participants with observed education, as a check that handling the 24
# missing values inside the likelihood changed nothing. This is "complete
# on education", not complete-case: the cognitive indicators still carry
# their missingness and are still handled by FIML in both fits.
cat("\n\n--- check: complete on education ---\n")
dat_cc <- dat[!is.na(dat$edu_z), , drop = FALSE]
fit_cc <- lavaan::sem(MOD_FREE, data = dat_cc, estimator = "MLR",
                      missing = "fiml", std.lv = TRUE, fixed.x = FALSE)
if (lavaan::lavInspect(fit_cc, "converged")) {
  s2 <- lavaan::standardizedSolution(fit_cc)
  for (tm in c("BAG_s", "BAG_f")) {
    v1 <- std$est.std[std$op == "~" & std$lhs == "g" & std$rhs == tm]
    v2 <- s2$est.std[s2$op == "~" & s2$lhs == "g" & s2$rhs == tm]
    cat(sprintf("  %-6s  all cases %+0.4f (n = %d)   education observed %+0.4f (n = %d)\n",
                tm, v1, N, v2, nrow(dat_cc)))
  }
} else cat("  the complete-on-education model did not converge\n")

write.csv(do.call(rbind, rows),
          file.path(RES, sprintf("h4_estimates_%s.csv", TAG)),
          row.names = FALSE)


# ---------------------------------------------------------------------
# 5. How to read this
# ---------------------------------------------------------------------
cat("\n\nHOW TO READ IT\n")
cat("  The Wald test answers the specific question of whether the two\n")
cat("  coefficients are equal. It is the direct answer to that question,\n")
cat("  not a verdict on H4 as a whole: H4 is a three-part exploratory\n")
cat("  hypothesis and this script covers two of them.\n")
cat("\n  The individual coefficients are informative too. Each estimates\n")
cat("  one modality's association with latent cognition conditional on\n")
cat("  the other modality and the covariates. How their precision\n")
cat("  compares with that of the contrast is in the decomposition printed\n")
cat("  above, computed from this fit rather than assumed.\n")
cat("\n  H4 is exploratory. No support or rejection verdict is assigned, no\n")
cat("  multiplicity correction is applied to the correlation or to the\n")
cat("  equality test, and nothing here enters the H1, H2 or H3 families.\n")
cat("\n  Part (iii), the six predictor-level differences by paired\n")
cat("  bootstrap, is not in this script.\n")

cat("\nwrote\n")
cat("  ", file.path(RES, sprintf("h4_estimates_%s.csv", TAG)), "\n")
cat("  ", file.path(RES, sprintf("h4_agreement_%s.csv", TAG)), "\n")
cat("  ", LOG, "\n")

# Clear the error handler before releasing, so that a later error in the
# session cannot call release() a second time against a closed connection.
options(error = NULL)
release()
