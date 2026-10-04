# ====================================================================
# Diagnostic for the imputation. Does not impute or model anything.
# Reports whether the chained-equations chains mixed. Kept separate so a plotting
# failure cannot cost a twenty-minute imputation run.
# ====================================================================

# =====================================================================
# Convergence diagnostics for the multiple imputation
#
# Reads the objects written by 03_impute.R and reports whether the
# chained-equations chains mixed. Nothing is imputed or modelled here.
#
# HOW TO RUN
#   Rscript code/sem/03b_convergence.R
#
# WHY THIS IS A SEPARATE SCRIPT
#   The imputation takes around twenty minutes; these diagnostics take
#   seconds. Keeping them apart means a plotting failure cannot cost the
#   imputation, and that the diagnostics can be regenerated without
#   repeating the run. An earlier version had them in the same script and
#   a single plotting error discarded the whole computation.
#
# WHAT IS DIAGNOSED
#   Only the five variables whose imputed values are kept. The cognitive
#   indicators, the psychosocial indicators and the brain age gap are
#   imputed as auxiliary variables and then restored to missing, so their
#   mixing cannot affect any reported estimate.
#
#   Two diagnostics, both computed from imp$chainMean and imp$chainVar,
#   the per-iteration summaries mice records for each chain.
#
#   Trace plots: chain mean and chain SD against iteration, one line per
#   imputation. Converged chains look like overlapping noise with no
#   trend. A chain still drifting at the last iteration means N_ITER in
#   03_impute.R is too low.
#
#   WHICH OF THE TWO IS THE PRE-SPECIFIED CRITERION
#     The analysis plan specifies 50 imputed datasets, 20 iterations per
#     chain, and convergence assessed by trace plots of the chain means
#     and standard deviations. The trace plots are therefore the criterion
#     of record. The statistic below is a supplementary diagnostic added
#     afterwards; it is reported for completeness and does not replace the
#     pre-specified check.
#
#     It also has to be read with its limits in mind. The Gelman-Rubin
#     statistic was designed for a few long chains. Here there are 50
#     chains of 20 iterations, so the second half gives ten points per
#     chain. In that regime the statistic is noisy and biased upward, and
#     values slightly above 1.05 are expected under perfect mixing. A
#     value of that size is a reason to look at the corresponding trace
#     panel, not a verdict on its own.
#
#   A Gelman-Rubin-style potential scale reduction factor:
#       W    = mean within-chain variance
#       B    = n / (m - 1) * sum((chain mean - grand mean)^2)
#       Vhat = (n - 1)/n * W + B/n
#       PSRF = sqrt(Vhat / W)
#   Reported over all iterations and over the second half, since early
#   drift inflates the all-iteration value even when the chains have
#   since settled. The second-half value is the one tested.
#
#   Called PSRF rather than R-hat deliberately. The formula is the
#   Gelman-Rubin one, but it is applied to mice chain summaries, not to
#   draws from a Bayesian posterior. It diagnoses how well the
#   chained-equations chains mixed, and the thesis describes it in those
#   terms rather than borrowing the stronger guarantees R-hat carries in
#   its usual setting.
#
# VARIABLES WITH NOTHING TO IMPUTE
#   A variable that was already complete is never imputed, so mice records
#   no chain for it and its chainMean is entirely NaN. In this dataset
#   ipaq_ord is complete. Such variables are reported as "not imputed" and
#   skipped rather than plotted; passing an all-NaN matrix to a plot is
#   what made the first version of this code fail.
#
# INPUT   data/clean/mice_<TAG>_functional.rds
#         data/clean/mice_<TAG>_structural.rds
# OUTPUT  results/sem/mice_convergence_<TAG>_functional.png
#         results/sem/mice_convergence_<TAG>_structural.png
#         results/sem/mice_psrf_<TAG>.csv
#         results/sem/mice_convergence_log_<TAG>.txt
# =====================================================================

TAG  <- Sys.getenv("FUNCTIONAL_TAG", "pca")
ROOT <- "C:/Users/inaye/Desktop/MASTERTHESIS"
IND  <- file.path(ROOT, "data", "clean")
RES  <- file.path(ROOT, "results", "sem")
dir.create(RES, recursive = TRUE, showWarnings = FALSE)

PSRF_MAX   <- 1.05

# Must match the labelling in 03_impute.R: a run with N_ITER other than the
# pre-specified 20 is stored separately. Set N_ITER here to diagnose such a
# run instead of the registered one.
N_ITER <- as.integer(Sys.getenv("N_ITER", "20"))
RUN    <- if (N_ITER == 20L) TAG else sprintf("%s_it%d", TAG, N_ITER)
TRACE_VARS <- c("income_log", "sleep_quality", "perceived_stress",
                "education_ord", "ipaq_ord")

while (sink.number() > 0) sink()
try(sink(type = "message"), silent = TRUE)

con <- file(file.path(RES, sprintf("mice_convergence_log_%s.txt", RUN)), open = "wt")
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

cat("run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", R.version.string, "\n")
cat("tag:", TAG, "\n")


# ---------------------------------------------------------------------
# 1. Load
# ---------------------------------------------------------------------
load_obj <- function(modality) {
  f <- file.path(IND, sprintf("mice_%s_%s.rds", RUN, modality))
  if (!file.exists(f))
    stop("not found: ", f, "\n  Run 03_impute.R first.")
  obj <- readRDS(f)
  if (is.null(obj$mids))
    stop(f, ": no $mids field, cannot diagnose convergence")
  cat(sprintf("\n%s: %d imputations, %d iterations\n", modality,
              obj$mids$m, obj$mids$iteration))
  if (!is.null(obj$provenance$seed))
    cat(sprintf("  imputed %s with seed %s, mice %s\n",
                obj$provenance$run_time, obj$provenance$seed,
                obj$provenance$mice_version))
  obj
}

objs <- list(functional = load_obj("functional"),
             structural = load_obj("structural"))


# ---------------------------------------------------------------------
# 2. Chain arrays
# ---------------------------------------------------------------------
# chainMean and chainVar are arrays indexed [variable, iteration, chain].
# A variable with nothing imputed has an all-NaN slice; those are
# identified here once and handled consistently by both diagnostics.

usable <- function(obj) {
  cm <- obj$mids$chainMean
  if (is.null(cm)) return(character(0))
  present <- intersect(TRACE_VARS, dimnames(cm)[[1]])
  keep <- present[vapply(present,
                         function(v) any(is.finite(cm[v, , ])), logical(1))]
  skipped <- setdiff(TRACE_VARS, keep)
  if (length(skipped))
    cat(sprintf("  %s: not imputed, no chain to diagnose: %s\n",
                obj$modality, paste(skipped, collapse = ", ")))
  keep
}

cat("\n=== variables with chains ===\n")
vars_of <- lapply(objs, usable)


# ---------------------------------------------------------------------
# 3. Trace plots
# ---------------------------------------------------------------------
trace_plot <- function(obj, vars) {
  f <- file.path(RES, sprintf("mice_convergence_%s_%s.png", RUN, obj$modality))
  if (!length(vars)) {
    cat(sprintf("  %-11s no variable has a chain; no plot written\n",
                obj$modality))
    return(invisible(FALSE))
  }
  cm <- obj$mids$chainMean
  cv <- obj$mids$chainVar

  png(f, width = 1200, height = 210 * length(vars), res = 105)
  # The device must close even if a panel fails, otherwise the file is left
  # truncated and the next plot inherits an open device.
  on.exit({ try(dev.off(), silent = TRUE) }, add = TRUE)

  op <- par(mfrow = c(length(vars), 2), mar = c(3.4, 4, 2.4, 0.8),
            mgp = c(2.1, 0.7, 0))
  on.exit(try(par(op), silent = TRUE), add = TRUE, after = FALSE)

  for (v in vars) {
    for (what in c("mean", "sd")) {
      X <- if (what == "mean") cm[v, , ] else sqrt(pmax(cv[v, , ], 0))
      X <- as.matrix(X)
      # Drop chains that are not fully finite rather than letting a single
      # NaN column collapse the axis limits.
      X <- X[, apply(X, 2, function(z) all(is.finite(z))), drop = FALSE]
      if (!ncol(X) || !any(is.finite(X))) {
        plot.new()
        title(main = sprintf("%s  (%s): no finite values", v, what))
        next
      }
      matplot(seq_len(nrow(X)), X, type = "l", lty = 1,
              col = grDevices::rainbow(ncol(X)),
              xlab = "iteration", ylab = what,
              main = sprintf("%s  (%s)", v, what))
    }
  }
  cat(sprintf("  %-11s %s  [%d variables, %d chains]\n",
              obj$modality, f, length(vars), dim(cm)[3]))
  invisible(TRUE)
}

cat("\n=== trace plots ===\n")
for (m in names(objs)) trace_plot(objs[[m]], vars_of[[m]])


# ---------------------------------------------------------------------
# 4. Mixing diagnostic
# ---------------------------------------------------------------------
psrf <- function(X) {
  X <- as.matrix(X)
  X <- X[, apply(X, 2, function(z) all(is.finite(z))), drop = FALSE]
  n <- nrow(X); m <- ncol(X)
  if (m < 2L || n < 2L) return(NA_real_)
  xbar_j <- colMeans(X)
  B <- n / (m - 1) * sum((xbar_j - mean(xbar_j))^2)
  W <- mean(apply(X, 2, stats::var))
  if (!is.finite(W) || W <= 0) return(NA_real_)
  sqrt(((n - 1) / n * W + B / n) / W)
}

psrf_table <- function(obj, vars) {
  if (!length(vars)) return(NULL)
  cm <- obj$mids$chainMean
  n  <- dim(cm)[2]
  half <- seq.int(floor(n / 2) + 1L, n)
  data.frame(
    modality   = obj$modality,
    variable   = vars,
    psrf_all   = vapply(vars, function(v) psrf(cm[v, , , drop = TRUE]), numeric(1)),
    psrf_half  = vapply(vars, function(v) psrf(cm[v, half, , drop = TRUE]), numeric(1)),
    iterations = n,
    chains     = dim(cm)[3],
    row.names  = NULL, stringsAsFactors = FALSE)
}

# Matrix, loaded as a mice dependency, makes rbind an S4 generic whose
# dispatch indexes list(...)[[1]], so do.call(rbind, list()) errors instead
# of returning NULL. Guard every rbind over a collected list.
parts <- Filter(Negate(is.null),
                lapply(names(objs), function(m) psrf_table(objs[[m]], vars_of[[m]])))
rh <- if (length(parts)) do.call(rbind, parts) else NULL

cat("\n=== mixing ===\n")
if (is.null(rh)) {
  cat("  no chains to summarise; inspect the trace plots directly\n")
} else {
  write.csv(rh, file.path(RES, sprintf("mice_psrf_%s.csv", RUN)), row.names = FALSE)
  cat(sprintf("  Gelman-Rubin-style PSRF from chain means, %d iterations, %d chains\n",
              rh$iterations[1], rh$chains[1]))
  cat(sprintf("  threshold %.2f, applied to the second half of the iterations\n\n",
              PSRF_MAX))
  print(rh[order(-rh$psrf_half),
           c("modality", "variable", "psrf_all", "psrf_half")], row.names = FALSE)
  cat("\n  wrote", file.path(RES, sprintf("mice_psrf_%s.csv", RUN)), "\n")

  bad <- rh[is.finite(rh$psrf_half) & rh$psrf_half > PSRF_MAX, ]
  if (nrow(bad)) {
    cat("\n  Supplementary diagnostic above", PSRF_MAX, "for:\n")
    for (i in seq_len(nrow(bad)))
      cat(sprintf("    %s / %s: %.3f\n",
                  bad$modality[i], bad$variable[i], bad$psrf_half[i]))
    cat("\n  This does NOT by itself mean the imputations are unusable. The\n")
    cat("  pre-specified criterion is the trace plots; read them first. Treat\n")
    cat("  a value in this range as a prompt to look, not as a verdict, for\n")
    cat("  the reason given in the header: with 50 short chains the statistic\n")
    cat("  is noisy and biased upward.\n")
    cat("\n  If the plots also show drift, raise N_ITER in 03_impute.R and\n")
    cat("  rerun. Runs with N_ITER other than 20 are written under their own\n")
    cat("  file names, so the pre-specified run is never overwritten.\n")
  } else {
    cat("\n  Every value on the second half is at or below", PSRF_MAX, "\n")
  }
  cat("\n  The pre-specified convergence check is the trace plots. Inspect\n")
  cat("  them before proceeding to 04_h1_base.R.\n")
}

close_log()
