# =====================================================================
# MULTIPLICITY FOR THE CONFIRMATORY (JOINT) MODELS OF H2 AND H3
#
# The confirmatory H2 and H3 models enter all six exposures together
# (05b_h2_joint.R, 06b_h3_joint.R). This script applies the
# Benjamini-Hochberg correction to those models:
#
#   H2  12 tests: six exposures x two modalities, outcome = brain-age gap
#   H3  12 tests: six exposures x two modalities, outcome = latent g
#   q = .05, two separate families (H2 and H3 ask different questions).
#
# One test per exposure. Physical activity enters the models as two
# indicators (moderate and high vs low); the tested contrast is high vs
# low (ipaq_high). ipaq_mod keeps low as the reference and is not tested.
#
# H1 (Bonferroni over two modalities) is unchanged and stays in
# 08_multiplicity.R. That script corrects the single-exposure H2/H3 models,
# which are reported as supplementary; this script corrects the joint
# models, which are confirmatory. Both lead to the same conclusion.
#
# The q values written here are the ones shown in the H2 and H3 figures
# (fig_h2_exposures_warm.py, fig_h3_exposures_warm.py compute the same
# correction from the same rows).
#
# HOW TO RUN
#   Rscript code/sem/08b_multiplicity_joint.R
#   (after 05b_h2_joint.R and 06b_h3_joint.R)
#
# READS   results/sem/h2_joint_estimates_<TAG>.csv
#         results/sem/h3_joint_estimates_<TAG>.csv
# WRITES  results/sem/multiplicity_joint_<TAG>.csv
#         results/sem/multiplicity_joint_log_<TAG>.txt
#         (aggregate estimates only; no participant-level output)
# =====================================================================

TAG  <- Sys.getenv("FUNCTIONAL_TAG", "pca_it50")
ROOT <- Sys.getenv("PROJECT_ROOT", "C:/Users/inaye/Desktop/MASTERTHESIS")
RES  <- file.path(ROOT, "results", "sem")
Q_LEVEL <- 0.05

SIX <- c(income_z = "household income", sleep_z = "sleep quality",
         ipaq_high = "physical activity (high vs low)",
         PWB = "psychological well-being", SOC = "social connectedness",
         stress_z = "perceived stress")

while (sink.number() > 0) sink()
con <- file(file.path(RES, sprintf("multiplicity_joint_log_%s.txt", TAG)), open = "wt")
sink(con, split = TRUE)
on.exit({ while (sink.number() > 0) sink(); try(close(con), silent = TRUE) }, add = TRUE)

cat("=== multiplicity, joint (confirmatory) models ===\n")
cat("  run ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("  tag ", TAG, "\n")

read_family <- function(hyp, file, type, outcome) {
  f <- file.path(RES, file)
  if (!file.exists(f)) stop("not found: ", f, "\n  Run the joint model script first.")
  e <- read.csv(f, stringsAsFactors = FALSE)
  cat(sprintf("\n  %s: %s\n    md5 %s\n", hyp, basename(f), unname(tools::md5sum(f))))
  x <- e[e$type == type & e$op == "~" & e$lhs == outcome & e$rhs %in% names(SIX), ]
  if (nrow(x) != 12)
    stop(hyp, ": expected 12 exposure rows (6 x 2 modalities), found ", nrow(x))
  data.frame(family = hyp, modality = x$modality, term = x$rhs,
             label = unname(SIX[x$rhs]), estimate = x$est, se = x$se,
             ci.lower = x$ci.lower, ci.upper = x$ci.upper, p = x$pvalue,
             stringsAsFactors = FALSE)
}

# H2: estimates in years of brain-age gap per SD of exposure.
# H3: standardised beta on latent general cognition.
# The p value of a coefficient does not depend on its scaling.
h2 <- read_family("H2", sprintf("h2_joint_estimates_%s.csv", TAG),
                  "unstandardized_years", "BAG")
h3 <- read_family("H3", sprintf("h3_joint_estimates_%s.csv", TAG),
                  "standardized", "g")

correct <- function(x) {
  x$q <- p.adjust(x$p, method = "BH")
  x$decision <- ifelse(x$q < Q_LEVEL, "supported", "not supported")
  x[order(x$p), ]
}
out <- rbind(correct(h2), correct(h3))

for (hyp in c("H2", "H3")) {
  y <- out[out$family == hyp, ]
  cat(sprintf("\n=== %s: 12 tests, Benjamini-Hochberg at q = %.2f ===\n", hyp, Q_LEVEL))
  show <- data.frame(modality = y$modality, exposure = y$label,
                     estimate = round(y$estimate, 3),
                     ci = sprintf("[%.3f, %.3f]", y$ci.lower, y$ci.upper),
                     p = signif(y$p, 3), q = signif(y$q, 3), decision = y$decision)
  print(show, row.names = FALSE)
}

write.csv(out, file.path(RES, sprintf("multiplicity_joint_%s.csv", TAG)), row.names = FALSE)
cat("\nwrote ", file.path(RES, sprintf("multiplicity_joint_%s.csv", TAG)), "\n")
