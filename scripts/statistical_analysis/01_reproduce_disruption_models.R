
# ==============================================================
# Objective 3: Independent reproduction of Table 3
# Negative-binomial regression
# ==============================================================
# Input: data/derived/01_model_analysis_data.tsv
# Output: results/model_verification/
#
# This script does not modify the original dataset.
# ==============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(MASS))

input_file <- "data/derived/01_model_analysis_data.tsv"
output_dir <- "results/model_verification"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------
# 1. Read data
# --------------------------------------------------------------

d <- read.delim(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n===== DATASET =====\n")
cat("Observations:", nrow(d), "\n")
cat("Unique isolates:", length(unique(d$sample_id)), "\n")

required <- c(
  "sample_id",
  "lineage",
  "n_primary_disrupted",
  "syn_snps",
  "mean_depth",
  "study_accession"
)

stopifnot(all(required %in% names(d)))
stopifnot(nrow(d) == 3227)
stopifnot(!anyDuplicated(d$sample_id))

# --------------------------------------------------------------
# 2. Prepare variables
# --------------------------------------------------------------

d$lineage <- factor(
  d$lineage,
  levels = paste0("lineage", c(4, 1, 2, 3, 5, 6))
)

d$y <- as.numeric(d$n_primary_disrupted)
d$syn <- as.numeric(d$syn_snps)
d$depth <- as.numeric(d$mean_depth)

stopifnot(!anyNA(d[, c("lineage", "y", "syn", "depth")]))
stopifnot(all(d$y >= 0))
stopifnot(all(d$syn >= 0))
stopifnot(all(d$depth >= 0))

# Centre sequencing depth on the full cohort mean.
d$depth_c <- d$depth - mean(d$depth)

cat("\n===== LINEAGE COUNTS =====\n")
print(table(d$lineage))

cat("\nZero synonymous SNP counts:",
    sum(d$syn == 0), "\n")

# --------------------------------------------------------------
# 3. Crude model
# --------------------------------------------------------------

cat("\n===== MODEL 1: CRUDE =====\n")

m1 <- glm.nb(
  y ~ lineage,
  data = d,
  control = glm.control(maxit = 100)
)

print(summary(m1))

# --------------------------------------------------------------
# 4. Divergence-scaled model
# --------------------------------------------------------------

adjusted <- subset(d, syn > 0)

cat("\n===== MODEL 2: DIVERGENCE-SCALED =====\n")
cat("Observations:", nrow(adjusted), "\n")

stopifnot(nrow(adjusted) == 3225)

m2 <- glm.nb(
  y ~ lineage + depth_c + offset(log(syn)),
  data = adjusted,
  control = glm.control(maxit = 100)
)

print(summary(m2))

# --------------------------------------------------------------
# 5. Within-BioProject sensitivity model
# --------------------------------------------------------------

within <- subset(
  adjusted,
  study_accession == "PRJEB36076"
)

cat("\n===== MODEL 3: PRJEB36076 =====\n")
cat("Observations:", nrow(within), "\n")

print(table(within$lineage))

stopifnot(nrow(within) == 449)

m3 <- glm.nb(
  y ~ lineage + depth_c + offset(log(syn)),
  data = within,
  control = glm.control(maxit = 100)
)

print(summary(m3))

# --------------------------------------------------------------
# 6. Extract rate ratios and Wald confidence intervals
# --------------------------------------------------------------

extract_rr <- function(model, model_name) {

  coefficients <- summary(model)$coefficients

  terms <- grep(
    "^lineage",
    rownames(coefficients),
    value = TRUE
  )

  results <- data.frame(
    model = model_name,
    lineage = sub("^lineagelineage", "L", terms),
    estimate = coefficients[terms, "Estimate"],
    se = coefficients[terms, "Std. Error"],
    p_value = coefficients[terms, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )

  results$RR <- exp(results$estimate)

  results$CI_lower <- exp(
    results$estimate - 1.96 * results$se
  )

  results$CI_upper <- exp(
    results$estimate + 1.96 * results$se
  )

  results
}

results <- rbind(
  extract_rr(m1, "Crude"),
  extract_rr(m2, "Divergence_scaled"),
  extract_rr(m3, "Within_PRJEB36076")
)

cat("\n===== REPRODUCED RATE RATIOS =====\n")

print(
  results[, c(
    "model",
    "lineage",
    "RR",
    "CI_lower",
    "CI_upper",
    "p_value"
  )],
  row.names = FALSE,
  digits = 5
)

# --------------------------------------------------------------
# 7. Model diagnostics
# --------------------------------------------------------------

diagnostics <- data.frame(
  model = c(
    "Crude",
    "Divergence_scaled",
    "Within_PRJEB36076"
  ),

  n = c(
    nobs(m1),
    nobs(m2),
    nobs(m3)
  ),

  theta = c(
    m1$theta,
    m2$theta,
    m3$theta
  ),

  residual_deviance_df = c(
    deviance(m1) / df.residual(m1),
    deviance(m2) / df.residual(m2),
    deviance(m3) / df.residual(m3)
  ),

  converged = c(
    isTRUE(m1$converged),
    isTRUE(m2$converged),
    isTRUE(m3$converged)
  )
)

cat("\n===== MODEL DIAGNOSTICS =====\n")
print(diagnostics, row.names = FALSE)

# --------------------------------------------------------------
# 8. Save results
# --------------------------------------------------------------

write.table(
  results,
  file.path(output_dir, "Table3_reproduced_rate_ratios.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  diagnostics,
  file.path(output_dir, "Table3_model_diagnostics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "R_session_info.txt")
)

cat("\n===== COMPLETE =====\n")
cat("Results saved to:", output_dir, "\n")

