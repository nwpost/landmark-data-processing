# density-based-estimation.R
#
# Purpose:
#   Run density-weighted missing-landmark estimation on one TPS landmark dataset.
#
# General use:
#   Option 1: edit `input_tps` and `output_dir` in the configuration section.
#   Option 2: pass paths from the command line:
#       Rscript density-based-estimation.R path/to/input_landmarks.tps results/density_weighting_run
#
# Input:
#   A TPS file containing 2D or 3D landmark coordinates. Specimen IDs should be
#   stored with ID= lines. Missing landmarks may be present.
#
# Output:
#   The script writes diagnostic tables, figures, and complete estimated TPS
#   configurations to `output_dir`.
#
# Notes:
#   - This script intentionally outputs complete estimated configurations.
#   - It does NOT restore any original missing landmark values.
#   - Density weights are proportional to local neighborhood radius:
#       smaller radius = denser morphospace region = lower weight
#       larger radius  = more isolated morphology  = higher weight

# -----------------------------------------------------------------------------
# 1. Configuration
# -----------------------------------------------------------------------------

set.seed(1)

args <- commandArgs(trailingOnly = TRUE)

input_tps <- if (length(args) >= 1) {
  args[[1]]
} else {
  "path/to/your/data.tps"
}

output_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path("results", "density_weighting_run")
}

# Analysis settings.
min_shared_landmarks <- 3
smacof_dimensions <- 10
nmds_dimensions <- 2
nmds_trymax <- 50
random_seed <- 123

# -----------------------------------------------------------------------------
# 2. Load packages and create output folders
# -----------------------------------------------------------------------------

library(geomorph)
library(smacof)
library(vegan)
library(FNN)
library(missMDA)
library(MASS)

if (!file.exists(input_tps)) {
  stop(
    "Input TPS file not found: ", input_tps, "\n\n",
    "Edit `input_tps` near the top of this script or run with command-line arguments, e.g.:\n",
    "  Rscript density-based-estimation.R path/to/input_landmarks.tps results/density_weighting_run",
    call. = FALSE
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "processed"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 3. Helper functions
# -----------------------------------------------------------------------------

hopkins_statistic <- function(coords, m = NULL, seed = 123) {
  # Computes Hopkins statistic using nearest-neighbor distances.
  # H near 0.5 is approximately random under this implementation.
  # Larger H values indicate stronger clustering tendency.
  
  set.seed(seed)
  
  coords <- as.matrix(coords)
  n <- nrow(coords)
  d <- ncol(coords)
  
  if (n < 4) {
    return(list(H = NA_real_, H_complement = NA_real_, m = NA_integer_))
  }
  
  if (is.null(m)) {
    m <- max(1, min(n - 1, floor(0.1 * n)))
  }
  
  sampled_indices <- sample(seq_len(n), size = m, replace = FALSE)
  sampled_points <- coords[sampled_indices, , drop = FALSE]
  
  mins <- apply(coords, 2, min, na.rm = TRUE)
  maxs <- apply(coords, 2, max, na.rm = TRUE)
  
  random_points <- matrix(
    stats::runif(m * d, min = 0, max = 1),
    nrow = m,
    ncol = d
  )
  
  for (j in seq_len(d)) {
    random_points[, j] <- mins[j] + random_points[, j] * (maxs[j] - mins[j])
  }
  
  # Distance from random points to nearest observed specimen.
  u <- FNN::get.knnx(coords, random_points, k = 1)$nn.dist[, 1]
  
  # Distance from sampled observed specimens to nearest other observed specimen.
  w <- FNN::get.knnx(coords, sampled_points, k = 2)$nn.dist[, 2]
  
  H <- sum(u^d) / (sum(u^d) + sum(w^d))
  
  list(
    H = H,
    H_complement = 1 - H,
    m = m
  )
}

# -----------------------------------------------------------------------------
# 4. Load unaligned landmark data
# -----------------------------------------------------------------------------

combined <- geomorph::readland.tps(
  input_tps,
  specID = "ID"
)

dim(combined)  # p × k × n

p <- dim(combined)[1]  # landmarks
k <- dim(combined)[2]  # dimensions
n <- dim(combined)[3]  # specimens
specimen_ids <- dimnames(combined)[[3]]

if (is.null(specimen_ids) || any(is.na(specimen_ids)) || any(specimen_ids == "")) {
  specimen_ids <- paste0("specimen_", seq_len(n))
  dimnames(combined)[[3]] <- specimen_ids
}

cat("Input TPS:", input_tps, "
")
cat("Output directory:", output_dir, "
")
cat("Dimensions:", paste(dim(combined), collapse = " × "), "
")
cat("Specimens with missing landmarks:",
    sum(apply(combined, 3, function(x) any(!stats::complete.cases(x)))), "
")

utils::write.csv(
  data.frame(
    specimen_id = specimen_ids,
    n_missing_landmarks = apply(combined, 3, function(x) sum(!stats::complete.cases(x))),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "tables", "input_missingness_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 5. Diagnostic setup from exploratory script
# -----------------------------------------------------------------------------

cat("
==============================
")
cat(" DENSITY-WEIGHTING DIAGNOSTIC
")
cat("==============================

")

# Check specimens that have too few shared landmarks with every other specimen.
pair_shared <- matrix(0, n, n)
dimnames(pair_shared) <- list(specimen_ids, specimen_ids)

if (n > 1) {
  for (i in seq_len(n - 1)) {
    Xi <- combined[, , i]
    for (j in (i + 1):n) {
      Xj <- combined[, , j]
      mask <- rowSums(is.na(Xi) | is.na(Xj)) == 0
      pair_shared[i, j] <- pair_shared[j, i] <- sum(mask)
    }
  }
}

diag(pair_shared) <- p

max_shared_to_any_other <- apply(pair_shared, 1, function(x) max(x, na.rm = TRUE))

cat("1. Specimens with ZERO shared landmarks to all others:
")
print(names(max_shared_to_any_other[max_shared_to_any_other == 0]))

cat("
==============================
")
cat(" END OF DIAGNOSTIC SETUP
")
cat("==============================

")

# -----------------------------------------------------------------------------
# 6. Pairwise Procrustes distances + shared-landmark counts
# -----------------------------------------------------------------------------

pairwise_procrustes_with_shared <- function(coords_array) {
  stopifnot(length(dim(coords_array)) == 3)
  
  n <- dim(coords_array)[3]
  p <- dim(coords_array)[1]
  D <- matrix(NA_real_, n, n)
  shared <- matrix(NA_real_, n, n)
  
  diag(D) <- 0
  diag(shared) <- p
  
  if (n > 1) {
    for (i in seq_len(n - 1)) {
      Xi <- coords_array[, , i]
      
      for (j in (i + 1):n) {
        Xj <- coords_array[, , j]
        
        # Landmarks present in both specimens: no NA in either specimen.
        mask <- rowSums(is.na(Xi) | is.na(Xj)) == 0
        m_ij <- sum(mask)
        shared[i, j] <- shared[j, i] <- m_ij
        
        if (m_ij >= min_shared_landmarks) {
          Xi_sub <- Xi[mask, , drop = FALSE]
          Xj_sub <- Xj[mask, , drop = FALSE]
          
          # Center to centroid.
          Xi_c <- scale(Xi_sub, scale = FALSE)
          Xj_c <- scale(Xj_sub, scale = FALSE)
          
          # Scale to unit centroid size.
          Xi_size <- sqrt(sum(Xi_c^2))
          Xj_size <- sqrt(sum(Xj_c^2))
          
          if (Xi_size > 0 && Xj_size > 0) {
            Xi_c <- Xi_c / Xi_size
            Xj_c <- Xj_c / Xj_size
            
            # Optimal rotation via SVD.
            svd_out <- svd(t(Xi_c) %*% Xj_c)
            R <- svd_out$u %*% t(svd_out$v)
            
            # Enforce right-handed rotation if needed.
            if (det(R) < 0) {
              svd_out$u[, ncol(svd_out$u)] <- -svd_out$u[, ncol(svd_out$u)]
              R <- svd_out$u %*% t(svd_out$v)
            }
            
            Xj_rot <- Xj_c %*% R
            
            # Procrustes distance.
            D[i, j] <- D[j, i] <- sqrt(sum((Xi_c - Xj_rot)^2))
          }
        }
      }
    }
  }
  
  dimnames(D) <- list(dimnames(coords_array)[[3]], dimnames(coords_array)[[3]])
  dimnames(shared) <- dimnames(D)
  
  list(D = D, shared = shared)
}

pw <- pairwise_procrustes_with_shared(combined)
D_raw <- pw$D
shared_mat <- pw$shared

# Correct distances by shared landmark fraction:
# d*_ij = d_ij * sqrt(p / m_ij)
D_corrected <- D_raw
valid_pairs <- !is.na(D_raw) & !is.na(shared_mat) & shared_mat >= min_shared_landmarks
D_corrected[valid_pairs] <- D_raw[valid_pairs] * sqrt(p / shared_mat[valid_pairs])
D_corrected[!valid_pairs] <- NA_real_
diag(D_corrected) <- 0

utils::write.csv(
  as.data.frame(as.table(D_raw)),
  file.path(output_dir, "tables", "pairwise_procrustes_distances_raw.csv"),
  row.names = FALSE
)

utils::write.csv(
  as.data.frame(as.table(D_corrected)),
  file.path(output_dir, "tables", "pairwise_procrustes_distances_corrected.csv"),
  row.names = FALSE
)

utils::write.csv(
  as.data.frame(as.table(shared_mat)),
  file.path(output_dir, "tables", "pairwise_shared_landmarks.csv"),
  row.names = FALSE
)

cat("Raw distances summary:
")
print(summary(as.vector(D_raw)))
cat("
Corrected distances summary:
")
print(summary(as.vector(D_corrected)))
cat("
Shared landmarks summary:
")
print(summary(as.vector(shared_mat)))

# -----------------------------------------------------------------------------
# 7. SMACOF MDS: smooth corrected partial distances
# -----------------------------------------------------------------------------

# SMACOF is used to smooth corrected partial distances. Missing pairwise
# distances are assigned zero weight explicitly, rather than passed as NA.

set.seed(random_seed)

smacof_distance <- D_corrected
smacof_weight <- ifelse(is.na(smacof_distance), 0, 1)
smacof_distance[is.na(smacof_distance)] <- 0
diag(smacof_distance) <- 0
diag(smacof_weight) <- 0

smacof_fit <- smacof::smacofSym(
  smacof_distance,
  type      = "ratio",
  ndim      = smacof_dimensions,
  itmax     = 500,
  weightmat = smacof_weight,
  verbose   = FALSE
)

coords_smacof <- smacof_fit$conf          # n × ndim coordinates
rownames(coords_smacof) <- specimen_ids
D_smooth      <- dist(coords_smacof)      # NA-free, metric distance matrix
D_smooth_mat  <- as.matrix(D_smooth)

utils::write.csv(
  data.frame(specimen_id = specimen_ids, coords_smacof, check.names = FALSE),
  file.path(output_dir, "tables", "smacof_coordinates.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 8. Hopkins statistic: first diagnostic on the SMACOF morphospace
# -----------------------------------------------------------------------------

hopkins <- hopkins_statistic(coords_smacof, seed = random_seed)

cat("
============================================
")
cat(" Hopkins Statistic
")
cat("============================================

")
cat("Hopkins H:", hopkins$H, "
")
cat("Hopkins 1-H:", hopkins$H_complement, "
")
cat("Sampled points:", hopkins$m, "

")

utils::write.csv(
  data.frame(
    metric = c("Hopkins_H", "Hopkins_1_minus_H", "m_sampled_points"),
    value = c(hopkins$H, hopkins$H_complement, hopkins$m)
  ),
  file.path(output_dir, "tables", "hopkins_statistic.csv"),
  row.names = FALSE
)

# Pre-alignment NMDS on smoothed distances.
set.seed(1)
nmds_pre <- vegan::metaMDS(
  D_smooth,
  k = 2,
  autotransform = FALSE,
  trymax = nmds_trymax,
  trace = FALSE
)

coords_pre <- nmds_pre$points
cat("Pre-alignment NMDS stress (smoothed distances):", nmds_pre$stress, "
")

# -----------------------------------------------------------------------------
# 9. Density-based weights from SMACOF configuration
# -----------------------------------------------------------------------------

coords <- coords_smacof
k_local <- max(3, round(sqrt(nrow(coords))))
if (k_local >= nrow(coords)) {
  k_local <- nrow(coords) - 1
}

nn <- FNN::get.knn(coords, k = k_local)
local_radius <- rowMeans(nn$nn.dist)

# IMPORTANT:
#   local_radius is smaller in dense regions and larger in isolated regions.
#   Therefore weights should be proportional to local_radius, not 1/local_radius,
#   if the goal is to downweight dense regions and upweight isolated morphologies.
w_density <- local_radius
invalid <- !is.finite(w_density) | w_density <= 0
if (any(invalid)) {
  message("Replaced ", sum(invalid), " invalid weights with 1 (neutral weighting).")
  w_density[invalid] <- 1
}

# Normalize to mean 1.
w_density <- w_density / mean(w_density)

# Use the original downstream name.
w <- w_density

cat("Density weights summary:
")
print(summary(w))

utils::write.csv(
  data.frame(
    specimen_id = specimen_ids,
    local_radius = local_radius,
    density_weight = w,
    sqrt_density_weight = sqrt(w),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "tables", "density_weights.csv"),
  row.names = FALSE
)

utils::write.csv(
  data.frame(
    metric = c(
      "k_local",
      "min_local_radius",
      "max_local_radius",
      "min_weight",
      "max_weight",
      "mean_weight",
      "cor_weight_local_radius"
    ),
    value = c(
      k_local,
      min(local_radius, na.rm = TRUE),
      max(local_radius, na.rm = TRUE),
      min(w, na.rm = TRUE),
      max(w, na.rm = TRUE),
      mean(w, na.rm = TRUE),
      stats::cor(w, local_radius, use = "complete.obs")
    )
  ),
  file.path(output_dir, "tables", "density_weight_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 10. Weighted PCA-based imputation with missMDA
# -----------------------------------------------------------------------------

flatten_to_matrix <- function(arr) {
  stopifnot(length(dim(arr)) == 3)
  p <- dim(arr)[1]
  k <- dim(arr)[2]
  n <- dim(arr)[3]
  mat <- array(arr, dim = c(p * k, n))
  t(mat)  # n × (p*k)
}

X <- flatten_to_matrix(combined)  # n × q

sqrt_w <- sqrt(w)

# Row-wise weighting, broadcast along columns.
X_weighted <- X * sqrt_w[row(X)]

set.seed(random_seed)
nb_w <- missMDA::estim_ncpPCA(X_weighted, method.cv = "Kfold", nbsim = 10)

imputed_weighted <- missMDA::imputePCA(
  X_weighted,
  ncp = nb_w$ncp
)

# Undo row scaling.
X_imputed <- imputed_weighted$completeObs / sqrt_w[row(X)]

# Back to p × k × n array.
estimated_weighted <- array(
  data = t(X_imputed),
  dim  = c(p, k, n),
  dimnames = dimnames(combined)
)

dim(estimated_weighted)

# Unweighted PCA imputation for comparison.
set.seed(random_seed)
nb_unw <- missMDA::estim_ncpPCA(X, method.cv = "Kfold", nbsim = 10)
imputed_unw <- missMDA::imputePCA(X, ncp = nb_unw$ncp)

estimated_unweighted <- array(
  data = t(imputed_unw$completeObs),
  dim  = c(p, k, n),
  dimnames = dimnames(combined)
)

utils::write.csv(
  data.frame(
    analysis = c("density_weighted_PCA", "unweighted_PCA"),
    ncp = c(nb_w$ncp, nb_unw$ncp)
  ),
  file.path(output_dir, "tables", "pca_imputation_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 11. GPA + PCA for weighted vs unweighted datasets
# -----------------------------------------------------------------------------

gpa_unw <- geomorph::gpagen(estimated_unweighted, PrinAxes = FALSE)
gpa_w   <- geomorph::gpagen(estimated_weighted,   PrinAxes = FALSE)

pca_unw <- geomorph::gm.prcomp(gpa_unw$coords)
pca_w   <- geomorph::gm.prcomp(gpa_w$coords)

# -----------------------------------------------------------------------------
# 12. Plots from exploratory script
# -----------------------------------------------------------------------------

pdf(file.path(output_dir, "figures", "pca_weighted_unweighted_side_by_side.pdf"), width = 9, height = 4.5)
par(mfrow = c(1, 2))

plot(
  pca_unw$x[, 1], pca_unw$x[, 2],
  pch = 21, bg = "gray80",
  xlab = "PC1 (unweighted)", ylab = "PC2 (unweighted)",
  main = "Unweighted GPA + PCA"
)

plot(
  pca_w$x[, 1], pca_w$x[, 2],
  pch = 21, bg = "gray80",
  xlab = "PC1 (weighted)", ylab = "PC2 (weighted)",
  main = "Weighted GPA + PCA"
)

par(mfrow = c(1, 1))
dev.off()

pdf(file.path(output_dir, "figures", "consensus_difference_weighted_vs_unweighted.pdf"), width = 6, height = 6)
geomorph::plotRefToTarget(
  gpa_unw$consensus,
  gpa_w$consensus,
  method = "vector",
  mag = 2,
  main = "Consensus Difference: Weighted vs Unweighted"
)
dev.off()

# NMDS on tangent-space distances after GPA.
X_unw <- geomorph::two.d.array(gpa_unw$coords)
X_w   <- geomorph::two.d.array(gpa_w$coords)

D_unw <- dist(X_unw)
D_w   <- dist(X_w)

set.seed(1)
nmds_unw <- vegan::metaMDS(D_unw, k = nmds_dimensions, autotransform = FALSE, trymax = nmds_trymax, trace = FALSE)
nmds_w   <- vegan::metaMDS(D_w,   k = nmds_dimensions, autotransform = FALSE, trymax = nmds_trymax, trace = FALSE)

# Align weighted NMDS to unweighted NMDS before interpreting overlay/displacement.
# This preserves the original NMDS comparison but avoids arbitrary rotation or reflection.
nmds_proc <- vegan::procrustes(
  X = nmds_unw$points,
  Y = nmds_w$points,
  symmetric = TRUE
)

coords_unw <- nmds_proc$X
coords_w   <- nmds_proc$Yrot

utils::write.csv(
  data.frame(
    specimen_id = rep(specimen_ids, times = 2),
    configuration = rep(c("unweighted", "density_weighted_aligned"), each = n),
    NMDS1 = c(coords_unw[, 1], coords_w[, 1]),
    NMDS2 = c(coords_unw[, 2], coords_w[, 2]),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "tables", "nmds_weighted_unweighted_coordinates.csv"),
  row.names = FALSE
)

pdf(file.path(output_dir, "figures", "nmds_weighted_unweighted_side_by_side.pdf"), width = 9, height = 4.5)
par(mfrow = c(1, 2))

plot(
  nmds_unw$points[, 1], nmds_unw$points[, 2],
  pch = 21, bg = "gray80",
  xlab = "NMDS 1", ylab = "NMDS 2",
  main = paste0("Unweighted NMDS (stress = ", round(nmds_unw$stress, 3), ")")
)

plot(
  nmds_w$points[, 1], nmds_w$points[, 2],
  pch = 21, bg = "gray80",
  xlab = "NMDS 1", ylab = "NMDS 2",
  main = paste0("Weighted NMDS (stress = ", round(nmds_w$stress, 3), ")")
)

par(mfrow = c(1, 1))
dev.off()

pdf(file.path(output_dir, "figures", "nmds_weighted_unweighted_overlay_aligned.pdf"), width = 7, height = 5)
plot(
  coords_unw[, 1], coords_unw[, 2],
  pch = 21, bg = "gray80",
  xlab = "NMDS1", ylab = "NMDS2",
  main = "Weighted (red) vs Unweighted (gray)"
)

points(
  coords_w[, 1], coords_w[, 2],
  pch = 21,
  bg = adjustcolor("red", alpha.f = 0.6)
)

segments(
  x0 = coords_unw[, 1], y0 = coords_unw[, 2],
  x1 = coords_w[, 1],   y1 = coords_w[, 2],
  col = adjustcolor("red", alpha.f = 0.4)
)

legend(
  "topright",
  legend = c("Unweighted", "Weighted"),
  pt.bg  = c("gray80", adjustcolor("red", alpha.f = 0.6)),
  pch = 21,
  bty = "n"
)
dev.off()

# Write weighted-imputed and unweighted-imputed TPS files.
geomorph::writeland.tps(
  A    = estimated_weighted,
  file = file.path(output_dir, "processed", "estimated_density_weighted_complete_unaligned.tps"),
  scale  = NULL,
  specID = TRUE
)

geomorph::writeland.tps(
  A    = estimated_unweighted,
  file = file.path(output_dir, "processed", "estimated_unweighted_complete_unaligned.tps"),
  scale  = NULL,
  specID = TRUE
)

geomorph::writeland.tps(
  A    = gpa_w$coords,
  file = file.path(output_dir, "processed", "estimated_density_weighted_complete_GPA_aligned.tps"),
  scale  = NULL,
  specID = TRUE
)

geomorph::writeland.tps(
  A    = gpa_unw$coords,
  file = file.path(output_dir, "processed", "estimated_unweighted_complete_GPA_aligned.tps"),
  scale  = NULL,
  specID = TRUE
)

# -----------------------------------------------------------------------------
# 13. Diagnostics from exploratory script
# -----------------------------------------------------------------------------

cat("
============================================
")
cat(" 1. Distance Matrix Diagnostics
")
cat("============================================

")

cat("Raw distances summary:
")
print(summary(as.vector(D_raw)))

cat("
Corrected distances summary:
")
print(summary(as.vector(D_corrected)))

cat("
Shared landmarks per pair:
")
print(summary(as.vector(shared_mat)))

pdf(file.path(output_dir, "figures", "hist_shared_landmarks.pdf"), width = 6, height = 4)
hist(shared_mat,
     breaks = 40,
     main = "Histogram: Number of Shared Landmarks",
     xlab = "Shared landmarks",
     col = "gray80")
dev.off()

cat("
SMACOF stress:", smacof_fit$stress, "
")

D_smooth_vals <- as.vector(D_smooth_mat)
cat("Smoothed distances summary:
")
print(summary(D_smooth_vals))

cat("
============================================
")
cat(" 2. Weighting Diagnostics (Density-Based)
")
cat("============================================

")

cat("Density-based weight summary:
")
print(summary(w))

pdf(file.path(output_dir, "figures", "density_based_weights.pdf"), width = 9, height = 4.5)
par(mfrow = c(1, 2))
hist(
  w,
  breaks = 25,
  main = "Density-Based Weights",
  xlab = "Weight",
  col = "gray80"
)

plot(
  w,
  pch = 21, bg = "gray70",
  ylab = "Weight",
  xlab = "Specimen index",
  main = "Specimen Weights (density-based)"
)
abline(h = 1, col = "red", lty = 2)
par(mfrow = c(1, 1))
dev.off()

cat("
Correlation (weight vs local radius):
")
cat("Expected: strongly positive
")
print(cor(w, local_radius))

cat("
============================================
")
cat(" 3. PCA Imputation Diagnostics
")
cat("============================================

")

cat("Remaining NA in weighted-imputed matrix:", sum(is.na(estimated_weighted)), "
")
cat("Remaining NA in unweighted-imputed matrix:", sum(is.na(estimated_unweighted)), "
")

cat("
Unweighted PCA variance (first 5 PCs):
")
print(summary(pca_unw)$importance[, 1:min(5, ncol(pca_unw$x))])

cat("
Weighted PCA variance (first 5 PCs):
")
print(summary(pca_w)$importance[, 1:min(5, ncol(pca_w$x))])

cat("
============================================
")
cat(" 4. GPA Diagnostics
")
cat("============================================

")

cat("Centroid sizes (unweighted):
")
print(summary(gpa_unw$Csize))

cat("
Centroid sizes (weighted):
")
print(summary(gpa_w$Csize))

cat("
============================================
")
cat(" 5. Shape Space Diagnostics (NMDS)
")
cat("============================================

")

cat("Unweighted NMDS stress:", round(nmds_unw$stress, 3), "
")
cat("Weighted NMDS stress:",   round(nmds_w$stress, 3), "
")
cat("NMDS Procrustes sum of squares:", nmds_proc$ss, "

")

cat("
============================================
")
cat(" 6. Cluster Structure Comparison
")
cat("============================================

")

kde_unw <- MASS::kde2d(nmds_unw$points[, 1], nmds_unw$points[, 2])
kde_w   <- MASS::kde2d(nmds_w$points[, 1],   nmds_w$points[, 2])

cluster_ratio <- max(kde_unw$z) / max(kde_w$z)

cat("Cluster density ratio (unweighted / weighted):", cluster_ratio, "
")
cat("(>1 means weighting reduced dominance of dense clusters)

")

utils::write.csv(
  data.frame(
    metric = c(
      "Hopkins_H",
      "Hopkins_1_minus_H",
      "Hopkins_sampled_points",
      "smacof_stress",
      "pre_alignment_nmds_stress",
      "unweighted_nmds_stress",
      "weighted_nmds_stress",
      "nmds_procrustes_ss",
      "cluster_density_ratio_unweighted_over_weighted",
      "remaining_NA_weighted_imputed",
      "remaining_NA_unweighted_imputed"
    ),
    value = c(
      hopkins$H,
      hopkins$H_complement,
      hopkins$m,
      smacof_fit$stress,
      nmds_pre$stress,
      nmds_unw$stress,
      nmds_w$stress,
      nmds_proc$ss,
      cluster_ratio,
      sum(is.na(estimated_weighted)),
      sum(is.na(estimated_unweighted))
    )
  ),
  file.path(output_dir, "tables", "density_weighting_diagnostics_summary.csv"),
  row.names = FALSE
)

utils::capture.output(
  utils::sessionInfo(),
  file = file.path(output_dir, "sessionInfo_density_weighting.txt")
)

cat("
============================================
")
cat(" Diagnostics complete.
")
cat(" Output written to: ", output_dir, "
", sep = "")
cat(" To use another dataset, replace `input_tps` and `output_dir` at the top or pass them as command-line arguments.
")
cat("============================================

")
