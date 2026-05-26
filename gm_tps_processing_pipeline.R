# gm_tps_processing_pipeline.R
#
# Purpose:
#   General TPS landmark processing pipeline for preparing specimen-level or taxon-level geometric morphometric landmark data for downstream analysis.

# This script consolidates the following operations into one generalized workflow:
#   1. Read and merge one or more TPS files.
#   2. Optionally use geomorph::bilat.symmetry() to estimate one-sided missing bilateral landmarks.
#   3. Estimate remaining missing landmarks, run GPA, and restore the missing data mask.
#   4. Optionally collapse specimens into taxon/group means using a group map.
#   5. Optionally realign the final group level dataset and restore group level missing data.
#
# Input requirements:
#   - TPS landmark files with ID= specimen identifiers.
#   - Optional group map CSV with columns: group, specimen_id.
#   - Optional bilateral landmark-pair CSV/XLSX with columns named below in the configuration.
#
# Notes:
#   - Missing landmarks are temporarily estimated only to permit bilateral symmetry estimation, GPA, and/or mean-shape calculation.
#   - If bilat.symmetry is used, only one-sided missing bilateral landmarks are filled from the symmetric component. Observed landmarks are not replaced, and landmarks missing on both sides of a pair remain missing.
#   - The remaining missing-data mask is restored after global GPA alignment.
#   - If taxon/group means are calculated, a landmark is restored to NA for a group when that landmark is missing in all specimens assigned to that group.

# -----------------------------------------------------------------------------
# 1. Configuration
# -----------------------------------------------------------------------------

# Directory containing raw TPS files. If `input_tps_files` is NULL, all .tps files in this directory will be read.
input_tps_dir <- file.path("data", "raw_tps")

# Optional explicit file list. Use NULL to read all .tps files in input_tps_dir.
# Example:
# input_tps_files <- c(
#   file.path("data", "raw_tps", "specimen_set_1.tps"),
#   file.path("data", "raw_tps", "specimen_set_2.tps")
# )
input_tps_files <- NULL

# Output directory.
output_dir <- file.path("results", "gm_tps_processing_pipeline")

# Optional symmetry estimation of one-sided missing bilateral landmarks.
# Use NULL to skip. Accepted file types: .csv, .xlsx, .xls.
# The pair file should identify corresponding left/right landmarks across the symmetry plane and include the columns named below.
# Example:
# bilateral_pairs_file <- file.path("metadata", "bilateral_pairs.csv")
bilateral_pairs_file <- NULL

# Names of columns in bilateral_pairs_file.
bilateral_landmark_col <- "landmark"
bilateral_pair_col <- "pair"

# Temporary missing-data estimation used only to create a complete array for geomorph::bilat.symmetry(). The original missing data mask is preserved, and only one sided missing bilateral landmarks are filled from the symmetric component.
# Options: "pca" or "TPS".
bilateral_imputation_method <- "TPS"

# Number of permutations used by bilat.symmetry(). This is not used here for hypothesis testing but a low value keeps the preprocessing step fast while still producing the symmetric component.
bilateral_symmetry_iter <- 1

# Missing-data estimation used before global GPA alignment.
# Options: "pca" or "TPS".
alignment_imputation_method <- "pca"

# Missing-data estimation used only for group/taxon mean-shape calculation.
# TPS interpolation is often appropriate for representative mean shapes.
# Options: "pca" or "TPS".
group_mean_imputation_method <- "TPS"

# Optional group/taxon map. CSV must include columns: group, specimen_id.
# Use NULL to skip group-mean creation and output only specimen-level aligned data.
# Example:
# group_map_file <- file.path("metadata", "group_map.csv")
group_map_file <- NULL

# If TRUE, after group means are created, estimate missing landmarks, GPA-align the
# group-level dataset, and restore the group-level missing-data mask.
realign_group_means <- TRUE

# Clean final group/taxon names. Example: remove an extant mean suffix.
clean_final_names <- TRUE
final_name_pattern <- "-mfavg$"
final_name_replacement <- ""

# -----------------------------------------------------------------------------
# 2. Packages and output folders
# -----------------------------------------------------------------------------

library(geomorph)
library(abind)
library(missMDA)
library(readxl)

set.seed(123)

# Output directory structure:
#   output_dir/
#     Main results directory for the full pipeline run.
#
#   output_dir/intermediate/
#     TPS files produced during processing, including the merged raw TPS, the post-bilateral-symmetry TPS, the globally GPA-aligned TPS with missing values restored, and the optional group/taxon mean TPS.
#
#   output_dir/group_means/
#     Individual TPS files for each group/taxon mean when `group_map_file` is supplied. This folder is created even when group means are skipped.
#
#   output_dir/ready_for_analysis/
#     Final TPS outputs intended for downstream analyses. The timestamped file is preserved for record keeping, and `ready_for_analysis_latest.tps` is overwritten on each run for convenience.
#
#   output_dir/tables/
#     CSV summaries documenting missingness, bilateral-symmetry filling, imputation settings, group composition, final output paths, and session info.

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "intermediate"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "group_means"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "ready_for_analysis"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 3. Helper functions
# -----------------------------------------------------------------------------

message_step <- function(text) {
  message("
--- ", text, " ---")
}

read_table_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  
  if (ext == "csv") {
    return(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  
  if (ext %in% c("xlsx", "xls")) {
    return(as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE))
  }
  
  stop("Unsupported table file type: ", path, call. = FALSE)
}

get_tps_paths <- function(input_tps_dir, input_tps_files = NULL) {
  if (!is.null(input_tps_files)) {
    paths <- input_tps_files
  } else {
    paths <- list.files(
      input_tps_dir,
      pattern = "[.]tps$",
      full.names = TRUE,
      ignore.case = TRUE
    )
  }
  
  if (length(paths) == 0) {
    stop("No TPS files found. Check `input_tps_dir` or `input_tps_files`.", call. = FALSE)
  }
  
  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths) > 0) {
    stop("Missing TPS file(s):
", paste(missing_paths, collapse = "
"), call. = FALSE)
  }
  
  paths
}

read_tps_array <- function(path) {
  geomorph::readland.tps(path, specID = "ID")
}

merge_tps_arrays <- function(paths) {
  arrays <- vector("list", length(paths))
  
  for (i in seq_along(paths)) {
    message("Reading TPS file: ", paths[i])
    arrays[[i]] <- read_tps_array(paths[i])
  }
  
  dims <- lapply(arrays, dim)
  pk <- do.call(rbind, lapply(dims, function(x) x[1:2]))
  
  if (length(unique(pk[, 1])) != 1 || length(unique(pk[, 2])) != 1) {
    stop(
      "All TPS files must have the same landmark count and dimensionality.",
      call. = FALSE
    )
  }
  
  do.call(abind::abind, c(arrays, along = 3))
}

write_tps <- function(arr, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  geomorph::writeland.tps(A = arr, file = path, scale = NULL, specID = TRUE)
}

count_missing_landmarks <- function(arr) {
  apply(arr, 3, function(x) sum(!stats::complete.cases(x)))
}

write_missingness_table <- function(arr, path) {
  utils::write.csv(
    data.frame(
      specimen_id = dimnames(arr)[[3]],
      n_missing_landmarks = count_missing_landmarks(arr),
      stringsAsFactors = FALSE
    ),
    path,
    row.names = FALSE
  )
}

read_bilateral_pairs <- function(path, landmark_col, pair_col) {
  pairs_df <- read_table_file(path)
  
  if (!all(c(landmark_col, pair_col) %in% names(pairs_df))) {
    stop(
      "Bilateral pair file must contain columns: ",
      landmark_col, " and ", pair_col,
      call. = FALSE
    )
  }
  
  pairs_df <- pairs_df[, c(landmark_col, pair_col)]
  names(pairs_df) <- c("landmark", "pair")
  pairs_df <- pairs_df[!is.na(pairs_df$pair), ]
  pairs_df$landmark <- as.integer(pairs_df$landmark)
  pairs_df$pair <- as.integer(pairs_df$pair)
  
  pair_id <- vapply(
    seq_len(nrow(pairs_df)),
    function(i) paste(sort(c(pairs_df$landmark[i], pairs_df$pair[i])), collapse = "_"),
    character(1)
  )
  
  pairs_df <- pairs_df[!duplicated(pair_id), ]
  as.matrix(pairs_df[, c("landmark", "pair")])
}

estimate_missing_pca <- function(arr) {
  n_specimens <- dim(arr)[3]
  reshaped <- array(arr, dim = c(dim(arr)[1] * dim(arr)[2], n_specimens))
  reshaped <- t(reshaped)
  
  nb_comp <- missMDA::estim_ncpPCA(reshaped, method.cv = "Kfold", nbsim = 10)
  imputed <- missMDA::imputePCA(reshaped, ncp = nb_comp$ncp)
  
  estimated <- array(
    data = t(imputed$completeObs),
    dim = dim(arr),
    dimnames = dimnames(arr)
  )
  
  list(array = estimated, ncp = nb_comp$ncp)
}

estimate_missing_array <- function(arr, method = c("pca", "TPS")) {
  method <- match.arg(method)
  
  if (!any(is.na(arr))) {
    return(list(array = arr, method = method, ncp = NA_integer_))
  }
  
  if (method == "pca") {
    out <- estimate_missing_pca(arr)
    return(list(array = out$array, method = method, ncp = out$ncp))
  }
  
  if (method == "TPS") {
    estimated <- geomorph::estimate.missing(arr, method = "TPS")
    return(list(array = estimated, method = method, ncp = NA_integer_))
  }
}

fill_unilateral_missing_with_bilat_symmetry <- function(landmarks,
                                                        landmark_pairs,
                                                        imputation_method = "TPS",
                                                        iter = 1,
                                                        seed = 123) {
  # Uses geomorph::bilat.symmetry() for object symmetry and fills only landmarks where one member of a bilateral pair is missing and the other member is observed. Cases where both paired landmarks are missing are left as missing.
  # Observed landmarks are not replaced.
  
  if (is.null(dimnames(landmarks)[[3]])) {
    dimnames(landmarks)[[3]] <- paste0("specimen_", seq_len(dim(landmarks)[3]))
  }
  
  missing_before <- count_missing_landmarks(landmarks)
  estimated <- estimate_missing_array(landmarks, method = imputation_method)
  specimen_ids <- dimnames(landmarks)[[3]]
  
  bilat_result <- geomorph::bilat.symmetry(
    A = estimated$array,
    ind = specimen_ids,
    object.sym = TRUE,
    land.pairs = landmark_pairs,
    iter = iter,
    seed = seed,
    RRPP = TRUE,
    print.progress = FALSE
  )
  
  if (is.null(bilat_result$symm.shape)) {
    stop(
      "geomorph::bilat.symmetry() did not return `symm.shape`. ",
      "Check your geomorph version and bilateral pair file.",
      call. = FALSE
    )
  }
  
  symmetric_component <- bilat_result$symm.shape
  dimnames(symmetric_component) <- dimnames(landmarks)
  
  filled <- landmarks
  fill_mask <- array(FALSE, dim = dim(landmarks), dimnames = dimnames(landmarks))
  
  for (i in seq_len(dim(landmarks)[3])) {
    shape <- landmarks[, , i]
    
    for (j in seq_len(nrow(landmark_pairs))) {
      a <- landmark_pairs[j, 1]
      b <- landmark_pairs[j, 2]
      
      a_missing <- any(is.na(shape[a, ]))
      b_missing <- any(is.na(shape[b, ]))
      
      if (a_missing && !b_missing) {
        filled[a, , i] <- symmetric_component[a, , i]
        fill_mask[a, , i] <- TRUE
      } else if (!a_missing && b_missing) {
        filled[b, , i] <- symmetric_component[b, , i]
        fill_mask[b, , i] <- TRUE
      }
    }
  }
  
  missing_after <- count_missing_landmarks(filled)
  
  list(
    array = filled,
    fill_mask = fill_mask,
    missing_before = missing_before,
    missing_after = missing_after,
    imputation_method = estimated$method,
    imputation_ncp = estimated$ncp,
    bilat_result = bilat_result
  )
}

align_and_restore_missing <- function(arr, imputation_method = "pca") {
  missing_mask <- is.na(arr)
  estimated <- estimate_missing_array(arr, method = imputation_method)
  gpa <- geomorph::gpagen(estimated$array, PrinAxes = FALSE, print.progress = FALSE)
  
  restored <- arr
  restored[!missing_mask] <- gpa$coords[!missing_mask]
  
  list(
    aligned = restored,
    estimated_complete = estimated$array,
    gpa = gpa,
    imputation_method = estimated$method,
    ncp = estimated$ncp
  )
}

make_group_means <- function(aligned_arr,
                             group_map_file,
                             group_mean_imputation_method = "TPS",
                             output_group_dir = NULL) {
  group_map <- utils::read.csv(group_map_file, stringsAsFactors = FALSE, check.names = FALSE)
  
  if (!all(c("group", "specimen_id") %in% names(group_map))) {
    stop("Group map must contain columns named `group` and `specimen_id`.", call. = FALSE)
  }
  
  specimen_ids <- dimnames(aligned_arr)[[3]]
  groups <- split(group_map$specimen_id, group_map$group)
  
  missing_in_data <- setdiff(unique(group_map$specimen_id), specimen_ids)
  if (length(missing_in_data) > 0) {
    warning(
      "The following specimen IDs are in the group map but not in the TPS data:
",
      paste(missing_in_data, collapse = ", ")
    )
  }
  
  unassigned <- setdiff(specimen_ids, unique(group_map$specimen_id))
  if (length(unassigned) > 0) {
    warning(
      "The following TPS specimens are not assigned to any group and will be omitted from group means:
",
      paste(unassigned, collapse = ", ")
    )
  }
  
  estimated <- estimate_missing_array(aligned_arr, method = group_mean_imputation_method)
  complete_arr <- estimated$array
  
  group_names <- names(groups)
  p <- dim(aligned_arr)[1]
  k <- dim(aligned_arr)[2]
  
  group_array <- array(
    NA_real_,
    dim = c(p, k, length(group_names)),
    dimnames = list(
      dimnames(aligned_arr)[[1]],
      dimnames(aligned_arr)[[2]],
      group_names
    )
  )
  
  group_summary <- data.frame()
  
  for (g in group_names) {
    ids <- groups[[g]]
    idx <- which(specimen_ids %in% ids)
    
    if (length(idx) == 0) {
      warning("Skipping group with no matching specimens: ", g)
      next
    }
    
    group_complete <- complete_arr[, , idx, drop = FALSE]
    group_original <- aligned_arr[, , idx, drop = FALSE]
    
    if (length(idx) == 1) {
      mean_shape <- group_complete[, , 1]
    } else {
      mean_shape <- geomorph::mshape(group_complete)
    }
    
    group_array[, , g] <- mean_shape
    
    # Restore a landmark to NA when it is missing in all specimens in the group.
    landmark_missing_all <- apply(is.na(group_original), 1, all)
    group_array[landmark_missing_all, , g] <- NA
    
    if (!is.null(output_group_dir)) {
      one_group <- array(
        group_array[, , g],
        dim = c(p, k, 1),
        dimnames = list(
          dimnames(aligned_arr)[[1]],
          dimnames(aligned_arr)[[2]],
          g
        )
      )
      
      write_tps(one_group, file.path(output_group_dir, paste0(g, ".tps")))
    }
    
    group_summary <- rbind(
      group_summary,
      data.frame(
        group = g,
        n_specimens = length(idx),
        specimen_ids = paste(specimen_ids[idx], collapse = ";"),
        n_group_missing_landmarks = sum(landmark_missing_all),
        stringsAsFactors = FALSE
      )
    )
  }
  
  list(array = group_array, summary = group_summary)
}

clean_names <- function(arr, pattern, replacement) {
  names_old <- dimnames(arr)[[3]]
  dimnames(arr)[[3]] <- gsub(pattern, replacement, names_old)
  arr
}

# -----------------------------------------------------------------------------
# 4. Read and merge TPS files
# -----------------------------------------------------------------------------

message_step("Reading and merging TPS files")

tps_paths <- get_tps_paths(input_tps_dir, input_tps_files)
combined_unaligned <- merge_tps_arrays(paths = tps_paths)

combined_unaligned_path <- file.path(output_dir, "intermediate", "01_combined_unaligned.tps")
write_tps(combined_unaligned, combined_unaligned_path)
write_missingness_table(
  combined_unaligned,
  file.path(output_dir, "tables", "01_combined_unaligned_missingness.csv")
)

message("Combined specimen-level array dimensions: ", paste(dim(combined_unaligned), collapse = " × "))

# -----------------------------------------------------------------------------
# 5. Optional bilat.symmetry fill for one-sided missing bilateral landmarks
# -----------------------------------------------------------------------------

message_step("Optional bilat.symmetry fill")

if (!is.null(bilateral_pairs_file)) {
  landmark_pairs <- read_bilateral_pairs(
    bilateral_pairs_file,
    landmark_col = bilateral_landmark_col,
    pair_col = bilateral_pair_col
  )
  
  bilateral_result <- fill_unilateral_missing_with_bilat_symmetry(
    landmarks = combined_unaligned,
    landmark_pairs = landmark_pairs,
    imputation_method = bilateral_imputation_method,
    iter = bilateral_symmetry_iter,
    seed = 123
  )
  
  combined_bilat <- bilateral_result$array
  
  bilateral_summary <- data.frame(
    specimen_id = dimnames(combined_unaligned)[[3]],
    missing_before = bilateral_result$missing_before,
    missing_after = bilateral_result$missing_after,
    landmarks_filled = bilateral_result$missing_before - bilateral_result$missing_after,
    bilateral_imputation_method = bilateral_result$imputation_method,
    bilateral_imputation_ncp = bilateral_result$imputation_ncp,
    bilat_symmetry_iter = bilateral_symmetry_iter,
    stringsAsFactors = FALSE
  )
  
  utils::write.csv(
    bilateral_summary,
    file.path(output_dir, "tables", "02_bilat_symmetry_fill_summary.csv"),
    row.names = FALSE
  )
  
  combined_bilat_path <- file.path(output_dir, "intermediate", "02_combined_after_bilat_symmetry.tps")
  write_tps(combined_bilat, combined_bilat_path)
  write_missingness_table(
    combined_bilat,
    file.path(output_dir, "tables", "02_combined_after_bilat_symmetry_missingness.csv")
  )
} else {
  message("No bilateral pair file supplied; skipping bilat.symmetry fill.")
  combined_bilat <- combined_unaligned
  combined_bilat_path <- combined_unaligned_path
  
  utils::write.csv(
    data.frame(
      specimen_id = dimnames(combined_unaligned)[[3]],
      missing_before = count_missing_landmarks(combined_unaligned),
      missing_after = count_missing_landmarks(combined_unaligned),
      landmarks_filled = 0,
      bilateral_imputation_method = NA_character_,
      bilateral_imputation_ncp = NA_integer_,
      bilat_symmetry_iter = NA_integer_,
      stringsAsFactors = FALSE
    ),
    file.path(output_dir, "tables", "02_bilat_symmetry_fill_summary.csv"),
    row.names = FALSE
  )
}

# -----------------------------------------------------------------------------
# 6. Estimate remaining missing landmarks, GPA-align, and restore missing mask
# -----------------------------------------------------------------------------

message_step("Estimating remaining missing data, aligning, and restoring missing mask")

alignment <- align_and_restore_missing(
  combined_bilat,
  imputation_method = alignment_imputation_method
)

combined_aligned <- alignment$aligned
combined_aligned_path <- file.path(output_dir, "intermediate", "03_combined_aligned_restored_missing.tps")
write_tps(combined_aligned, combined_aligned_path)
write_missingness_table(
  combined_aligned,
  file.path(output_dir, "tables", "03_combined_aligned_missingness.csv")
)

utils::write.csv(
  data.frame(
    step = "specimen_level_alignment",
    imputation_method = alignment$imputation_method,
    ncp = alignment$ncp,
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "tables", "03_alignment_imputation_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 7. Optional group/taxon mean creation
# -----------------------------------------------------------------------------

message_step("Optional group/taxon mean creation")

if (!is.null(group_map_file)) {
  if (!file.exists(group_map_file)) {
    stop("Group map file not found: ", group_map_file, call. = FALSE)
  }
  
  group_result <- make_group_means(
    aligned_arr = combined_aligned,
    group_map_file = group_map_file,
    group_mean_imputation_method = group_mean_imputation_method,
    output_group_dir = file.path(output_dir, "group_means")
  )
  
  group_means <- group_result$array
  
  if (clean_final_names) {
    group_means <- clean_names(group_means, final_name_pattern, final_name_replacement)
  }
  
  group_means_path <- file.path(output_dir, "intermediate", "04_group_means_restored_missing.tps")
  write_tps(group_means, group_means_path)
  
  utils::write.csv(
    group_result$summary,
    file.path(output_dir, "tables", "04_group_mean_summary.csv"),
    row.names = FALSE
  )
  
  write_missingness_table(
    group_means,
    file.path(output_dir, "tables", "04_group_means_missingness.csv")
  )
  
  final_input <- group_means
} else {
  message("No group map supplied; skipping group means and using specimen-level aligned data as final input.")
  group_means_path <- NA_character_
  final_input <- combined_aligned
}

# -----------------------------------------------------------------------------
# 8. Final alignment and ready-for-analysis TPS
# -----------------------------------------------------------------------------

message_step("Writing final ready-for-analysis TPS")

if (!is.null(group_map_file) && realign_group_means) {
  missing_mask_final <- is.na(final_input)
  
  if (any(missing_mask_final)) {
    estimated_final <- geomorph::estimate.missing(final_input, method = "TPS")
  } else {
    estimated_final <- final_input
  }
  
  final_gpa <- geomorph::gpagen(estimated_final, print.progress = FALSE)
  final_array <- final_gpa$coords
  final_array[missing_mask_final] <- NA
} else {
  final_array <- final_input
}

if (clean_final_names) {
  final_array <- clean_names(final_array, final_name_pattern, final_name_replacement)
}

timestamp <- format(Sys.time(), "%m%d%Y-%H%M%S")
final_tps_path <- file.path(
  output_dir,
  "ready_for_analysis",
  paste0("ready_for_analysis_", timestamp, ".tps")
)

write_tps(final_array, final_tps_path)
write_missingness_table(
  final_array,
  file.path(output_dir, "tables", "05_final_ready_missingness.csv")
)

# Also write a stable filename that is overwritten on rerun.
final_tps_latest_path <- file.path(output_dir, "ready_for_analysis", "ready_for_analysis_latest.tps")
write_tps(final_array, final_tps_latest_path)

# -----------------------------------------------------------------------------
# 9. Pipeline summary
# -----------------------------------------------------------------------------

message_step("Writing pipeline summary")

pipeline_summary <- data.frame(
  item = c(
    "n_input_tps_files",
    "combined_unaligned_path",
    "combined_after_bilat_symmetry_path",
    "combined_aligned_path",
    "group_means_path",
    "final_ready_tps_path",
    "final_ready_tps_latest_path",
    "alignment_imputation_method",
    "alignment_pca_ncp",
    "group_map_file",
    "bilateral_pairs_file"
  ),
  value = c(
    length(tps_paths),
    combined_unaligned_path,
    combined_bilat_path,
    combined_aligned_path,
    group_means_path,
    final_tps_path,
    final_tps_latest_path,
    alignment$imputation_method,
    alignment$ncp,
    ifelse(is.null(group_map_file), NA_character_, group_map_file),
    ifelse(is.null(bilateral_pairs_file), NA_character_, bilateral_pairs_file)
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  pipeline_summary,
  file.path(output_dir, "tables", "00_pipeline_summary.csv"),
  row.names = FALSE
)

utils::capture.output(
  utils::sessionInfo(),
  file = file.path(output_dir, "sessionInfo_gm_tps_processing_pipeline.txt")
)

message("
Done.")
message("Final TPS file: ", final_tps_path)
message("Latest final TPS file: ", final_tps_latest_path)
