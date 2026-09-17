#!/usr/bin/env Rscript

# ==============================================================================
# Validation Tests for PTEN CNA Linkage Implementation
# ==============================================================================

cat("==============================================================================\n")
cat("Running Linkage Implementation Validation Tests for PTEN CNA\n")
cat("==============================================================================\n\n")

source("GRIN3D_exon_CNA_hotspots.R")

test_pass_count <- 0L
test_fail_count <- 0L

assert_true <- function(cond, message) {
  if (isTRUE(cond)) {
    cat(sprintf("[PASS] %s\n", message))
    test_pass_count <<- test_pass_count + 1L
  } else {
    cat(sprintf("[FAIL] %s\n", message))
    test_fail_count <<- test_fail_count + 1L
  }
}

assert_error <- function(expr, pattern = NULL, message) {
  has_pattern <- !is.null(pattern)
  res <- tryCatch(
    {
      eval(substitute(expr), envir = parent.frame())
      FALSE
    },
    error = function(e) {
      if (!has_pattern || grepl(pattern, e$message, ignore.case = TRUE)) {
        TRUE
      } else {
        cat(sprintf("   Unexpected error message: %s\n", e$message))
        FALSE
      }
    }
  )
  if (res) {
    cat(sprintf("[PASS] %s\n", message))
    test_pass_count <<- test_pass_count + 1L
  } else {
    cat(sprintf("[FAIL] %s\n", message))
    test_fail_count <<- test_fail_count + 1L
  }
}

# ------------------------------------------------------------------------------
# 1. Default linkage.method = "complete" in all function signatures
# ------------------------------------------------------------------------------
cat("\n--- Test Suite 1: Default Function Signatures ---\n")
assert_true(
  identical(formals(build_cna_tree_sets)$linkage.method, "complete"),
  "build_cna_tree_sets defaults to linkage.method = 'complete'"
)
assert_true(
  identical(formals(build_cna_candidate_clusters)$linkage.method, "complete"),
  "build_cna_candidate_clusters defaults to linkage.method = 'complete'"
)
assert_true(
  identical(formals(run_cna_null_simulations)$linkage.method, "complete"),
  "run_cna_null_simulations defaults to linkage.method = 'complete'"
)
assert_true(
  identical(formals(run_grin3d_exon_cna_hotspots)$linkage.method, "complete"),
  "run_grin3d_exon_cna_hotspots defaults to linkage.method = 'complete'"
)

# ------------------------------------------------------------------------------
# 2. Valid and Invalid Linkage Methods
# ------------------------------------------------------------------------------
cat("\n--- Test Suite 2: Linkage Method Validation ---\n")
valid_methods <- c("complete", "average", "single", "ward.D2")
for (m in valid_methods) {
  assert_true(
    tryCatch({ validate_linkage_method(m); TRUE }, error = function(e) FALSE),
    sprintf("Valid method '%s' is accepted by validate_linkage_method", m)
  )
}

assert_error(
  validate_linkage_method("centroid"),
  "Invalid linkage.method",
  "Invalid method 'centroid' throws clear error"
)
assert_error(
  validate_linkage_method("median"),
  "Invalid linkage.method",
  "Invalid method 'median' throws clear error"
)
assert_error(
  validate_linkage_method(123),
  "Invalid linkage.method",
  "Non-character linkage.method throws clear error"
)
assert_error(
  validate_linkage_method(NULL),
  "Invalid linkage.method",
  "NULL linkage.method throws clear error"
)
assert_error(
  validate_linkage_method(c("complete", "single")),
  "Invalid linkage.method",
  "Vector of methods throws clear error"
)

# ------------------------------------------------------------------------------
# 3. Duplicate Event IDs, Exon IDs, and Residue Coordinates Fail Clearly
# ------------------------------------------------------------------------------
cat("\n--- Test Suite 3: Input Integrity Checks ---\n")

# Duplicate residue coordinates check
dup_coord_data <- data.frame(
  residue = c(1L, 2L, 2L),
  x = c(0, 1, 2),
  y = c(0, 1, 2),
  z = c(0, 1, 2)
)
assert_error(
  prepare_cna_coordinates(
    dup_coord_data,
    list(residue = "residue", x = "x", y = "y", z = "z", confidence = NULL),
    min.confidence = 70
  ),
  "residue",
  "Duplicate residue coordinates fail with clear error"
)

# Duplicate exon IDs check
exon_base <- utils::read.csv("input_files/PTEN_exon_to_protein_map.csv", stringsAsFactors = FALSE)
dup_exon_data <- rbind(exon_base[1L, ], exon_base)
assert_error(
  prepare_coding_exons(
    dup_exon_data,
    list(
      gene_id = "gene_id", gene_name = "gene_name", transcript_id = "transcript_id",
      protein_id = "protein_id", exon_id = "exon_id", exon_rank = "exon_rank",
      chrom = "chrom", strand = "strand", genomic_start = "cds_genomic_start",
      genomic_end = "cds_genomic_end", residue_start = "aa_start", residue_end = "aa_end"
    ),
    coordinates = utils::read.csv("input_files/PTEN_alphafold_coordinates.csv", stringsAsFactors = FALSE)[, c("residue", "x", "y", "z")]
  ),
  "unique",
  "Duplicate exon IDs fail with clear error"
)

# Duplicate event IDs check
cna_base <- utils::read.csv("input_files/PTEN_CNA_lesions.csv", stringsAsFactors = FALSE)
cna_dup_events <- cna_base
cna_dup_events$event_id <- c("event_1", "event_1", paste0("event_", seq_len(nrow(cna_base) - 2L) + 1L))
assert_error(
  prepare_genomic_cna_events(
    cna_dup_events,
    list(subject = "ID", event = "event_id", chrom = "chrom", start = "loc.start", end = "loc.end", type = "lsn.type")
  ),
  "event",
  "Duplicate event IDs fail with clear error"
)

# ------------------------------------------------------------------------------
# 4. Exon Set Deduplication Across 1D and 3D Trees
# ------------------------------------------------------------------------------
cat("\n--- Test Suite 4: Exon Set Deduplication Across Trees ---\n")

rds <- readRDS("results/GRIN3D_exon_CNA_results.rds")
homdel_map <- rds$event_exon_mapping[rds$event_exon_mapping$cna_type == "HOMDEL", ]
candidates_complete <- build_cna_candidate_clusters(
  homdel_map, rds$exons, rds$exon_distance_matrix,
  min.subjects = 2L, min.events = 2L, linkage.method = "complete"
)

# Verify no duplicate candidate_keys
assert_true(
  !anyDuplicated(candidates_complete$candidate_key),
  "Candidate clusters contain no duplicate exon sets (unique candidate_key)"
)

# Verify 1D;3D combined tree sources
combined_sources <- candidates_complete[candidates_complete$tree_sources == "1D;3D", ]
assert_true(
  nrow(combined_sources) > 0L,
  "Sets detected by both 1D and 3D trees are merged into tree_sources = '1D;3D'"
)

# ------------------------------------------------------------------------------
# 5. Baseline Preservation: Default call matches explicit linkage.method = "complete"
# ------------------------------------------------------------------------------
cat("\n--- Test Suite 5: Baseline Equivalence ---\n")

candidates_default <- build_cna_candidate_clusters(
  homdel_map, rds$exons, rds$exon_distance_matrix,
  min.subjects = 2L, min.events = 2L
)
assert_true(
  identical(candidates_default, candidates_complete),
  "build_cna_candidate_clusters default call is identical to explicit linkage.method = 'complete'"
)

# Tree sets baseline equivalence
items <- sort(unique(homdel_map$exon_order))
trees_default <- build_cna_tree_sets(items, rds$exon_distance_matrix)
trees_complete <- build_cna_tree_sets(items, rds$exon_distance_matrix, linkage.method = "complete")
assert_true(
  identical(trees_default, trees_complete),
  "build_cna_tree_sets default call is identical to explicit linkage.method = 'complete'"
)

cat("\n==============================================================================\n")
cat(sprintf("Validation Results: %d passed, %d failed.\n", test_pass_count, test_fail_count))
cat("==============================================================================\n")

if (test_fail_count > 0L) {
  quit(status = 1L)
} else {
  cat("All validation tests passed successfully!\n")
}
