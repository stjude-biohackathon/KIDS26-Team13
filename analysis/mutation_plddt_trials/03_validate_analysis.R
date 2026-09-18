#!/usr/bin/env Rscript

find_repo_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "development-code", "GRIN3D_mutation_hotspots.R"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Could not find the repository root")
    path <- parent
  }
}

repo.root <- find_repo_root()
analysis.root <- file.path(repo.root, "analysis", "mutation_plddt_trials")
summary.dir <- file.path(analysis.root, "results", "summary")
config <- read.csv(file.path(analysis.root, "configurations.csv"), stringsAsFactors = FALSE)
checks <- list(); n <- 0L
add_check <- function(name, pass, detail) {
  n <<- n + 1L
  checks[[n]] <<- data.frame(check = name, passed = isTRUE(pass), detail = detail, stringsAsFactors = FALSE)
}

result.files <- file.path(config$results_dir, "GRIN3D_mutation_hotspot_results.rds")
manifest.files <- file.path(config$results_dir, "stability_run_manifest.csv")
add_check("all_result_objects_exist", all(file.exists(result.files)),
          paste(sum(file.exists(result.files)), "of", length(result.files), "result objects found"))
add_check("all_manifests_exist", all(file.exists(manifest.files)),
          paste(sum(file.exists(manifest.files)), "of", length(manifest.files), "manifests found"))

runs <- setNames(lapply(result.files, readRDS), config$label)
manifests <- do.call(rbind, lapply(manifest.files, function(f) read.csv(f, stringsAsFactors = FALSE)))
add_check("configured_trial_counts_used", identical(as.integer(manifests$n_sim), as.integer(config$n_sim)),
          "Manifest trial counts match configurations.csv")

fixed.columns <- c("min_subjects", "min_events", "random_seed", "require_complete_event_mapping",
                   "avoid_within_subject_overlap", "opportunity_weights")
for (column in fixed.columns) {
  values <- unique(manifests[[column]])
  add_check(paste0("fixed_", column), length(values) == 1L,
            paste(column, "values:", paste(values, collapse = ";")))
}

group.keys <- ifelse(is.na(config$plddt_threshold), "NA", as.character(config$plddt_threshold))
groups <- split(seq_len(nrow(config)), paste(config$protein, group.keys, sep = "_"))
for (indices in groups) {
  if (length(indices) < 2L) next
  memberships <- lapply(config$label[indices], function(label) sort(runs[[label]]$clusters$candidate_key))
  unchanged <- all(vapply(memberships[-1L], identical, logical(1L), memberships[[1L]]))
  key <- paste(config$protein[indices[[1L]]],
               ifelse(is.na(config$plddt_threshold[indices[[1L]]]), "no_pLDDT",
                      paste0("pLDDT_", config$plddt_threshold[indices[[1L]]])), sep = "_")
  add_check(paste0("trial_membership_invariant_", key), unchanged,
            paste("Compared", length(indices), "trial counts at a fixed structural threshold"))
}

for (protein in unique(config$protein)) {
  thresholds <- sort(unique(config$plddt_threshold[config$protein == protein]))
  add_check(paste0("plddt_available_", protein),
            isTRUE(all.equal(as.numeric(thresholds), c(0, 70, 90))),
            paste(protein, "pLDDT thresholds:", paste(thresholds, collapse = ", ")))
}

for (protein in unique(config$protein)) {
  subset <- config[config$protein == protein & config$n_sim == max(config$n_sim), ]
  subset <- subset[order(subset$plddt_threshold), ]
  eligible <- vapply(subset$label, function(label) nrow(runs[[label]]$coordinates), integer(1L))
  mapped <- vapply(subset$label, function(label) nrow(runs[[label]]$mapped_events), integer(1L))
  add_check(paste0("eligible_residues_nonincreasing_", protein), all(diff(eligible) <= 0),
            paste(eligible, collapse = " -> "))
  add_check(paste0("mapped_events_nonincreasing_", protein), all(diff(mapped) <= 0),
            paste(mapped, collapse = " -> "))
}

validation <- do.call(rbind, checks)
dir.create(summary.dir, recursive = TRUE, showWarnings = FALSE)
write.csv(validation, file.path(summary.dir, "validation_checks.csv"), row.names = FALSE)
print(validation, row.names = FALSE)
if (any(!validation$passed)) stop("One or more sensitivity-analysis validations failed")
message("All validation checks passed.")
