#!/usr/bin/env Rscript

find_repo_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "development-code", "GRIN3D_mutation_hotspots.R"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) stop("Could not find the repository root")
    path <- parent
  }
}

parse_integer_grid <- function(value, fallback) {
  if (!nzchar(value)) return(as.integer(fallback))
  result <- suppressWarnings(as.integer(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (anyNA(result) || any(result < 100L)) stop("Trial values must be integers >= 100")
  sort(unique(result))
}

parse_numeric_grid <- function(value, fallback) {
  if (!nzchar(value)) return(as.numeric(fallback))
  result <- suppressWarnings(as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (anyNA(result) || any(!is.finite(result))) stop("pLDDT values must be finite numbers")
  sort(unique(result))
}

repo.root <- find_repo_root()
analysis.root <- file.path(repo.root, "analysis", "mutation_plddt_trials")
run.root <- file.path(analysis.root, "results", "runs")
dir.create(run.root, recursive = TRUE, showWarnings = FALSE)

source(file.path(repo.root, "development-code", "GRIN3D_mutation_hotspots.R"))
source(file.path(repo.root, "development-code", "GRIN3D_mutation_stability.R"))

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args)) args[[1L]] else "quick"
profiles <- list(
  quick = list(trials = c(100L, 200L), plddt = c(0, 70, 90)),
  dense = list(trials = seq.int(100L, 1000L, by = 100L), plddt = c(0, 70, 90)),
  standard = list(trials = c(200L, 500L, 1000L), plddt = c(0, 50, 70, 90)),
  custom = list(trials = c(100L, 200L), plddt = c(0, 70, 90))
)
if (!profile %in% names(profiles)) stop("Profile must be quick, dense, standard, or custom")
trials <- parse_integer_grid(Sys.getenv("GRIN3D_TRIALS"), profiles[[profile]]$trials)
plddt <- parse_numeric_grid(Sys.getenv("GRIN3D_PLDDT"), profiles[[profile]]$plddt)

inputs <- data.frame(
  protein = c("TP53", "PTEN", "SUZ12", "EZH2", "LEF1"),
  lesions = file.path(repo.root, "examples", c("TP53_mutations", "PTEN_mutations", "SUZ12_mutations", "EZH2_mutations", "LEF1_mutations"),
                      "input_files", c("tp53_start_end.csv", "pten_start_end.csv", "suz12_start_end.csv", "ezh2_start_end.csv", "lef1_start_end.csv")),
  coordinates = file.path(repo.root, "examples", c("TP53_mutations", "PTEN_mutations", "SUZ12_mutations", "EZH2_mutations", "LEF1_mutations"),
                          "input_files", c("tp53_alphafold_v6_coordinates.csv", "pten_coordinates.csv", "suz12_coordinates.csv", "ezh2_coordinates.csv", "lef1_coordinates.csv")),
  stringsAsFactors = FALSE
)
if (any(!file.exists(unlist(inputs[c("lesions", "coordinates")])))) stop("One or more input files are missing")
inputs$has_plddt <- vapply(inputs$coordinates, function(f) "plddt" %in% names(read.csv(f, nrows = 1L)), logical(1L))

grid.rows <- lapply(seq_len(nrow(inputs)), function(i) {
  thresholds <- if (inputs$has_plddt[[i]]) plddt else NA_real_
  expand.grid(
    protein = inputs$protein[[i]], n_sim = trials,
    plddt_threshold = thresholds, stringsAsFactors = FALSE
  )
})
config <- do.call(rbind, grid.rows)
config$has_plddt <- inputs$has_plddt[match(config$protein, inputs$protein)]
config$lesion_file <- inputs$lesions[match(config$protein, inputs$protein)]
config$coordinate_file <- inputs$coordinates[match(config$protein, inputs$protein)]
threshold.label <- ifelse(is.na(config$plddt_threshold), "none", format(config$plddt_threshold, trim = TRUE))
config$label <- paste0(config$protein, "_trials", config$n_sim, "_plddt", threshold.label)
config$random_seed <- 20260916L
config$min_subjects <- 2L
config$min_events <- 2L
config$calibration_fraction <- 0.5
config$alpha <- 0.05
config$profile <- profile
config$status <- ifelse(
  file.exists(file.path(file.path(run.root, config$label), "GRIN3D_mutation_hotspot_results.rds")),
  "available", "pending"
)
config$results_dir <- file.path(run.root, config$label)

protein.filter <- trimws(strsplit(Sys.getenv("GRIN3D_PROTEINS"), ",", fixed = TRUE)[[1L]])
protein.filter <- protein.filter[nzchar(protein.filter)]
worker.mode <- length(protein.filter) > 0L
if (worker.mode && any(!protein.filter %in% inputs$protein)) stop("Unknown GRIN3D_PROTEINS value")
if (!worker.mode) write.csv(config, file.path(analysis.root, "configurations.csv"), row.names = FALSE)
if (identical(tolower(Sys.getenv("GRIN3D_PLAN_ONLY")), "true")) {
  message("Configuration plan written; no analyses were run.")
  quit(save = "no", status = 0L)
}

run.indices <- if (worker.mode) which(config$protein %in% protein.filter) else seq_len(nrow(config))
for (i in run.indices) {
  row <- config[i, ]
  result.file <- file.path(row$results_dir, "GRIN3D_mutation_hotspot_results.rds")
  if (file.exists(result.file)) {
    message("Reusing completed run: ", row$label)
    config$status[[i]] <- "reused"
    next
  }
  input.row <- inputs[inputs$protein == row$protein, ]
  coordinate.columns <- list(
    residue = "residue", x = "x", y = "y", z = "z",
    confidence = if (row$has_plddt) "plddt" else NULL,
    opportunity = NULL
  )
  message("Starting ", i, "/", nrow(config), ": ", row$label)
  grin3d_run_stability_configuration(
    label = row$label,
    protein = row$protein,
    lesion.file = input.row$lesions,
    coordinate.file = input.row$coordinates,
    output.root = run.root,
    n.sim = row$n_sim,
    min.subjects = row$min_subjects,
    min.events = row$min_events,
    random.seed = row$random_seed,
    lesion.columns = list(subject = "ID", event = "id", start = "start", end = "end"),
    coordinate.columns = coordinate.columns,
    min.confidence = if (row$has_plddt) row$plddt_threshold else 0,
    calibration.fraction = row$calibration_fraction,
    alpha = row$alpha,
    progress.every = 100L
  )
  config$status[[i]] <- "completed"
  if (!worker.mode) write.csv(config, file.path(analysis.root, "configurations.csv"), row.names = FALSE)
}

if (!worker.mode) write.csv(config, file.path(analysis.root, "configurations.csv"), row.names = FALSE)
message("Sensitivity grid complete. Run 02_summarize_visualize.R next.")
