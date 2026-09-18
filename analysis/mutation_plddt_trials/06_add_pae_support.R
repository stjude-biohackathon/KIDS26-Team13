#!/usr/bin/env Rscript

# Add AlphaFold predicted-aligned-error (PAE) support to completed hotspot runs.
# This is a post-calibration structural-confidence layer: it does not rerun the
# null simulations and does not change cluster membership or statistical p-values.

find_repo_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "development-code", "GRIN3D_mutation_hotspots.R"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Could not find the repository root")
    path <- parent
  }
}

parse_residues <- function(x) {
  values <- suppressWarnings(as.integer(strsplit(as.character(x), ";", fixed = TRUE)[[1L]]))
  sort(unique(values[!is.na(values)]))
}

repo.root <- find_repo_root()
analysis.root <- file.path(repo.root, "analysis", "mutation_plddt_trials")
summary.dir <- file.path(analysis.root, "results", "summary")
figure.dir <- file.path(analysis.root, "results", "figures")
pae.dir <- file.path(analysis.root, "input_files", "pae")
dir.create(summary.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pae.dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(repo.root, "development-code", "GRIN3D_mutation_hotspots.R"))

proteins <- data.frame(
  protein = c("TP53", "PTEN", "SUZ12", "EZH2", "LEF1"),
  accession = c("P04637", "P60484", "Q15022", "Q15910", "Q9UJU2"),
  model_version = 6L,
  stringsAsFactors = FALSE
)
proteins$url <- sprintf(
  "https://alphafold.ebi.ac.uk/files/AF-%s-F1-predicted_aligned_error_v%d.json",
  proteins$accession, proteins$model_version
)
proteins$pae_file <- file.path(
  pae.dir,
  sprintf("AF-%s-F1-predicted_aligned_error_v%d.json", proteins$accession, proteins$model_version)
)

for (i in seq_len(nrow(proteins))) {
  if (!file.exists(proteins$pae_file[[i]])) {
    message("Downloading PAE for ", proteins$protein[[i]], " from AlphaFold DB")
    tryCatch(
      utils::download.file(proteins$url[[i]], proteins$pae_file[[i]], mode = "wb", quiet = TRUE),
      error = function(error) stop(
        "Could not download ", proteins$url[[i]], ". ",
        "Download it manually to ", proteins$pae_file[[i]], ". ", conditionMessage(error)
      )
    )
  }
}
proteins$downloaded_on <- format(Sys.Date(), "%Y-%m-%d")
proteins$source <- "AlphaFold Protein Structure Database"
write.csv(proteins, file.path(pae.dir, "pae_sources.csv"), row.names = FALSE)

pae.matrices <- setNames(lapply(seq_len(nrow(proteins)), function(i) {
  read_alphafold_pae(proteins$pae_file[[i]])
}), proteins$protein)

config.file <- file.path(analysis.root, "configurations.csv")
if (!file.exists(config.file)) stop("Missing configurations.csv; run 01_run_sensitivity.R first")
config <- read.csv(config.file, stringsAsFactors = FALSE, check.names = FALSE)

pae.columns <- c(
  "pae_n_pairs", "pae_median", "pae_p90", "pae_max",
  "pae_reliable_pair_fraction", "pae_supported", "pae_interpretation",
  "pae_cutoff", "pae_minimum_reliable_fraction"
)
all.clusters <- list()
run.summaries <- list()
cutoff.rows <- list()
cluster.index <- summary.index <- cutoff.index <- 0L

for (i in seq_len(nrow(config))) {
  cfg <- config[i, , drop = FALSE]
  result.file <- file.path(cfg$results_dir[[1L]], "GRIN3D_mutation_hotspot_results.rds")
  if (!file.exists(result.file)) next
  result <- readRDS(result.file)
  pae <- pae.matrices[[cfg$protein[[1L]]]]
  if (is.null(pae)) stop("No PAE matrix configured for ", cfg$protein[[1L]])
  if (any(result$coordinates$residue > nrow(pae))) {
    stop("PAE matrix is shorter than the coordinate numbering for ", cfg$label[[1L]])
  }

  base.clusters <- result$clusters[, setdiff(names(result$clusters), pae.columns), drop = FALSE]
  result$clusters <- annotate_clusters_with_pae(
    base.clusters, pae, cutoff = 10, minimum.reliable.fraction = 0.80
  )
  result$pae <- pae
  result$settings$pae.file <- proteins$pae_file[proteins$protein == cfg$protein[[1L]]][[1L]]
  result$settings$analysis.options$pae.cutoff <- 10
  result$settings$analysis.options$pae.minimum.reliable.fraction <- 0.80
  saveRDS(result, result.file)
  write.csv(
    result$clusters,
    file.path(cfg$results_dir[[1L]], "mutation_hotspot_clusters.csv"),
    row.names = FALSE
  )

  annotated <- result$clusters
  annotated$protein <- cfg$protein[[1L]]
  annotated$run_label <- cfg$label[[1L]]
  annotated$n_sim <- cfg$n_sim[[1L]]
  annotated$plddt_threshold <- cfg$plddt_threshold[[1L]]
  cluster.index <- cluster.index + 1L
  all.clusters[[cluster.index]] <- annotated

  multi <- annotated$pae_n_pairs > 0L
  summary.index <- summary.index + 1L
  run.summaries[[summary.index]] <- data.frame(
    protein = cfg$protein[[1L]], run_label = cfg$label[[1L]],
    n_sim = cfg$n_sim[[1L]], plddt_threshold = cfg$plddt_threshold[[1L]],
    n_clusters = nrow(annotated), n_multiresidue_clusters = sum(multi),
    n_pae_supported = sum(annotated$pae_supported %in% TRUE),
    pae_supported_fraction = if (any(multi)) mean(annotated$pae_supported[multi]) else NA_real_,
    n_significant = sum(annotated$significant_any),
    n_significant_pae_supported = sum(annotated$significant_any & annotated$pae_supported %in% TRUE),
    stringsAsFactors = FALSE
  )

  for (cutoff in c(5, 10, 15)) {
    sensitivity <- annotate_clusters_with_pae(
      base.clusters, pae, cutoff = cutoff, minimum.reliable.fraction = 0.80
    )
    cutoff.index <- cutoff.index + 1L
    cutoff.rows[[cutoff.index]] <- data.frame(
      protein = cfg$protein[[1L]], run_label = cfg$label[[1L]],
      n_sim = cfg$n_sim[[1L]], plddt_threshold = cfg$plddt_threshold[[1L]],
      cluster_id = sensitivity$cluster_id, residues = sensitivity$residues,
      significant_any = sensitivity$significant_any,
      pae_cutoff = cutoff,
      pae_reliable_pair_fraction = sensitivity$pae_reliable_pair_fraction,
      pae_supported = sensitivity$pae_supported,
      stringsAsFactors = FALSE
    )
  }
}

cluster.table <- do.call(rbind, all.clusters)
run.summary <- do.call(rbind, run.summaries)
cutoff.table <- do.call(rbind, cutoff.rows)
write.csv(cluster.table, file.path(summary.dir, "cluster_pae_support.csv"), row.names = FALSE)
write.csv(run.summary, file.path(summary.dir, "pae_run_summary.csv"), row.names = FALSE)
write.csv(cutoff.table, file.path(summary.dir, "pae_cutoff_sensitivity.csv"), row.names = FALSE)

# Plot one full-protein PAE heatmap per protein. Rectangles show up to ten
# lowest-p-value multi-residue clusters from the unfiltered 1,000-trial run.
for (protein in proteins$protein) {
  pae <- pae.matrices[[protein]]
  candidates <- cluster.table[
    cluster.table$protein == protein & cluster.table$n_sim == max(config$n_sim) &
      cluster.table$plddt_threshold == min(config$plddt_threshold) &
      cluster.table$pae_n_pairs > 0L,
    , drop = FALSE
  ]
  candidates <- head(candidates[order(candidates$p_any_joint), , drop = FALSE], 10L)
  png(file.path(figure.dir, paste0(protein, "_pae_cluster_heatmap.png")),
      width = 1800, height = 1600, res = 180)
  old <- par(mar = c(5, 5, 4, 7), xpd = FALSE)
  palette <- hcl.colors(100, "YlOrRd")
  image(
    seq_len(nrow(pae)), seq_len(ncol(pae)), pae,
    zlim = c(0, max(30, max(pae))), col = palette, useRaster = TRUE,
    xlab = "Aligned residue", ylab = "Predicted residue",
    main = paste(protein, "AlphaFold PAE with candidate clusters")
  )
  if (nrow(candidates)) {
    colors <- grDevices::rainbow(nrow(candidates), s = 0.75, v = 0.70)
    for (j in seq_len(nrow(candidates))) {
      residues <- parse_residues(candidates$residues[[j]])
      rect(min(residues), min(residues), max(residues), max(residues),
           border = colors[[j]], lwd = 2)
    }
    legend("right", inset = c(-0.23, 0), xpd = NA, bty = "n", cex = 0.72,
           legend = paste(candidates$cluster_id, paste0("p=", signif(candidates$p_any_joint, 3))),
           col = colors, lwd = 2, title = "Lowest-p clusters")
  }
  par(old)
  dev.off()
}

message("PAE annotations added to completed runs and summary tables.")
message("PAE figures written to: ", figure.dir)
