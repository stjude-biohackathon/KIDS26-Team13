# Utilities for evaluating GRIN3D mutation-hotspot stability across runs.
# Source GRIN3D_mutation_hotspots.R before calling the execution helper.

grin3d_stability_parse_residues <- function(x) {
  values <- suppressWarnings(as.integer(strsplit(as.character(x), ";", fixed = TRUE)[[1L]]))
  sort(unique(values[!is.na(values)]))
}

# Convert an AlphaFold PDB into the coordinate schema used by the mutation
# prototype. AlphaFold stores pLDDT in the PDB B-factor field.
grin3d_alphafold_pdb_to_coordinates <- function(pdb.file, output.file) {
  lines <- readLines(pdb.file, warn = FALSE)
  keep <- substr(lines, 1L, 6L) == "ATOM  " &
    trimws(substr(lines, 13L, 16L)) == "CA" &
    substr(lines, 17L, 17L) %in% c(" ", "A")
  atoms <- lines[keep]
  if (!length(atoms)) stop("No C-alpha ATOM records found in ", pdb.file)
  coordinates <- data.frame(
    residue = suppressWarnings(as.integer(trimws(substr(atoms, 23L, 26L)))),
    x = suppressWarnings(as.numeric(trimws(substr(atoms, 31L, 38L)))),
    y = suppressWarnings(as.numeric(trimws(substr(atoms, 39L, 46L)))),
    z = suppressWarnings(as.numeric(trimws(substr(atoms, 47L, 54L)))),
    plddt = suppressWarnings(as.numeric(trimws(substr(atoms, 61L, 66L)))),
    stringsAsFactors = FALSE
  )
  if (anyNA(coordinates) || anyDuplicated(coordinates$residue)) {
    stop("PDB has missing or duplicated C-alpha coordinates")
  }
  coordinates <- coordinates[order(coordinates$residue), , drop = FALSE]
  expected <- seq.int(min(coordinates$residue), max(coordinates$residue))
  if (!identical(coordinates$residue, expected) || expected[[1L]] != 1L) {
    stop("AlphaFold residues are not a complete sequence beginning at position 1")
  }
  dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(coordinates, output.file, row.names = FALSE)
  invisible(coordinates)
}

grin3d_cluster_similarity_matrix <- function(reference, comparison) {
  if (!all(c("cluster_id", "residues") %in% names(reference)) ||
      !all(c("cluster_id", "residues") %in% names(comparison))) {
    stop("Both cluster tables require cluster_id and residues columns")
  }
  reference.sets <- lapply(reference$residues, grin3d_stability_parse_residues)
  comparison.sets <- lapply(comparison$residues, grin3d_stability_parse_residues)
  jaccard <- matrix(
    0, nrow = nrow(reference), ncol = nrow(comparison),
    dimnames = list(reference$cluster_id, comparison$cluster_id)
  )
  containment <- jaccard
  for (i in seq_along(reference.sets)) {
    for (j in seq_along(comparison.sets)) {
      shared <- length(intersect(reference.sets[[i]], comparison.sets[[j]]))
      union.size <- length(union(reference.sets[[i]], comparison.sets[[j]]))
      minimum.size <- min(length(reference.sets[[i]]), length(comparison.sets[[j]]))
      jaccard[i, j] <- if (union.size) shared / union.size else 0
      containment[i, j] <- if (minimum.size) shared / minimum.size else 0
    }
  }
  list(jaccard = jaccard, containment = containment)
}

# Greedy one-to-one matching prioritizes Jaccard, then containment, then the
# smallest difference in cluster size. IDs are never assumed to persist.
grin3d_match_clusters <- function(
    reference,
    comparison,
    reference.significant.only = TRUE,
    comparison.significant.only = FALSE) {
  if (reference.significant.only) {
    reference <- reference[reference$significant_any %in% TRUE, , drop = FALSE]
  }
  if (comparison.significant.only) {
    comparison <- comparison[comparison$significant_any %in% TRUE, , drop = FALSE]
  }
  if (!nrow(reference) || !nrow(comparison)) {
    return(data.frame(
      reference_cluster_id = character(), comparison_cluster_id = character(),
      jaccard = numeric(), containment = numeric(), stringsAsFactors = FALSE
    ))
  }
  similarity <- grin3d_cluster_similarity_matrix(reference, comparison)
  candidates <- expand.grid(
    reference_index = seq_len(nrow(reference)),
    comparison_index = seq_len(nrow(comparison)),
    KEEP.OUT.ATTRS = FALSE
  )
  candidates$jaccard <- similarity$jaccard[cbind(
    candidates$reference_index, candidates$comparison_index
  )]
  candidates$containment <- similarity$containment[cbind(
    candidates$reference_index, candidates$comparison_index
  )]
  candidates$size_difference <- abs(
    reference$n_residues[candidates$reference_index] -
      comparison$n_residues[candidates$comparison_index]
  )
  candidates <- candidates[order(
    -candidates$jaccard, -candidates$containment, candidates$size_difference
  ), , drop = FALSE]
  used.reference <- integer(); used.comparison <- integer(); selected <- integer()
  for (k in seq_len(nrow(candidates))) {
    i <- candidates$reference_index[[k]]
    j <- candidates$comparison_index[[k]]
    if (i %in% used.reference || j %in% used.comparison) next
    selected <- c(selected, k)
    used.reference <- c(used.reference, i)
    used.comparison <- c(used.comparison, j)
  }
  matched <- candidates[selected, , drop = FALSE]
  ref <- reference[matched$reference_index, , drop = FALSE]
  cmp <- comparison[matched$comparison_index, , drop = FALSE]
  data.frame(
    reference_cluster_id = ref$cluster_id,
    comparison_cluster_id = cmp$cluster_id,
    reference_residues = ref$residues,
    comparison_residues = cmp$residues,
    reference_n_residues = ref$n_residues,
    comparison_n_residues = cmp$n_residues,
    reference_tree_sources = ref$tree_sources,
    comparison_tree_sources = cmp$tree_sources,
    reference_significant = ref$significant_any,
    comparison_significant = cmp$significant_any,
    reference_p_any = ref$p_any_joint,
    comparison_p_any = cmp$p_any_joint,
    jaccard = round(matched$jaccard, 6L),
    containment = round(matched$containment, 6L),
    membership_interpretation = ifelse(
      matched$jaccard >= 0.8, "highly_stable",
      ifelse(
        matched$containment >= 0.8, "stable_nested_core",
        ifelse(matched$jaccard >= 0.5, "boundary_sensitive", "unstable")
      )
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

grin3d_plot_jaccard_heatmap <- function(
    reference,
    comparison,
    reference.label = "Reference",
    comparison.label = "Comparison",
    reference.significant.only = TRUE) {
  if (reference.significant.only) {
    reference <- reference[reference$significant_any %in% TRUE, , drop = FALSE]
  }
  if (!nrow(reference) || !nrow(comparison)) {
    plot.new(); title("No clusters available for Jaccard heatmap"); return(invisible(NULL))
  }
  similarity <- grin3d_cluster_similarity_matrix(reference, comparison)$jaccard
  old <- par(mar = c(8, 8, 3, 2))
  on.exit(par(old), add = TRUE)
  image(
    x = seq_len(nrow(similarity)), y = seq_len(ncol(similarity)),
    z = similarity, zlim = c(0, 1), col = hcl.colors(20, "YlOrRd", rev = TRUE),
    axes = FALSE, xlab = reference.label, ylab = comparison.label,
    main = "Cluster residue-set Jaccard similarity"
  )
  axis(1, at = seq_len(nrow(similarity)), labels = rownames(similarity), las = 2, cex.axis = 0.65)
  axis(2, at = seq_len(ncol(similarity)), labels = colnames(similarity), las = 2, cex.axis = 0.55)
  box()
  invisible(similarity)
}

grin3d_plot_cluster_matches <- function(
    matches,
    reference.label = "Reference",
    comparison.label = "Comparison",
    minimum.jaccard = 0) {
  matches <- matches[matches$jaccard >= minimum.jaccard, , drop = FALSE]
  if (!nrow(matches)) {
    plot.new(); title("No matched clusters to plot"); return(invisible(NULL))
  }
  matches <- matches[order(matches$jaccard), , drop = FALSE]
  n <- nrow(matches)
  y.left <- seq_len(n)
  right.order <- order(matches$comparison_cluster_id)
  y.right <- match(seq_len(n), right.order)
  colors <- hcl.colors(100, "YlOrRd", rev = TRUE)[
    pmax(1L, pmin(100L, 1L + floor(matches$jaccard * 99)))
  ]
  old <- par(mar = c(2, 10, 3, 10), xpd = NA)
  on.exit(par(old), add = TRUE)
  plot(c(0, 1), c(0.5, n + 0.5), type = "n", axes = FALSE, xlab = "", ylab = "",
       main = "Best one-to-one cluster matches")
  segments(0, y.left, 1, y.right, col = colors, lwd = 1 + 5 * matches$jaccard)
  points(rep(0, n), y.left, pch = 21, bg = colors)
  points(rep(1, n), y.right, pch = 21, bg = colors)
  text(0, y.left, matches$reference_cluster_id, pos = 2, cex = 0.75)
  text(1, y.right, matches$comparison_cluster_id, pos = 4, cex = 0.75)
  mtext(reference.label, side = 3, at = 0, line = 0.3, font = 2)
  mtext(comparison.label, side = 3, at = 1, line = 0.3, font = 2)
  invisible(matches)
}

grin3d_inputs_available <- function(lesion.file, coordinate.file, protein) {
  missing <- c(lesion.file, coordinate.file)[!file.exists(c(lesion.file, coordinate.file))]
  if (length(missing)) {
    message(protein, " cells are ready but skipped. Missing: ", paste(missing, collapse = "; "))
    return(FALSE)
  }
  TRUE
}

grin3d_run_stability_configuration <- function(
    label,
    protein,
    lesion.file,
    coordinate.file,
    output.root,
    n.sim = 200L,
    min.subjects = 2L,
    min.events = 2L,
    random.seed = 20260828L,
    require.complete.event.mapping = TRUE,
    avoid.within.subject.overlap = TRUE,
    lesion.columns = list(subject = "ID", event = "id", start = "start", end = "end"),
    coordinate.columns = list(
      residue = "residue", x = "x", y = "y", z = "z",
      confidence = NULL, opportunity = NULL
    ),
    min.confidence = 70,
    calibration.fraction = 0.50,
    alpha = 0.05,
    progress.every = 100L) {
  if (!exists("run_grin3d_mutation_hotspots", mode = "function")) {
    stop("Source GRIN3D_mutation_hotspots.R before running a configuration")
  }
  results.dir <- file.path(output.root, label)
  started <- proc.time()[["elapsed"]]
  result <- run_grin3d_mutation_hotspots(
    lesion.file = lesion.file, coordinate.file = coordinate.file,
    results.dir = results.dir, protein = protein,
    lesion.columns = lesion.columns, coordinate.columns = coordinate.columns,
    min.subjects = min.subjects, min.events = min.events,
    min.confidence = min.confidence,
    require.complete.event.mapping = require.complete.event.mapping,
    avoid.within.subject.overlap = avoid.within.subject.overlap,
    alpha = alpha, n.sim = n.sim, calibration.fraction = calibration.fraction,
    random.seed = random.seed, progress.every = progress.every
  )
  elapsed <- proc.time()[["elapsed"]] - started
  manifest <- data.frame(
    label = label, protein = protein, n_sim = n.sim,
    min_subjects = min.subjects, min_events = min.events,
    random_seed = random.seed,
    require_complete_event_mapping = require.complete.event.mapping,
    avoid_within_subject_overlap = avoid.within.subject.overlap,
    confidence_filter = !is.null(coordinate.columns$confidence),
    min_confidence = if (is.null(coordinate.columns$confidence)) NA_real_ else min.confidence,
    opportunity_weights = !is.null(coordinate.columns$opportunity),
    elapsed_seconds = round(elapsed, 3L),
    n_clusters = nrow(result$clusters),
    n_significant = sum(result$clusters$significant_any),
    results_dir = normalizePath(results.dir, winslash = "/", mustWork = FALSE),
    stringsAsFactors = FALSE
  )
  utils::write.csv(manifest, file.path(results.dir, "stability_run_manifest.csv"), row.names = FALSE)
  list(label = label, result = result, manifest = manifest)
}

grin3d_compare_stability_runs <- function(
    reference.run,
    comparison.run,
    output.file = NULL,
    reference.significant.only = TRUE) {
  matches <- grin3d_match_clusters(
    reference.run$result$clusters,
    comparison.run$result$clusters,
    reference.significant.only = reference.significant.only
  )
  matches$reference_run <- reference.run$label
  matches$comparison_run <- comparison.run$label
  matches <- matches[, c(
    "reference_run", "comparison_run",
    setdiff(names(matches), c("reference_run", "comparison_run"))
  )]
  if (!is.null(output.file)) {
    dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(matches, output.file, row.names = FALSE)
  }
  matches
}
