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

parse_residues <- function(x) {
  z <- suppressWarnings(as.integer(strsplit(as.character(x), ";", fixed = TRUE)[[1L]]))
  sort(unique(z[!is.na(z)]))
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || length(unique(x[keep])) < 2L || length(unique(y[keep])) < 2L) return(NA_real_)
  suppressWarnings(cor(x[keep], y[keep], method = "spearman"))
}

best_matches <- function(reference, comparison, reference.label, comparison.label) {
  if (!nrow(reference)) return(data.frame())
  cmp.sets <- lapply(comparison$residues, parse_residues)
  rows <- lapply(seq_len(nrow(reference)), function(i) {
    a <- parse_residues(reference$residues[[i]])
    if (!length(cmp.sets)) {
      return(data.frame(reference_cluster_id = reference$cluster_id[[i]], comparison_cluster_id = NA_character_,
                        shared_residues = 0L, union_residues = length(a), jaccard = 0, dice = 0,
                        overlap_coefficient = 0, reference_recall = 0, comparison_precision = 0,
                        size_difference = NA_integer_, left_boundary_shift = NA_integer_,
                        right_boundary_shift = NA_integer_, centroid_shift = NA_real_,
                        reference_p_any = reference$p_any_joint[[i]], comparison_p_any = NA_real_,
                        delta_p_any = NA_real_, log10_p_ratio = NA_real_,
                        reference_significant = reference$significant_any[[i]], comparison_significant = FALSE,
                        significance_flip = reference$significant_any[[i]], stringsAsFactors = FALSE))
    }
    shared <- vapply(cmp.sets, function(b) length(intersect(a, b)), integer(1L))
    unions <- vapply(cmp.sets, function(b) length(union(a, b)), integer(1L))
    jaccard <- shared / unions
    containment <- shared / pmin(length(a), lengths(cmp.sets))
    size.diff <- abs(length(a) - lengths(cmp.sets))
    j <- order(-jaccard, -containment, size.diff)[[1L]]
    b <- cmp.sets[[j]]
    p.ref <- reference$p_any_joint[[i]]
    p.cmp <- comparison$p_any_joint[[j]]
    data.frame(
      reference_cluster_id = reference$cluster_id[[i]], comparison_cluster_id = comparison$cluster_id[[j]],
      shared_residues = shared[[j]], union_residues = unions[[j]], jaccard = jaccard[[j]],
      dice = 2 * shared[[j]] / (length(a) + length(b)), overlap_coefficient = containment[[j]],
      reference_recall = shared[[j]] / length(a), comparison_precision = shared[[j]] / length(b),
      size_difference = length(b) - length(a), left_boundary_shift = min(b) - min(a),
      right_boundary_shift = max(b) - max(a), centroid_shift = mean(b) - mean(a),
      reference_p_any = p.ref, comparison_p_any = p.cmp, delta_p_any = p.cmp - p.ref,
      log10_p_ratio = log10(p.cmp) - log10(p.ref),
      reference_significant = reference$significant_any[[i]], comparison_significant = comparison$significant_any[[j]],
      significance_flip = reference$significant_any[[i]] != comparison$significant_any[[j]],
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result$reference_run <- reference.label
  result$comparison_run <- comparison.label
  result[, c("reference_run", "comparison_run", setdiff(names(result), c("reference_run", "comparison_run")))]
}

overlap_edges <- function(reference, comparison, reference.label, comparison.label) {
  if (!nrow(reference) || !nrow(comparison)) return(data.frame())
  a <- lapply(reference$residues, parse_residues)
  b <- lapply(comparison$residues, parse_residues)
  rows <- list(); k <- 0L
  for (i in seq_along(a)) for (j in seq_along(b)) {
    shared <- length(intersect(a[[i]], b[[j]]))
    if (!shared) next
    k <- k + 1L
    rows[[k]] <- data.frame(
      reference_run = reference.label, comparison_run = comparison.label,
      reference_cluster_id = reference$cluster_id[[i]], comparison_cluster_id = comparison$cluster_id[[j]],
      shared_residues = shared, jaccard = shared / length(union(a[[i]], b[[j]])),
      reference_recall = shared / length(a[[i]]), comparison_precision = shared / length(b[[j]]),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) data.frame() else do.call(rbind, rows)
}

repo.root <- find_repo_root()
analysis.root <- file.path(repo.root, "analysis", "mutation_plddt_trials")
summary.dir <- file.path(analysis.root, "results", "summary")
figure.dir <- file.path(analysis.root, "results", "figures")
dir.create(summary.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure.dir, recursive = TRUE, showWarnings = FALSE)
config <- read.csv(file.path(analysis.root, "configurations.csv"), stringsAsFactors = FALSE, check.names = FALSE)
required.config <- c("protein", "label", "n_sim", "plddt_threshold", "coordinate_file", "results_dir")
if (!all(required.config %in% names(config))) stop("configurations.csv is incomplete; rerun 01_run_sensitivity.R")

runs <- list()
for (i in seq_len(nrow(config))) {
  f <- file.path(config$results_dir[[i]], "GRIN3D_mutation_hotspot_results.rds")
  if (file.exists(f)) runs[[config$label[[i]]]] <- readRDS(f)
}
config <- config[config$label %in% names(runs), , drop = FALSE]
if (!nrow(config)) stop("No completed runs were found")

references <- do.call(rbind, lapply(split(config, config$protein), function(x) {
  threshold.rank <- ifelse(is.na(x$plddt_threshold), -Inf, x$plddt_threshold)
  x <- x[order(threshold.rank, -x$n_sim), , drop = FALSE]
  x[1L, , drop = FALSE]
}))

metrics <- list(); matches <- list(); edges <- list(); m <- e <- 0L
for (i in seq_len(nrow(config))) {
  cfg <- config[i, ]; result <- runs[[cfg$label]]
  ref.cfg <- references[references$protein == cfg$protein, ]; ref <- runs[[ref.cfg$label]]
  match.table <- best_matches(ref$clusters, result$clusters, ref.cfg$label, cfg$label)
  match.table$protein <- cfg$protein; match.table$n_sim <- cfg$n_sim
  match.table$plddt_threshold <- cfg$plddt_threshold
  m <- m + 1L; matches[[m]] <- match.table
  edge.table <- overlap_edges(ref$clusters, result$clusters, ref.cfg$label, cfg$label)
  if (nrow(edge.table)) {
    edge.table$protein <- cfg$protein; edge.table$n_sim <- cfg$n_sim
    edge.table$plddt_threshold <- cfg$plddt_threshold
    e <- e + 1L; edges[[e]] <- edge.table
  }
  manifest.file <- file.path(cfg$results_dir, "stability_run_manifest.csv")
  manifest <- read.csv(manifest.file, stringsAsFactors = FALSE)
  total.coordinates <- nrow(read.csv(cfg$coordinate_file, stringsAsFactors = FALSE))
  metrics[[i]] <- data.frame(
    protein = cfg$protein, label = cfg$label, reference_run = ref.cfg$label,
    n_sim = cfg$n_sim, plddt_threshold = cfg$plddt_threshold,
    n_calibration = result$n_calibration, n_evaluation = result$n_evaluation,
    minimum_empirical_p = 1 / (result$n_evaluation + 1),
    eligible_residues = nrow(result$coordinates), total_coordinate_residues = total.coordinates,
    coordinate_retention = nrow(result$coordinates) / total.coordinates,
    mapped_events = nrow(result$mapped_events), excluded_events = nrow(result$excluded_events),
    mapped_subjects = length(unique(result$mapped_events$subject_id)),
    n_clusters = nrow(result$clusters), n_significant = sum(result$clusters$significant_any),
    mean_best_jaccard = mean(match.table$jaccard), median_best_jaccard = median(match.table$jaccard),
    reference_recovery_jaccard_0_8 = mean(match.table$jaccard >= 0.8),
    pvalue_spearman = safe_spearman(match.table$reference_p_any, match.table$comparison_p_any),
    significance_flips = sum(match.table$significance_flip), elapsed_seconds = manifest$elapsed_seconds[[1L]],
    stringsAsFactors = FALSE
  )
}

run.metrics <- do.call(rbind, metrics)
match.table <- do.call(rbind, matches)
edge.table <- if (length(edges)) do.call(rbind, edges) else data.frame()

write.csv(run.metrics, file.path(summary.dir, "run_metrics.csv"), row.names = FALSE)
write.csv(match.table, file.path(summary.dir, "reference_cluster_matches.csv"), row.names = FALSE)
write.csv(edge.table, file.path(summary.dir, "cluster_overlap_edges.csv"), row.names = FALSE)

# Summarize each reference cluster across all completed configurations. The
# score prioritizes reproducible membership, recurring significance, and
# subject support; it is descriptive and must not be read as a p-value.
stability.rows <- list(); k <- 0L
for (protein in unique(config$protein)) {
  ref.label <- references$label[references$protein == protein][[1L]]
  ref.clusters <- runs[[ref.label]]$clusters
  protein.matches <- match.table[match.table$protein == protein, , drop = FALSE]
  for (i in seq_len(nrow(ref.clusters))) {
    z <- protein.matches[protein.matches$reference_cluster_id == ref.clusters$cluster_id[[i]], , drop = FALSE]
    k <- k + 1L
    stability.rows[[k]] <- data.frame(
      protein = protein, reference_run = ref.label,
      reference_cluster_id = ref.clusters$cluster_id[[i]], residues = ref.clusters$residues[[i]],
      residue_min = ref.clusters$residue_min[[i]], residue_max = ref.clusters$residue_max[[i]],
      n_residues = ref.clusters$n_residues[[i]], n_subjects = ref.clusters$n_subjects[[i]],
      n_events = ref.clusters$n_events[[i]], tree_sources = ref.clusters$tree_sources[[i]],
      reference_p_any = ref.clusters$p_any_joint[[i]],
      reference_significant = ref.clusters$significant_any[[i]],
      mean_jaccard = mean(z$jaccard), min_jaccard = min(z$jaccard),
      highly_stable_fraction = mean(z$jaccard >= 0.8),
      exact_membership_fraction = mean(z$jaccard == 1),
      significance_fraction = mean(z$comparison_significant),
      mean_abs_delta_p = mean(abs(z$delta_p_any), na.rm = TRUE),
      max_abs_boundary_shift = suppressWarnings(max(abs(c(z$left_boundary_shift, z$right_boundary_shift)), na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }
}
cluster.stability <- do.call(rbind, stability.rows)
cluster.stability$max_abs_boundary_shift[!is.finite(cluster.stability$max_abs_boundary_shift)] <- NA_real_
cluster.stability$support_percentile <- ave(
  cluster.stability$n_subjects, cluster.stability$protein,
  FUN = function(x) rank(x, ties.method = "average") / length(x)
)
cluster.stability$priority_score <- round(
  0.50 * cluster.stability$mean_jaccard +
    0.25 * cluster.stability$significance_fraction +
    0.25 * cluster.stability$support_percentile,
  6L
)
cluster.stability <- cluster.stability[order(
  cluster.stability$protein, -cluster.stability$priority_score,
  cluster.stability$reference_p_any
), , drop = FALSE]
write.csv(cluster.stability, file.path(summary.dir, "reference_cluster_stability.csv"), row.names = FALSE)

threshold_name <- function(x) ifelse(is.na(x), "no pLDDT", paste0("pLDDT >= ", format(x, trim = TRUE)))

plot_metric_lines <- function(data, column, ylab, main, ylim = NULL) {
  thresholds <- unique(data$plddt_threshold)
  threshold.keys <- ifelse(is.na(thresholds), "NA", as.character(thresholds))
  colors <- seq_along(thresholds)
  if (is.null(ylim)) ylim <- range(data[[column]], finite = TRUE)
  if (!all(is.finite(ylim)) || diff(ylim) == 0) ylim <- c(ylim[[1L]] - 0.5, ylim[[1L]] + 0.5)
  plot(range(data$n_sim), ylim, type = "n", log = "x", xlab = "Null-simulation trials",
       ylab = ylab, main = main)
  for (j in seq_along(thresholds)) {
    key <- ifelse(is.na(data$plddt_threshold), "NA", as.character(data$plddt_threshold))
    z <- data[key == threshold.keys[[j]], , drop = FALSE]
    z <- z[order(z$n_sim), , drop = FALSE]
    lines(z$n_sim, z[[column]], type = "b", pch = 19L, col = colors[[j]], lwd = 2)
  }
  legend("topright", legend = threshold_name(thresholds), col = colors, pch = 19L, lty = 1L, cex = 0.75)
}

for (protein in unique(run.metrics$protein)) {
  z <- run.metrics[run.metrics$protein == protein, , drop = FALSE]
  png(file.path(figure.dir, paste0(protein, "_global_sensitivity.png")), width = 1800, height = 1400, res = 180)
  old <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  plot_metric_lines(z, "mean_best_jaccard", "Mean best-match Jaccard", paste(protein, "membership stability"), c(0, 1))
  plot_metric_lines(z, "coordinate_retention", "Coordinate retention", "Structure retained", c(0, 1))
  plot_metric_lines(z, "mapped_events", "Mapped events", "Mutation mapping retention")
  plot_metric_lines(z, "n_significant", "Significant clusters", "Significance calls",
                    c(0, max(1, z$n_significant)))
  par(old); dev.off()
}

# Heatmap of the best-match Jaccard for the highest-priority reference clusters.
for (protein in unique(config$protein)) {
  priority <- cluster.stability[cluster.stability$protein == protein, , drop = FALSE]
  priority <- head(priority, 25L)
  z <- match.table[match.table$protein == protein &
                     match.table$reference_cluster_id %in% priority$reference_cluster_id, , drop = FALSE]
  run.order <- config$label[config$protein == protein][order(
    config$plddt_threshold[config$protein == protein], config$n_sim[config$protein == protein], na.last = FALSE
  )]
  matrix.j <- matrix(NA_real_, nrow = nrow(priority), ncol = length(run.order),
                     dimnames = list(priority$reference_cluster_id, run.order))
  for (i in seq_len(nrow(z))) matrix.j[z$reference_cluster_id[[i]], z$comparison_run[[i]]] <- z$jaccard[[i]]
  png(file.path(figure.dir, paste0(protein, "_cluster_jaccard_heatmap.png")),
      width = max(1800, 190 * ncol(matrix.j)), height = max(1200, 65 * nrow(matrix.j)), res = 180)
  old <- par(mar = c(12, 7, 4, 2))
  image(seq_len(ncol(matrix.j)), seq_len(nrow(matrix.j)), t(matrix.j[nrow(matrix.j):1, , drop = FALSE]),
        zlim = c(0, 1), col = hcl.colors(20, "YlOrRd", rev = TRUE), axes = FALSE,
        xlab = "Configuration", ylab = "", main = paste(protein, "best-match Jaccard"))
  axis(1, at = seq_len(ncol(matrix.j)), labels = colnames(matrix.j), las = 2, cex.axis = 0.65)
  axis(2, at = seq_len(nrow(matrix.j)), labels = rev(rownames(matrix.j)), las = 2, cex.axis = 0.7)
  box(); par(old); dev.off()
}

# A compact atlas puts every trial count in one view. Numeric labels make the
# Jaccard mapping explicit, while separate panels prevent pLDDT and trial-count
# effects from being visually conflated.
for (protein in unique(config$protein)) {
  priority <- head(cluster.stability[cluster.stability$protein == protein, ], 20L)
  protein.config <- config[config$protein == protein, , drop = FALSE]
  threshold.keys <- unique(ifelse(is.na(protein.config$plddt_threshold), "NA", as.character(protein.config$plddt_threshold)))
  png(file.path(figure.dir, paste0(protein, "_all_trials_cluster_overlap_atlas.png")),
      width = 2600, height = max(1400, 1050 * length(threshold.keys)), res = 180)
  old <- par(mfrow = c(length(threshold.keys), 1L), mar = c(5, 7, 4, 5))
  palette <- hcl.colors(101, "YlOrRd", rev = TRUE)
  for (threshold.key in threshold.keys) {
    keep <- if (threshold.key == "NA") is.na(protein.config$plddt_threshold) else
      as.character(protein.config$plddt_threshold) == threshold.key
    panel.config <- protein.config[keep, , drop = FALSE]
    panel.config <- panel.config[order(panel.config$n_sim), , drop = FALSE]
    q <- match.table[match.table$protein == protein &
                       match.table$reference_cluster_id %in% priority$reference_cluster_id &
                       match.table$comparison_run %in% panel.config$label, , drop = FALSE]
    matrix.j <- matrix(NA_real_, nrow = nrow(priority), ncol = nrow(panel.config),
                       dimnames = list(priority$reference_cluster_id, panel.config$n_sim))
    for (i in seq_len(nrow(q))) matrix.j[q$reference_cluster_id[[i]], as.character(q$n_sim[[i]])] <- q$jaccard[[i]]
    plot.matrix <- matrix.j[nrow(matrix.j):1L, , drop = FALSE]
    image(seq_len(ncol(plot.matrix)), seq_len(nrow(plot.matrix)), t(plot.matrix),
          zlim = c(0, 1), col = palette, axes = FALSE, xlab = "Null-simulation trials", ylab = "",
          main = paste(protein, if (threshold.key == "NA") "no pLDDT available" else paste0("pLDDT >= ", threshold.key)))
    axis(1, at = seq_len(ncol(plot.matrix)), labels = colnames(plot.matrix))
    axis(2, at = seq_len(nrow(plot.matrix)), labels = rownames(plot.matrix), las = 2, cex.axis = 0.65)
    for (x in seq_len(ncol(plot.matrix))) for (y in seq_len(nrow(plot.matrix))) {
      value <- plot.matrix[y, x]
      if (is.finite(value)) text(x, y, sprintf("%.2f", value), cex = 0.55,
                                 col = if (value >= 0.75 || value <= 0.15) "white" else "black")
    }
    legend("right", inset = c(-0.12, 0), xpd = NA,
           legend = sprintf("%.2f", seq(0, 1, by = 0.25)),
           fill = palette[1L + 100L * seq(0, 1, by = 0.25)],
           title = "Jaccard", cex = 0.65, bty = "n")
    box()
  }
  par(old); dev.off()
}

# Symmetric configuration similarity summarizes the entire candidate-cluster
# landscape. Each cell averages directed best-match Jaccard in both directions.
configuration.overlaps <- list(); overlap.index <- 0L
for (protein in unique(config$protein)) {
  labels <- config$label[config$protein == protein]
  short.labels <- vapply(labels, function(label) {
    cfg <- config[config$label == label, ]
    paste0("t", cfg$n_sim, "_p", ifelse(is.na(cfg$plddt_threshold), "NA", cfg$plddt_threshold))
  }, character(1L))
  similarity <- matrix(NA_real_, length(labels), length(labels), dimnames = list(short.labels, short.labels))
  cluster.sets <- setNames(lapply(labels, function(label) {
    lapply(runs[[label]]$clusters$residues, parse_residues)
  }), labels)
  membership.signature <- vapply(labels, function(label) {
    paste(vapply(cluster.sets[[label]], paste, collapse = ";", character(1L)), collapse = "|")
  }, character(1L))
  membership.group <- match(membership.signature, unique(membership.signature))
  similarity.cache <- new.env(parent = emptyenv())
  for (i in seq_along(labels)) for (j in seq_along(labels)) {
    cache.key <- paste(membership.group[[i]], membership.group[[j]], sep = "::")
    if (!exists(cache.key, envir = similarity.cache, inherits = FALSE)) {
      a.sets <- cluster.sets[[labels[[i]]]]; b.sets <- cluster.sets[[labels[[j]]]]
      directed <- function(from, to) mean(vapply(from, function(x) {
        max(vapply(to, function(y) length(intersect(x, y)) / length(union(x, y)), numeric(1L)))
      }, numeric(1L)))
      assign(cache.key, mean(c(directed(a.sets, b.sets), directed(b.sets, a.sets))),
             envir = similarity.cache)
    }
    similarity[i, j] <- get(cache.key, envir = similarity.cache, inherits = FALSE)
  }
  write.csv(cbind(configuration = rownames(similarity), as.data.frame(similarity, check.names = FALSE)),
            file.path(summary.dir, paste0(protein, "_configuration_overlap_matrix.csv")), row.names = FALSE)
  png(file.path(figure.dir, paste0(protein, "_configuration_overlap_matrix.png")),
      width = max(1800, 145 * ncol(similarity)), height = max(1600, 145 * nrow(similarity)), res = 180)
  old <- par(mar = c(9, 9, 4, 2))
  image(seq_len(ncol(similarity)), seq_len(nrow(similarity)), t(similarity[nrow(similarity):1L, , drop = FALSE]),
        zlim = c(0, 1), col = hcl.colors(101, "YlOrRd", rev = TRUE), axes = FALSE,
        xlab = "Configuration", ylab = "", main = paste(protein, "configuration-level cluster overlap"))
  axis(1, at = seq_len(ncol(similarity)), labels = colnames(similarity), las = 2, cex.axis = 0.6)
  axis(2, at = seq_len(nrow(similarity)), labels = rev(rownames(similarity)), las = 2, cex.axis = 0.6)
  box(); par(old); dev.off()
}

# Residue maps make expansions, contractions, and lost clusters inspectable.
for (protein in unique(config$protein)) {
  priority <- head(cluster.stability[cluster.stability$protein == protein, ], 8L)
  z <- match.table[match.table$protein == protein, , drop = FALSE]
  pdf(file.path(figure.dir, paste0(protein, "_cluster_residue_maps.pdf")), width = 12, height = 7)
  for (cluster.id in priority$reference_cluster_id) {
    one <- z[z$reference_cluster_id == cluster.id, , drop = FALSE]
    one <- one[order(one$plddt_threshold, one$n_sim, na.last = FALSE), , drop = FALSE]
    sets <- lapply(one$comparison_run, function(label) {
      id <- one$comparison_cluster_id[one$comparison_run == label][[1L]]
      table <- runs[[label]]$clusters
      if (is.na(id) || !id %in% table$cluster_id) integer() else parse_residues(table$residues[table$cluster_id == id][[1L]])
    })
    reference.set <- parse_residues(priority$residues[priority$reference_cluster_id == cluster.id][[1L]])
    limits <- range(c(reference.set, unlist(sets)))
    old <- par(mar = c(5, 13, 4, 2))
    plot(limits, c(0.5, length(sets) + 0.5), type = "n", yaxt = "n", xlab = "Residue position", ylab = "",
         main = paste(protein, cluster.id, "membership across configurations"))
    axis(2, at = seq_along(sets), labels = paste0(one$comparison_run, "  J=", sprintf("%.2f", one$jaccard)),
         las = 2, cex.axis = 0.62)
    abline(v = range(reference.set), col = "grey80", lty = 3)
    for (i in seq_along(sets)) {
      segments(min(reference.set), i, max(reference.set), i, col = "grey85", lwd = 4)
      points(sets[[i]], rep(i, length(sets[[i]])), pch = 15, col = ifelse(sets[[i]] %in% reference.set, "#2166AC", "#B2182B"), cex = 0.75)
    }
    legend("topright", c("shared/reference residue", "added residue"), pch = 15,
           col = c("#2166AC", "#B2182B"), bty = "n", cex = 0.8)
    par(old)
  }
  dev.off()
}

# Bipartite overlap diagrams retain one-to-many edges, making candidate splits
# and merges visible rather than forcing every cluster into a one-to-one match.
for (protein in unique(config$protein)) {
  ref.label <- references$label[references$protein == protein][[1L]]
  priority.ids <- head(
    cluster.stability$reference_cluster_id[cluster.stability$protein == protein],
    20L
  )
  comparisons <- setdiff(config$label[config$protein == protein], ref.label)
  pdf(file.path(figure.dir, paste0(protein, "_cluster_overlap_transitions.pdf")), width = 12, height = 8)
  for (comparison in comparisons) {
    q <- edge.table[
      edge.table$protein == protein & edge.table$comparison_run == comparison &
        edge.table$reference_cluster_id %in% priority.ids & edge.table$jaccard >= 0.20,
      , drop = FALSE
    ]
    if (!nrow(q)) {
      plot.new(); title(main = paste(protein, comparison), sub = "No overlap edges with Jaccard >= 0.20")
      next
    }
    left.ids <- priority.ids[priority.ids %in% q$reference_cluster_id]
    right.ids <- unique(q$comparison_cluster_id[order(-q$jaccard)])
    left.y <- setNames(seq_along(left.ids), left.ids)
    right.y <- setNames(seq(1, length(left.ids), length.out = length(right.ids)), right.ids)
    old <- par(mar = c(2, 11, 4, 11), xpd = NA)
    plot(c(0, 1), c(0.5, length(left.ids) + 0.5), type = "n", axes = FALSE,
         xlab = "", ylab = "", main = paste(protein, "cluster overlap:", comparison))
    edge.colors <- hcl.colors(100, "YlOrRd", rev = TRUE)[pmax(1L, ceiling(q$jaccard * 100))]
    segments(0, left.y[q$reference_cluster_id], 1, right.y[q$comparison_cluster_id],
             col = adjustcolor(edge.colors, alpha.f = 0.7), lwd = 1 + 5 * q$jaccard)
    points(rep(0, length(left.ids)), left.y, pch = 21, bg = "#2166AC", cex = 1.2)
    points(rep(1, length(right.ids)), right.y, pch = 21, bg = "#B2182B", cex = 1.2)
    text(0, left.y, left.ids, pos = 2, cex = 0.68)
    text(1, right.y, right.ids, pos = 4, cex = 0.68)
    mtext(ref.label, side = 3, at = 0, line = 0.4, font = 2)
    mtext(comparison, side = 3, at = 1, line = 0.4, font = 2)
    par(old)
  }
  dev.off()
}

# P-value trajectories show Monte Carlo variation separately from membership.
for (protein in unique(config$protein)) {
  priority <- head(cluster.stability[cluster.stability$protein == protein, ], 6L)
  z <- match.table[match.table$protein == protein &
                     match.table$reference_cluster_id %in% priority$reference_cluster_id, , drop = FALSE]
  pdf(file.path(figure.dir, paste0(protein, "_pvalue_trajectories.pdf")), width = 9, height = 6)
  for (cluster.id in priority$reference_cluster_id) {
    one <- z[z$reference_cluster_id == cluster.id, , drop = FALSE]
    thresholds <- unique(one$plddt_threshold); colors <- seq_along(thresholds)
    y <- -log10(one$comparison_p_any)
    plot(range(one$n_sim), range(y, finite = TRUE), type = "n", log = "x",
         xlab = "Null-simulation trials", ylab = expression(-log[10](p[any])),
         main = paste(protein, cluster.id, "empirical p-value"))
    for (j in seq_along(thresholds)) {
      keep <- if (is.na(thresholds[[j]])) is.na(one$plddt_threshold) else one$plddt_threshold == thresholds[[j]]
      q <- one[keep, ]; q <- q[order(q$n_sim), ]
      lines(q$n_sim, -log10(q$comparison_p_any), type = "b", pch = 19, col = colors[[j]], lwd = 2)
    }
    abline(h = -log10(0.05), lty = 2, col = "grey40")
    legend("topright", threshold_name(thresholds), col = colors, pch = 19, lty = 1, bty = "n")
  }
  dev.off()
}

message("Summary tables written to: ", summary.dir)
message("Figures written to: ", figure.dir)
