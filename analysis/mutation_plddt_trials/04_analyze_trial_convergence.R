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

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || length(unique(x[keep])) < 2L || length(unique(y[keep])) < 2L) return(NA_real_)
  suppressWarnings(cor(x[keep], y[keep], method = "spearman"))
}

threshold_label <- function(x) ifelse(is.na(x), "no pLDDT", paste0("pLDDT >= ", x))
threshold_key <- function(x) ifelse(is.na(x), "NA", as.character(x))

repo.root <- find_repo_root()
analysis.root <- file.path(repo.root, "analysis", "mutation_plddt_trials")
summary.dir <- file.path(analysis.root, "results", "summary")
figure.dir <- file.path(analysis.root, "results", "figures")
config <- read.csv(file.path(analysis.root, "configurations.csv"), stringsAsFactors = FALSE)

result.files <- file.path(config$results_dir, "GRIN3D_mutation_hotspot_results.rds")
if (any(!file.exists(result.files))) stop("Complete the configured runs before convergence analysis")
runs <- setNames(lapply(result.files, readRDS), config$label)

keys <- paste(config$protein, threshold_key(config$plddt_threshold), sep = "::")
groups <- split(seq_len(nrow(config)), keys)
metric.rows <- list(); cluster.rows <- list(); earliest.rows <- list()
metric.index <- cluster.index <- earliest.index <- 0L

for (indices in groups) {
  group.config <- config[indices, , drop = FALSE]
  group.config <- group.config[order(group.config$n_sim), , drop = FALSE]
  reference.config <- group.config[group.config$n_sim == max(group.config$n_sim), , drop = FALSE]
  reference <- runs[[reference.config$label]]$clusters
  reference <- reference[order(reference$candidate_key), , drop = FALSE]
  reference.significant.keys <- sort(reference$candidate_key[reference$significant_any])
  significance.set.matches <- logical(nrow(group.config))

  for (i in seq_len(nrow(group.config))) {
    cfg <- group.config[i, ]
    comparison <- runs[[cfg$label]]$clusters
    comparison <- comparison[match(reference$candidate_key, comparison$candidate_key), , drop = FALSE]
    if (anyNA(comparison$candidate_key)) stop("Cluster membership differs within a fixed pLDDT group: ", cfg$label)
    delta <- comparison$p_any_joint - reference$p_any_joint
    significance.agreement <- mean(comparison$significant_any == reference$significant_any)
    comparison.significant.keys <- sort(comparison$candidate_key[comparison$significant_any])
    significance.set.matches[[i]] <- identical(comparison.significant.keys, reference.significant.keys)

    metric.index <- metric.index + 1L
    metric.rows[[metric.index]] <- data.frame(
      protein = cfg$protein, plddt_threshold = cfg$plddt_threshold,
      n_sim = cfg$n_sim, reference_n_sim = reference.config$n_sim,
      n_clusters = nrow(reference), n_significant = sum(comparison$significant_any),
      reference_n_significant = sum(reference$significant_any),
      mean_absolute_p_difference = mean(abs(delta)),
      median_absolute_p_difference = median(abs(delta)),
      root_mean_square_p_difference = sqrt(mean(delta^2)),
      maximum_absolute_p_difference = max(abs(delta)),
      pvalue_spearman = safe_spearman(comparison$p_any_joint, reference$p_any_joint),
      significance_agreement = significance.agreement,
      significance_flips = sum(comparison$significant_any != reference$significant_any),
      exact_significance_set_match = significance.set.matches[[i]],
      stringsAsFactors = FALSE
    )

    cluster.index <- cluster.index + 1L
    cluster.rows[[cluster.index]] <- data.frame(
      protein = cfg$protein, plddt_threshold = cfg$plddt_threshold,
      n_sim = cfg$n_sim, reference_n_sim = reference.config$n_sim,
      candidate_key = reference$candidate_key,
      cluster_id = comparison$cluster_id,
      residues = comparison$residues,
      n_subjects = comparison$n_subjects,
      n_events = comparison$n_events,
      p_any_joint = comparison$p_any_joint,
      reference_p_any_joint = reference$p_any_joint,
      absolute_p_difference = abs(delta),
      significant_any = comparison$significant_any,
      reference_significant_any = reference$significant_any,
      significance_flip = comparison$significant_any != reference$significant_any,
      stringsAsFactors = FALSE
    )
  }

  persistent <- vapply(seq_along(significance.set.matches), function(i) {
    all(significance.set.matches[seq.int(i, length(significance.set.matches))])
  }, logical(1L))
  earliest.index <- earliest.index + 1L
  earliest.rows[[earliest.index]] <- data.frame(
    protein = group.config$protein[[1L]],
    plddt_threshold = group.config$plddt_threshold[[1L]],
    earliest_persistent_significance_match = if (any(persistent)) min(group.config$n_sim[persistent]) else NA_integer_,
    reference_n_significant = sum(reference$significant_any),
    interpretation = if (sum(reference$significant_any) == 0L)
      "Trivial agreement: no significant clusters in the 1,000-trial reference"
    else "First trial count whose significance set matches 1,000 trials and remains matched",
    stringsAsFactors = FALSE
  )
}

metrics <- do.call(rbind, metric.rows)
clusters <- do.call(rbind, cluster.rows)
earliest <- do.call(rbind, earliest.rows)
metrics <- metrics[order(metrics$protein, metrics$plddt_threshold, metrics$n_sim, na.last = FALSE), ]
clusters <- clusters[order(clusters$protein, clusters$plddt_threshold, clusters$n_sim, clusters$candidate_key,
                           na.last = FALSE), ]
earliest <- earliest[order(earliest$protein, earliest$plddt_threshold, na.last = FALSE), ]

write.csv(metrics, file.path(summary.dir, "trial_convergence_metrics.csv"), row.names = FALSE)
write.csv(clusters, file.path(summary.dir, "cluster_pvalue_convergence.csv"), row.names = FALSE)
write.csv(earliest, file.path(summary.dir, "earliest_significance_convergence.csv"), row.names = FALSE)

plot_one_metric <- function(data, column, ylab, ylim = NULL) {
  if (is.null(ylim)) ylim <- range(data[[column]], finite = TRUE)
  if (!all(is.finite(ylim)) || diff(ylim) == 0) ylim <- c(ylim[[1L]] - 0.05, ylim[[1L]] + 0.05)
  plot(data$n_sim, data[[column]], type = "b", pch = 19, lwd = 2,
       xlab = "Null-simulation trials", ylab = ylab, ylim = ylim)
  grid(col = "grey90")
}

for (protein in unique(metrics$protein)) {
  protein.data <- metrics[metrics$protein == protein, , drop = FALSE]
  thresholds <- unique(threshold_key(protein.data$plddt_threshold))
  png(file.path(figure.dir, paste0(protein, "_trial_convergence_overview.png")),
      width = 2400, height = max(1000, 700 * length(thresholds)), res = 180)
  old <- par(mfrow = c(length(thresholds), 3L), mar = c(4, 4, 3, 1))
  for (key in thresholds) {
    keep <- threshold_key(protein.data$plddt_threshold) == key
    z <- protein.data[keep, , drop = FALSE]
    plot_one_metric(z, "mean_absolute_p_difference", "Mean |p - p1000|")
    title(main = paste(protein, threshold_label(z$plddt_threshold[[1L]])))
    plot_one_metric(z, "pvalue_spearman", "P-value Spearman", c(0, 1))
    abline(h = 0.9, lty = 2, col = "grey40")
    plot_one_metric(z, "significance_agreement", "Significance agreement", c(0, 1))
  }
  par(old); dev.off()
}

# Show which exact clusters cross the significance threshold as trials grow.
for (protein in unique(clusters$protein)) {
  protein.data <- clusters[clusters$protein == protein, , drop = FALSE]
  thresholds <- unique(threshold_key(protein.data$plddt_threshold))
  png(file.path(figure.dir, paste0(protein, "_significance_by_trial.png")),
      width = 2400, height = max(1000, 750 * length(thresholds)), res = 180)
  old <- par(mfrow = c(length(thresholds), 1L), mar = c(5, 8, 4, 2))
  for (key in thresholds) {
    z <- protein.data[threshold_key(protein.data$plddt_threshold) == key, , drop = FALSE]
    candidate.ids <- unique(z$candidate_key[z$significant_any | z$reference_significant_any])
    if (!length(candidate.ids)) {
      plot.new(); title(main = paste(protein, threshold_label(z$plddt_threshold[[1L]])),
                        sub = "No significant clusters at any evaluated trial count")
      next
    }
    trials <- sort(unique(z$n_sim))
    matrix.significant <- matrix(0, nrow = length(candidate.ids), ncol = length(trials),
                                 dimnames = list(candidate.ids, trials))
    for (i in seq_len(nrow(z))) {
      if (z$candidate_key[[i]] %in% candidate.ids)
        matrix.significant[z$candidate_key[[i]], as.character(z$n_sim[[i]])] <- as.integer(z$significant_any[[i]])
    }
    plot.matrix <- matrix.significant[nrow(matrix.significant):1L, , drop = FALSE]
    image(seq_len(ncol(plot.matrix)), seq_len(nrow(plot.matrix)), t(plot.matrix),
          col = c("grey92", "#B2182B"), zlim = c(0, 1), axes = FALSE,
          xlab = "Null-simulation trials", ylab = "",
          main = paste(protein, threshold_label(z$plddt_threshold[[1L]])))
    axis(1, at = seq_len(ncol(plot.matrix)), labels = colnames(plot.matrix))
    labels <- rownames(plot.matrix)
    labels <- ifelse(nchar(labels) > 40L, paste0(substr(labels, 1L, 37L), "..."), labels)
    axis(2, at = seq_len(nrow(plot.matrix)), labels = labels, las = 2, cex.axis = 0.55)
    box()
  }
  par(old); dev.off()
}

message("Trial-convergence tables and figures written successfully.")
