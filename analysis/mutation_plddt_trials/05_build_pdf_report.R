#!/usr/bin/env Rscript

# Build a self-contained PDF report from completed sensitivity-analysis CSVs.
# This script uses base R graphics only and does not rerun hotspot simulations.

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
report.dir <- file.path(analysis.root, "report")
dir.create(report.dir, recursive = TRUE, showWarnings = FALSE)
report.file <- file.path(report.dir, "mutation_hotspot_stability_report.pdf")

required.files <- file.path(summary.dir, c(
  "run_metrics.csv", "reference_cluster_matches.csv",
  "reference_cluster_stability.csv", "trial_convergence_metrics.csv",
  "earliest_significance_convergence.csv", "validation_checks.csv"
))
if (any(!file.exists(required.files))) {
  stop("Missing summary outputs. Run scripts 02, 03, and 04 before building the report.")
}

run.metrics <- read.csv(required.files[[1L]], stringsAsFactors = FALSE, check.names = FALSE)
matches <- read.csv(required.files[[2L]], stringsAsFactors = FALSE, check.names = FALSE)
cluster.stability <- read.csv(required.files[[3L]], stringsAsFactors = FALSE, check.names = FALSE)
convergence <- read.csv(required.files[[4L]], stringsAsFactors = FALSE, check.names = FALSE)
earliest <- read.csv(required.files[[5L]], stringsAsFactors = FALSE, check.names = FALSE)
validation <- read.csv(required.files[[6L]], stringsAsFactors = FALSE, check.names = FALSE)
config <- read.csv(file.path(analysis.root, "configurations.csv"), stringsAsFactors = FALSE)

configured_file <- function(protein, column) {
  values <- unique(config[config$protein == protein, column])
  if (length(values) != 1L || !file.exists(values)) {
    stop("Could not resolve one existing ", column, " for ", protein)
  }
  values[[1L]]
}

colors <- list(
  navy = "#17365D", blue = "#2166AC", light.blue = "#67A9CF",
  red = "#B2182B", orange = "#EF8A62", green = "#1B7837",
  gold = "#D9A300", grey = "#5B6573", light.grey = "#EEF1F5",
  dark = "#1E2530"
)
threshold.colors <- c(`0` = colors$navy, `70` = colors$orange, `90` = colors$red, `NA` = colors$grey)
page.number <- 0L

page_setup <- function(title, subtitle = NULL) {
  page.number <<- page.number + 1L
  par(mar = c(3.1, 4.2, 4.5, 1.5), oma = c(0, 0, 0, 0), family = "sans")
  plot.new()
  plot.window(c(0, 1), c(0, 1))
  rect(0, 0.93, 1, 1, col = colors$navy, border = NA)
  text(0.035, 0.965, title, adj = c(0, 0.5), col = "white", cex = 1.35, font = 2)
  if (!is.null(subtitle)) text(0.035, 0.915, subtitle, adj = c(0, 1), col = colors$grey, cex = 0.78)
  text(0.98, 0.018, paste("GRIN3D mutation prototype evaluation  |  page", page.number),
       adj = c(1, 0), cex = 0.58, col = colors$grey)
}

wrapped_text <- function(text.value, x, y, width = 100L, cex = 0.88,
                         line.height = 0.035, col = colors$dark, font = 1L) {
  lines <- unlist(lapply(strsplit(text.value, "\n", fixed = TRUE)[[1L]], function(z) {
    if (!nzchar(z)) "" else strwrap(z, width = width)
  }))
  for (line in lines) {
    text(x, y, line, adj = c(0, 1), cex = cex, col = col, font = font)
    y <- y - line.height
  }
  invisible(y)
}

bullet_list <- function(items, x = 0.07, y = 0.82, width = 100L, cex = 0.83,
                        line.height = 0.032, gap = 0.014) {
  for (item in items) {
    text(x, y, intToUtf8(0x2022), adj = c(0, 1), cex = cex, col = colors$blue)
    lines <- strwrap(item, width = width)
    for (line in lines) {
      text(x + 0.025, y, line, adj = c(0, 1), cex = cex, col = colors$dark)
      y <- y - line.height
    }
    y <- y - gap
  }
  invisible(y)
}

text_table <- function(data, x = 0.05, y = 0.80, widths = NULL, row.height = 0.052,
                       cex = 0.72, header.col = colors$navy, wrap.widths = NULL,
                       line.height = 0.023) {
  data[] <- lapply(data, as.character)
  if (is.null(widths)) widths <- rep(0.9 / ncol(data), ncol(data))
  if (is.null(wrap.widths)) wrap.widths <- rep(Inf, ncol(data))
  if (length(wrap.widths) != ncol(data)) stop("wrap.widths must match the table column count")
  xpos <- x + c(0, cumsum(widths[-length(widths)]))
  rect(x - 0.01, y - 0.035, x + sum(widths), y + 0.022, col = header.col, border = NA)
  for (j in seq_len(ncol(data))) {
    text(xpos[[j]], y, names(data)[[j]], adj = c(0, 0.5), col = "white", font = 2, cex = cex)
  }
  y <- y - row.height
  for (i in seq_len(nrow(data))) {
    cells <- lapply(seq_len(ncol(data)), function(j) {
      if (is.infinite(wrap.widths[[j]])) data[[j]][[i]] else strwrap(data[[j]][[i]], width = wrap.widths[[j]])
    })
    cell.height <- max(row.height, max(lengths(cells)) * line.height + 0.014)
    if (i %% 2L == 0L) rect(x - 0.01, y - cell.height / 2, x + sum(widths), y + cell.height / 2,
                            col = colors$light.grey, border = NA)
    for (j in seq_len(ncol(data))) {
      lines <- cells[[j]]
      line.y <- y + ((length(lines) - 1L) * line.height / 2)
      for (line in lines) {
        text(xpos[[j]], line.y, line, adj = c(0, 0.5), cex = cex, col = colors$dark)
        line.y <- line.y - line.height
      }
    }
    y <- y - cell.height
  }
  invisible(y)
}

threshold_key <- function(x) ifelse(is.na(x), "NA", as.character(x))
threshold_label <- function(x) ifelse(is.na(x), "No pLDDT", paste0("pLDDT >= ", x))
get_threshold_color <- function(x) {
  key <- threshold_key(x)
  unname(ifelse(key == "NA", colors$grey, threshold.colors[key]))
}

draw_page_header <- function(title, subtitle = NULL) {
  mtext(title, side = 3, outer = TRUE, line = -1.8, font = 2, cex = 1.25, col = colors$navy)
  if (!is.null(subtitle)) mtext(subtitle, side = 3, outer = TRUE, line = -3.1,
                               cex = 0.72, col = colors$grey)
  mtext(paste("GRIN3D mutation prototype evaluation  |  page", page.number),
        side = 1, outer = TRUE, line = -1.2, adj = 1, cex = 0.55, col = colors$grey)
}

new_plot_page <- function(title, subtitle = NULL, mfrow = c(1, 1), mar = c(4, 4, 3, 1)) {
  page.number <<- page.number + 1L
  pending.plot.title <<- title
  pending.plot.subtitle <<- subtitle
  par(mfrow = mfrow, mar = mar, oma = c(2, 1, 4, 1), family = "sans")
}

pending.plot.title <- NULL
pending.plot.subtitle <- NULL
finish_plot_page <- function() draw_page_header(pending.plot.title, pending.plot.subtitle)

parse_residues <- function(x) {
  z <- suppressWarnings(as.integer(strsplit(as.character(x), ";", fixed = TRUE)[[1L]]))
  sort(unique(z[!is.na(z)]))
}

plot_jaccard_panel <- function(protein, threshold, n.clusters = 15L, show.x = TRUE) {
  stable <- cluster.stability[cluster.stability$protein == protein, , drop = FALSE]
  cluster.ids <- head(stable$reference_cluster_id, n.clusters)
  z <- matches[matches$protein == protein & matches$reference_cluster_id %in% cluster.ids, , drop = FALSE]
  keep <- if (is.na(threshold)) is.na(z$plddt_threshold) else z$plddt_threshold == threshold
  z <- z[keep, , drop = FALSE]
  trials <- sort(unique(z$n_sim))
  matrix.j <- matrix(NA_real_, nrow = length(cluster.ids), ncol = length(trials),
                     dimnames = list(cluster.ids, trials))
  for (i in seq_len(nrow(z))) matrix.j[z$reference_cluster_id[[i]], as.character(z$n_sim[[i]])] <- z$jaccard[[i]]
  matrix.j <- matrix.j[nrow(matrix.j):1L, , drop = FALSE]
  palette <- hcl.colors(101, "YlOrRd", rev = TRUE)
  image(seq_len(ncol(matrix.j)), seq_len(nrow(matrix.j)), t(matrix.j), zlim = c(0, 1),
        col = palette, axes = FALSE, xlab = if (show.x) "Null-simulation trials" else "", ylab = "",
        main = paste(protein, threshold_label(threshold)))
  axis(1, at = seq_len(ncol(matrix.j)), labels = colnames(matrix.j), cex.axis = 0.68)
  axis(2, at = seq_len(nrow(matrix.j)), labels = rownames(matrix.j), las = 2, cex.axis = 0.55)
  for (x in seq_len(ncol(matrix.j))) for (y in seq_len(nrow(matrix.j))) {
    value <- matrix.j[y, x]
    if (is.finite(value)) text(x, y, sprintf("%.2f", value), cex = 0.43,
                               col = if (value >= 0.75 || value <= 0.15) "white" else "black")
  }
  box()
}

plot_lines_by_threshold <- function(data, column, ylab, ylim = NULL, add.alpha = FALSE) {
  keys <- unique(threshold_key(data$plddt_threshold))
  if (is.null(ylim)) ylim <- range(data[[column]], finite = TRUE)
  if (!all(is.finite(ylim)) || diff(ylim) == 0) ylim <- c(ylim[[1L]] - 0.05, ylim[[1L]] + 0.05)
  plot(range(data$n_sim), ylim, type = "n", xlab = "Null-simulation trials", ylab = ylab)
  grid(col = "grey90")
  for (key in keys) {
    keep <- threshold_key(data$plddt_threshold) == key
    z <- data[keep, , drop = FALSE]; z <- z[order(z$n_sim), ]
    col <- get_threshold_color(z$plddt_threshold[[1L]])
    lines(z$n_sim, z[[column]], type = "b", pch = 19, lwd = 2, col = col)
  }
  if (add.alpha) abline(h = 0.05, lty = 2, col = colors$grey)
  legend("topright", legend = vapply(keys, function(key) threshold_label(if (key == "NA") NA else as.numeric(key)), character(1L)),
         col = vapply(keys, function(key) get_threshold_color(if (key == "NA") NA else as.numeric(key)), character(1L)),
         lty = 1, pch = 19, cex = 0.66, bty = "n")
}

grDevices::cairo_pdf(report.file, width = 11, height = 8.5, onefile = TRUE,
                     family = "sans")
on.exit(dev.off(), add = TRUE)

# 1. Title and executive conclusion -------------------------------------------------
page_setup("GRIN3D mutation-hotspot prototype evaluation",
           "TP53, SUZ12, and EZH2 | null-trial and AlphaFold pLDDT sensitivity")
text(0.06, 0.83, "Question", adj = c(0, 1), font = 2, cex = 1.05, col = colors$navy)
y <- wrapped_text(
  "Does the mutation prototype recover stable clusters in a dense positive-control protein and in proteins with sparser mutation patterns when null-simulation trials and structural-confidence filtering are varied?",
  0.06, 0.79, width = 102, cex = 0.9, line.height = 0.038
)
rect(0.055, y - 0.29, 0.945, y - 0.02, col = "#EDF4FA", border = colors$light.blue, lwd = 1.2)
text(0.075, y - 0.05, "Defensible conclusion", adj = c(0, 1), font = 2, cex = 1.02, col = colors$navy)
wrapped_text(
  "Observed cluster membership is deterministic and fully stable across 100-1,000 null trials at a fixed pLDDT threshold. Empirical p-values and significance calls are not uniformly converged by 1,000 trials. pLDDT filtering has a larger effect because it removes structural coordinates and mutation events, producing a different analyzed cohort. The prototype is suitable for cluster generation and exploratory sensitivity analysis, but final inference requires more null trials and independent-seed checks.",
  0.075, y - 0.09, width = 125, cex = 0.78, line.height = 0.029
)
facts <- c(
  paste(nrow(config), "completed configurations; trial counts:", paste(seq(100, 1000, 100), collapse = ", ")),
  "pLDDT thresholds: 0, 70, and 90 for TP53, SUZ12, and EZH2",
  paste(sum(validation$passed), "of", nrow(validation), "automated validation checks passed"),
  "Reference for each condition: its 1,000-trial result; membership reference across pLDDT: unfiltered 1,000-trial result"
)
bullet_list(facts, y = y - 0.34, width = 125, cex = 0.72, line.height = 0.027, gap = 0.009)

# 2. Study design -------------------------------------------------------------------
page_setup("Study design and comparison logic")
input.summary <- do.call(rbind, lapply(c("TP53", "SUZ12", "EZH2"), function(protein) {
  lesions <- read.csv(configured_file(protein, "lesion_file"), stringsAsFactors = FALSE)
  coordinates <- read.csv(configured_file(protein, "coordinate_file"), stringsAsFactors = FALSE)
  data.frame(
    Protein = protein, Events = nrow(lesions), Subjects = length(unique(lesions$ID)),
    `Protein residues` = nrow(coordinates), `pLDDT available` = ifelse("plddt" %in% names(coordinates), "Yes", "No"),
    check.names = FALSE
  )
}))
text_table(input.summary, y = 0.82, widths = c(0.16, 0.14, 0.14, 0.20, 0.18), cex = 0.78)
text(0.06, 0.56, "Varied settings", adj = c(0, 1), font = 2, cex = 1.0, col = colors$navy)
bullet_list(c(
  "Null simulations: 100 through 1,000 in increments of 100.",
  "Minimum pLDDT: 0 (unfiltered), 70 (confident), and 90 (very high confidence)."
), x = 0.07, y = 0.51, width = 78)
text(0.53, 0.56, "Held fixed", adj = c(0, 1), font = 2, cex = 1.0, col = colors$navy)
bullet_list(c(
  "Minimum support: 2 subjects and 2 events.",
  "Random seed 20260916; 50/50 calibration/evaluation split.",
  "Complete event mapping, within-subject overlap avoidance, and uniform positional null."
), x = 0.54, y = 0.51, width = 68)
rect(0.055, 0.10, 0.945, 0.22, col = "#FFF6DF", border = colors$gold)
wrapped_text(
  "Interpretation boundary: trial count changes the simulated null and empirical p-values, not the observed clustering tree. pLDDT changes the eligible structure and may exclude mutation events; pLDDT comparisons therefore do not contain identical mutation cohorts.",
  0.075, 0.195, width = 112, cex = 0.84, line.height = 0.034
)

# 3. Mutation landscapes -------------------------------------------------------------
new_plot_page("Mutation landscapes", "Mutation starts per residue; protein length makes density differences visible",
              mfrow = c(3, 1), mar = c(3.2, 4.2, 2.2, 3.2))
for (protein in c("TP53", "SUZ12", "EZH2")) {
  lesions <- read.csv(configured_file(protein, "lesion_file"), stringsAsFactors = FALSE)
  coordinates <- read.csv(configured_file(protein, "coordinate_file"), stringsAsFactors = FALSE)
  counts <- table(lesions$start)
  positions <- as.integer(names(counts)); heights <- as.integer(counts)
  plot(c(1, max(coordinates$residue)), c(0, max(heights) * 1.18), type = "n",
       xlab = "Residue", ylab = "Mutation events", main = protein)
  segments(positions, 0, positions, heights, col = colors$blue, lwd = 1.3)
  points(positions, heights, pch = 16, col = colors$red, cex = 0.7)
  if ("plddt" %in% names(coordinates)) {
    scaled <- coordinates$plddt / 100 * max(heights) * 0.20
    lines(coordinates$residue, scaled, col = adjustcolor(colors$green, alpha.f = 0.55), lwd = 1)
    legend("topright", c("mutation count", "pLDDT scaled to lower 20%"),
           col = c(colors$blue, colors$green), lty = 1, pch = c(16, NA), cex = 0.62, bty = "n")
  }
}
finish_plot_page()

# 4. Global Jaccard summary ----------------------------------------------------------
new_plot_page("Cluster-membership stability", "Mean best-match Jaccard against each protein's unfiltered 1,000-trial reference",
              mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1))
for (protein in c("TP53", "SUZ12", "EZH2")) {
  z <- run.metrics[run.metrics$protein == protein, , drop = FALSE]
  plot_lines_by_threshold(z, "mean_best_jaccard", "Mean best-match Jaccard", c(0, 1))
  title(main = protein)
  abline(h = c(0.5, 0.8), lty = c(3, 2), col = c(colors$grey, colors$green))
}
finish_plot_page()

# 5-7. Cluster-level Jaccard atlases -------------------------------------------------
new_plot_page("TP53 cluster-level Jaccard atlas", "Rows are the 15 highest-priority reference clusters; values are residue-set Jaccard",
              mfrow = c(3, 1), mar = c(3.5, 7, 2.6, 1))
for (threshold in c(0, 70, 90)) plot_jaccard_panel("TP53", threshold, n.clusters = 15L)
finish_plot_page()

new_plot_page("SUZ12 cluster-level Jaccard atlas", "All trial counts in each pLDDT panel; scores are printed in each cell",
              mfrow = c(3, 1), mar = c(3.5, 7, 2.6, 1))
for (threshold in c(0, 70, 90)) plot_jaccard_panel("SUZ12", threshold, n.clusters = 15L)
finish_plot_page()

new_plot_page("EZH2 cluster-level Jaccard atlas", "All trial counts in each pLDDT panel; scores are printed in each cell",
              mfrow = c(3, 1), mar = c(3.5, 7, 2.6, 1))
for (threshold in c(0, 70, 90)) plot_jaccard_panel("EZH2", threshold, n.clusters = 15L)
finish_plot_page()

# 8. Spearman convergence ------------------------------------------------------------
new_plot_page("Empirical p-value rank convergence", "Spearman correlation with the 1,000-trial result at the same pLDDT",
              mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1))
for (protein in c("TP53", "SUZ12", "EZH2")) {
  z <- convergence[convergence$protein == protein, , drop = FALSE]
  plot_lines_by_threshold(z, "pvalue_spearman", "Spearman correlation", c(0, 1))
  title(main = protein)
  abline(h = 0.9, lty = 2, col = colors$green)
}
finish_plot_page()

# 9. P-value error and significance calls -------------------------------------------
new_plot_page("Monte Carlo convergence", "P-value error decreases with trials, but corrected significance can change discretely",
              mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (protein in c("TP53", "SUZ12", "EZH2")) {
  z <- convergence[convergence$protein == protein, , drop = FALSE]
  plot_lines_by_threshold(z, "mean_absolute_p_difference", "Mean |p - p1000|")
  title(main = paste(protein, "p-value error"))
}
for (protein in c("TP53", "SUZ12", "EZH2")) {
  z <- convergence[convergence$protein == protein, , drop = FALSE]
  ylim <- c(0, max(1, z$n_significant))
  plot_lines_by_threshold(z, "n_significant", "Significant clusters", ylim)
  title(main = paste(protein, "significance calls"))
}
finish_plot_page()

# 10. pLDDT retention and membership ------------------------------------------------
new_plot_page("Effect of structural-confidence filtering", "Values at 1,000 trials; event loss and membership change are shown together",
              mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
for (protein in c("TP53", "SUZ12", "EZH2")) {
  z <- run.metrics[run.metrics$protein == protein & run.metrics$n_sim == 1000, , drop = FALSE]
  z <- z[order(z$plddt_threshold), ]
  event.retention <- z$mapped_events / z$mapped_events[z$plddt_threshold == 0]
  mat <- cbind(
    `Coordinates retained` = z$coordinate_retention,
    `Mutation events retained` = event.retention,
    `Mean best-match Jaccard` = z$mean_best_jaccard
  )
  matplot(z$plddt_threshold, mat, type = "b", pch = c(16, 17, 15), lty = 1,
          lwd = 2, col = c(colors$navy, colors$green, colors$red), ylim = c(0, 1),
          xlab = "Minimum pLDDT", ylab = "Fraction / similarity", main = protein)
  grid(col = "grey90")
  legend("bottomleft", colnames(mat), col = c(colors$navy, colors$green, colors$red),
         pch = c(16, 17, 15), lty = 1, cex = 0.70, bty = "n")
}
finish_plot_page()

# 11. Numeric results at maximum trial depth ----------------------------------------
page_setup("Results at 1,000 null-simulation trials",
           "Jaccard is measured against each protein's unfiltered 1,000-trial cluster landscape")
maximum.results <- run.metrics[run.metrics$n_sim == 1000, , drop = FALSE]
maximum.results <- maximum.results[order(match(maximum.results$protein, c("TP53", "SUZ12", "EZH2")),
                                         maximum.results$plddt_threshold), ]
display.maximum <- data.frame(
  Protein = maximum.results$protein,
  pLDDT = maximum.results$plddt_threshold,
  Residues = paste0(maximum.results$eligible_residues, "/", maximum.results$total_coordinate_residues),
  Events = paste0(maximum.results$mapped_events, "/",
                  maximum.results$mapped_events + maximum.results$excluded_events),
  Subjects = maximum.results$mapped_subjects,
  Clusters = maximum.results$n_clusters,
  Significant = maximum.results$n_significant,
  `Mean Jaccard` = sprintf("%.3f", maximum.results$mean_best_jaccard),
  check.names = FALSE
)
text_table(display.maximum, y = 0.82,
           widths = c(0.12, 0.09, 0.13, 0.12, 0.11, 0.11, 0.13, 0.14),
           row.height = 0.049, cex = 0.69)
rect(0.055, 0.12, 0.945, 0.23, col = "#EDF4FA", border = colors$light.blue)
wrapped_text(
  "Residue and event denominators are the unfiltered input totals. A filtered row is a different analyzed cohort, so changes in cluster or significance counts must be interpreted alongside retention and Jaccard—not as evidence that filtering itself strengthens association.",
  0.075, 0.205, width = 118, cex = 0.80, line.height = 0.031
)

# 12. Significance-set convergence --------------------------------------------------
page_setup("When do significance calls stabilize?",
           "Earliest trial count whose exact significant-cluster set matches 1,000 trials and remains matched")
display.earliest <- earliest
display.earliest$plddt_threshold <- ifelse(is.na(display.earliest$plddt_threshold), "Unavailable", display.earliest$plddt_threshold)
names(display.earliest)[names(display.earliest) == "protein"] <- "Protein"
names(display.earliest)[names(display.earliest) == "plddt_threshold"] <- "pLDDT"
names(display.earliest)[names(display.earliest) == "earliest_persistent_significance_match"] <- "Earliest match"
names(display.earliest)[names(display.earliest) == "reference_n_significant"] <- "Significant at 1,000"
display.earliest <- display.earliest[, c("Protein", "pLDDT", "Earliest match", "Significant at 1,000")]
text_table(display.earliest, y = 0.82, widths = c(0.16, 0.18, 0.20, 0.24), row.height = 0.047, cex = 0.72)
rect(0.055, 0.10, 0.945, 0.22, col = "#FFF6DF", border = colors$gold)
wrapped_text(
  "A match first reached at exactly 1,000 trials is not evidence of convergence because there is no later trial count to confirm persistence. Conditions with zero significant clusters have trivial agreement on an empty set and should not be described as positive convergence.",
  0.075, 0.195, width = 112, cex = 0.80, line.height = 0.032
)

# 13. Evidence synthesis -------------------------------------------------------------
page_setup("Evidence synthesis")
stability_status <- function(protein) {
  z <- earliest[earliest$protein == protein, ]
  z <- z[order(z$plddt_threshold), ]
  paste(z$earliest_persistent_significance_match, collapse = "/")
}
jaccard_status <- function(protein) {
  z <- run.metrics[run.metrics$protein == protein & run.metrics$n_sim == 1000 &
                     run.metrics$plddt_threshold %in% c(70, 90), ]
  z <- z[order(z$plddt_threshold), ]
  paste(sprintf("%.3f", z$mean_best_jaccard), collapse = " / ")
}
evidence.table <- data.frame(
  Protein = c("TP53", "SUZ12", "EZH2"),
  `Pattern represented` = c("Dense positive control", "Sparse long-protein pattern", "Sparse long-protein pattern"),
  `Trial membership` = c("Exact within pLDDT", "Exact within pLDDT", "Exact within pLDDT"),
  `pLDDT robustness` = paste("Mean J:", vapply(c("TP53", "SUZ12", "EZH2"), jaccard_status, character(1L))),
  `Statistical status` = paste("Earliest sets:", vapply(c("TP53", "SUZ12", "EZH2"), stability_status, character(1L))),
  check.names = FALSE
)
text_table(evidence.table, y = 0.82, widths = c(0.11, 0.23, 0.18, 0.23, 0.20),
           cex = 0.63, wrap.widths = c(10, 24, 20, 24, 24), line.height = 0.021)
text(0.06, 0.53, "What the evidence supports", adj = c(0, 1), font = 2, cex = 1.0, col = colors$navy)
bullet_list(c(
  "Cluster construction is reproducible with respect to null-trial count.",
  "Increasing trial count refines empirical p-values and can change corrected significance without changing cluster residues.",
  "pLDDT sensitivity is biologically and statistically consequential because filtering changes both structural coverage and the retained mutation cohort.",
  "TP53 provides a useful positive-control pattern; SUZ12 and EZH2 expose uncertainty that is less apparent in a dense hotspot example."
), y = 0.49, width = 112, cex = 0.80)

# 14. Recommendation ----------------------------------------------------------------
page_setup("Defensible recommendation")
text(0.06, 0.84, "Recommended use of the prototype", adj = c(0, 1), font = 2, cex = 1.05, col = colors$navy)
y <- bullet_list(c(
  "Use the observed cluster tree and residue memberships as reproducible candidate generation; report membership independently from p-values.",
  "Use pLDDT >= 70 as a sensitivity analysis rather than silently replacing the unfiltered analysis. Always report how many coordinates, events, and subjects were removed.",
  "Treat pLDDT >= 90 as a stringent stress test. It can remove a substantial fraction of coordinates and mutation events, materially changing the estimand.",
  "Do not make final significance claims from the 100-1,000 trial grid. Extend borderline conditions beyond 1,000 trials, ideally to the prototype's final-inference target, and repeat with independent seeds to quantify Monte Carlo variability.",
  "For each reported hotspot, show its residue set, Jaccard/containment under sensitivity settings, event retention, empirical p-value trajectory, and biological annotations as a separate evidence layer."
), y = 0.79, width = 112, cex = 0.82, line.height = 0.032, gap = 0.014)
rect(0.055, y - 0.13, 0.945, y - 0.01, col = "#EDF4FA", border = colors$light.blue)
wrapped_text(
  "Bottom line: the prototype is ready for documented exploratory hotspot discovery and sensitivity evaluation. It is not yet justified to label every corrected significance call as stable, especially where the significant set first appears only at the 1,000-trial endpoint or where pLDDT filtering removes a large fraction of the input events.",
  0.075, y - 0.035, width = 112, cex = 0.86, line.height = 0.034, font = 2
)

# 15. Reproducibility and provenance ------------------------------------------------
page_setup("Reproducibility and provenance")
text(0.06, 0.84, "Analysis artifacts", adj = c(0, 1), font = 2, cex = 1.0, col = colors$navy)
bullet_list(c(
  "TP53 structural confidence: AlphaFold DB AF-P04637-F1 model v6 (UniProt P04637); pLDDT read from the PDB B-factor field.",
  "01_run_sensitivity.R: declares and runs the trial/pLDDT grid.",
  "02_summarize_visualize.R: computes Jaccard, Dice, containment, overlap transitions, and run summaries.",
  "03_validate_analysis.R: verifies complete outputs, fixed settings, and within-threshold membership invariance.",
  "04_analyze_trial_convergence.R: evaluates p-value and significance convergence against 1,000 trials.",
  "05_build_pdf_report.R: builds this report from saved CSV and RDS results without rerunning simulations."
), y = 0.79, width = 125, cex = 0.72, line.height = 0.027, gap = 0.010)
text(0.06, 0.44, "Validation", adj = c(0, 1), font = 2, cex = 1.0, col = colors$navy)
validation.summary <- data.frame(
  Item = c("Configured runs", "Result objects", "Run manifests", "Validation checks", "Trial levels"),
  Result = c(nrow(config), sum(file.exists(file.path(config$results_dir, "GRIN3D_mutation_hotspot_results.rds"))),
             sum(file.exists(file.path(config$results_dir, "stability_run_manifest.csv"))),
             paste0(sum(validation$passed), "/", nrow(validation), " passed"),
             paste(sort(unique(config$n_sim)), collapse = ", ")),
  stringsAsFactors = FALSE
)
text_table(validation.summary, y = 0.39, widths = c(0.25, 0.62), row.height = 0.047, cex = 0.72)
wrapped_text(
  paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), ". Source directory: analysis/mutation_plddt_trials/. Statistical significance and biological annotation are intentionally maintained as separate evidence layers."),
  0.06, 0.09, width = 125, cex = 0.65, line.height = 0.026, col = colors$grey
)

dev.off()
on.exit(NULL, add = FALSE)
message("PDF report written to: ", report.file)
