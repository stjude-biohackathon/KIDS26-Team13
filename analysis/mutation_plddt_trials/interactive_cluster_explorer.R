build_cluster_explorer <- function(
    pvalue_clusters,
    run_metrics,
    template_file,
    max_clusters = 40L) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The jsonlite R package is required to build the cluster explorer.")
  }

  required <- c(
    "protein", "plddt_threshold", "n_sim", "candidate_key", "cluster_id",
    "residues", "n_subjects", "n_events", "p_any_joint", "significant_any"
  )
  missing <- setdiff(required, names(pvalue_clusters))
  if (length(missing)) {
    stop("Missing cluster columns: ", paste(missing, collapse = ", "))
  }

  reference <- pvalue_clusters[pvalue_clusters$n_sim == max(pvalue_clusters$n_sim), ]
  groups <- split(
    reference,
    interaction(reference$protein, reference$plddt_threshold, drop = TRUE)
  )
  selected <- do.call(rbind, lapply(groups, function(x) {
    residue_min <- vapply(
      strsplit(x$residues, ";", fixed = TRUE),
      function(z) min(as.integer(z)),
      integer(1L)
    )
    x <- x[order(x$p_any_joint, -x$n_subjects, residue_min), ]
    head(x, max_clusters)
  }))

  selected_keys <- unique(selected[c("protein", "plddt_threshold", "candidate_key")])
  history <- merge(
    pvalue_clusters,
    selected_keys,
    by = c("protein", "plddt_threshold", "candidate_key"),
    all = FALSE,
    sort = FALSE
  )

  compact <- function(x) {
    data.frame(
      protein = x$protein,
      threshold = x$plddt_threshold,
      trials = x$n_sim,
      key = x$candidate_key,
      cluster = x$cluster_id,
      residues = x$residues,
      subjects = x$n_subjects,
      events = x$n_events,
      p = x$p_any_joint,
      significant = x$significant_any,
      stringsAsFactors = FALSE
    )
  }

  lengths <- run_metrics[
    run_metrics$n_sim == max(run_metrics$n_sim) &
      run_metrics$plddt_threshold == min(run_metrics$plddt_threshold),
    c("protein", "total_coordinate_residues")
  ]
  lengths <- lengths[!duplicated(lengths$protein), ]
  names(lengths) <- c("protein", "length")

  payload <- list(
    history = compact(history),
    reference = compact(reference),
    lengths = lengths
  )
  payload_json <- jsonlite::toJSON(
    payload,
    dataframe = "rows",
    auto_unbox = TRUE,
    na = "null",
    digits = 8
  )

  template <- paste(readLines(template_file, warn = FALSE), collapse = "\n")
  sub("__CLUSTER_DATA__", payload_json, template, fixed = TRUE)
}
