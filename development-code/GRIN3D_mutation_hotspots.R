
# GRIN3D mutation-hotspot prototype scratch code
#
# ----------------------
# Required input columns
# ----------------------
# Mutation-lesion file:
#   ID, start, end
#   A unique event-ID column such as id is recommended. If it is omitted, the
#   code creates a unique event ID for each input row.
#
# AlphaFold C-alpha coordinate file:
#   residue, x, y, z
#   plddt and opportunity are optional. Set confidence = "plddt" to apply the
#   min.confidence filter. Set opportunity to a residue-weight column to use a
#   nonuniform positional null.
#   this file can be generated for selected proteins using GRIN3D_exon_to_protein_download_alphafold.R script
#
# This is an initial single-protein proof of concept. Without an opportunity
# column, each event is placed uniformly among structurally eligible windows.
# Supplying opportunity weights produces a weighted positional null instead.
# It describes how likely a mutation would be observed there under the 
# background model, before looking for hotspots.
#
# Column names can differ from these defaults. Supply the corresponding names
# through lesion.columns and coordinate.columns in the main runner.
#
# ---------------------------------
# Main analyses
# ----------------------------
# 1. Read mutation events with subject IDs, event IDs and amino-acid start and
#    end positions. Each mutation interval is retained as one event.
#
# 2. Read one AlphaFold C-alpha coordinate record per amino acid. Optionally,
#    use pLDDT to exclude unreliable residues. Residue-specific opportunity
#    weights can account for differences in background mutability;
#    without these weights, all eligible residues are sampled equally.
#
# 3. Map each mutation event to the protein structure. Events without the
#    required structural coordinates are identified and, when complete mapping
#    is required, excluded from the structural analysis.
#
# 4. Build separate 1D and 3D trees from the distinct altered residues:
#    - 1D tree: uses amino-acid sequence distance.
#    - 3D tree: uses C-alpha distance in the folded protein.
#    Each tree starts with one residue per branch. At each step, merge the two
#    branches with the smallest maximum pairwise distance in the merged group.
#    A merge may join two single residues, add a residue to a branch, or join
#    two larger branches.
#
# 5. Treat each single-residue leaf and merged branch as a candidate cluster.
#    If the same residue set appears in both trees, keep it once and record
#    both tree sources. Test candidates affecting at least the required number
#    of subjects and events. No fixed sequence or 3D distance cutoff is used.
#
# 6. Measure compactness for every candidate:
#    - 1D diameter: highest residue position minus lowest residue position.
#    - 3D diameter: largest pairwise C-alpha distance among its residues.
#
# 7. In each simulation, relocate all retained mutation events to eligible
#    protein positions. Keep each event's subject, ID and amino-acid interval
#    length; this also preserves the number of events per subject.
#
# 8. Rebuild both trees after each relocation. For each observed number of
#    affected subjects, record the smallest 1D and 3D diameters found anywhere
#    in the simulated tree search.
#
# 9. Calculate empirical size-specific p-values for compactness at the same subject count,
#    then correct for searching across subject counts and both trees. Use
#    separate subsets of the simulations for size calibration and search
#    correction.
# 9. Calculate size-specific p-values (p_1d_size, p_3d_size). For example,
#    if an observed cluster contains 24 events from 19 subjects and has a
#    3D diameter of 33 Å, ask how often a complete simulated tree search finds
#    a cluster anywhere in the protein with exactly 19 affected subjects and
#    a 3D diameter <= 33 Å. The simulated cluster need not contain 24 events.
#    Then correct for searching across subject counts within each tree
#    (p_1d_protein, p_3d_protein) and across both trees
#    (p_1d_joint, p_3d_joint, p_any_joint). Separate simulation subsets are
#    used for size calibration and search correction.
#
# 10. Save cluster, cluster-member, size-summary and event-mapping tables,
#     plus an RDS file containing the complete analysis.
#
# -----------------
# Expected output
# ---------------
# 1. mutation_hotspot_clusters.csv
#    Main hotspot-results table. Each row represents a unique residue cluster
#    identified by the 1D tree, the 3D tree or both trees.
#
#    The table reports:
#    - cluster and member-residue identifiers;
#    - originating tree or trees;
#    - affected subjects, events and residues;
#    - 1D sequence diameter and 3D structural diameter;
#    - size-specific, protein-wide and joint empirical p-values;
#    - significance indicators and hotspot classification.
#
# 2. mutation_hotspot_cluster_members.csv
#    Lists the mutation events and subjects contributing to each cluster,
#    including their original amino-acid intervals and the affected residues
#    contained within the cluster.
#
# 3. mutation_hotspot_size_summary.csv
#    Reports the most compact observed 1D and 3D clusters for each number of
#    affected subjects, together with their joint p-values.
#
# 4. mutation_event_mapping_summary.csv
#    Shows whether each input event was successfully mapped to structurally
#    eligible residues and identifies events excluded because of incomplete
#    structural mapping.
#
# 5. GRIN3D_mutation_hotspot_results.rds
#    Contains the complete analysis object, including input settings, mapped and
#    excluded events, cluster tables, null distributions and calibration results.
#
# ==============================================================================
# BIOHACKATHON PRIORITIES
# ==============================================================================
#
# 1. Validate the statistics
#    - Confirm that the null preserves subject IDs, events per subject and
#      mutation interval lengths.
#    - Check size-specific, protein-wide and joint 1D-3D p-values.
#
# 2. Test cluster robustness
#    - Compare distance definitions and tree-linkage rules.
#    - Benchmark the lesion-tree model against the original subject-kNN model.
#
# 3. Address structural uncertainty
#    - Assess regions using pLDDT, PAE and IDR annotations.
#    - Make 3D claims only where the structure is reliable.
#    - Continue to report 1D evidence in uncertain regions.
#
# 4. Add biological context (lower priority)
#    - Annotate domains, motifs, catalytic sites, PTMs, protein-interaction
#      interfaces and drug-binding sites.
#    - Keep annotations separate from statistical significance.
#
# 5. Test contrasting mutation patterns (high priority)
#    - Use TP53 as a dense-mutation example.
#    - Test proteins with sparse mutations that may converge in 3D (EX: SUZ12, EZH2).
#
# ==============================================================================
# 1. INPUT AND VALIDATION HELPERS
# ==============================================================================

# Display the default required and optional columns for both input files.
grin3d_mutation_input_spec <- function() {
  data.frame(
    input_file = c("Mutation lesions", "AlphaFold C-alpha coordinates"),
    required_columns = c("ID, start, end", "residue, x, y, z"),
    optional_columns = c("unique event ID, such as id", "plddt, opportunity"),
    stringsAsFactors = FALSE
  )
}

# Validate the names supplied in a column-mapping list before reading values.
validate_mutation_column_map <- function(
    columns,
    required.keys,
    optional.keys = character(),
    map.name) {
  if (!is.list(columns) || is.null(names(columns))) {
    stop(map.name, " must be a named list.")
  }
  missing.keys <- setdiff(required.keys, names(columns))
  if (length(missing.keys)) {
    stop(
      map.name, " is missing required mapping(s): ",
      paste(missing.keys, collapse = ", ")
    )
  }
  unknown.keys <- setdiff(names(columns), c(required.keys, optional.keys))
  if (length(unknown.keys)) {
    warning(
      map.name, " contains unused mapping(s): ",
      paste(unknown.keys, collapse = ", ")
    )
  }
  empty.required <- required.keys[vapply(
    columns[required.keys],
    function(value) {
      is.null(value) || length(value) != 1L || is.na(value) || !nzchar(value)
    },
    logical(1L)
  )]
  if (length(empty.required)) {
    stop(
      map.name, " has empty required mapping(s): ",
      paste(empty.required, collapse = ", ")
    )
  }
  invisible(TRUE)
}

# Validate scalar analysis settings before starting a long simulation run.
validate_mutation_analysis_settings <- function(
    min.subjects,
    min.events,
    alpha,
    n.sim,
    calibration.fraction,
    progress.every) {
  if (length(min.subjects) != 1L || is.na(min.subjects) || min.subjects < 1L) {
    stop("min.subjects must be one positive integer.")
  }
  if (length(min.events) != 1L || is.na(min.events) || min.events < 1L) {
    stop("min.events must be one positive integer.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1.")
  }
  if (length(n.sim) != 1L || is.na(n.sim) || n.sim < 100L) {
    stop("n.sim must be an integer of at least 100.")
  }
  if (
    length(calibration.fraction) != 1L ||
      !is.finite(calibration.fraction) ||
      calibration.fraction <= 0 || calibration.fraction >= 1
  ) {
    stop("calibration.fraction must be between 0 and 1.")
  }
  n.calibration <- floor(n.sim * calibration.fraction)
  if (n.calibration < 50L || n.sim - n.calibration < 50L) {
    stop(
      "n.sim and calibration.fraction must provide at least 50 calibration ",
      "and 50 evaluation simulations."
    )
  }
  if (
    length(progress.every) != 1L || is.na(progress.every) ||
      progress.every < 0L
  ) {
    stop("progress.every must be a nonnegative integer.")
  }
  invisible(TRUE)
}

# Read a CSV, TSV, TXT, RDS or RData input and return it as a data frame.
read_input_file <- function(file, object.name = NULL) {
  if (!file.exists(file)) {
    stop("Input file not found: ", file)
  }

  extension <- tolower(tools::file_ext(file))

  x <- switch(
    extension,
    csv = utils::read.csv(
      file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    tsv = utils::read.delim(
      file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    txt = utils::read.delim(
      file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    rds = readRDS(file),
    rdata = {
      env <- new.env(parent = emptyenv())
      loaded <- load(file, envir = env)

      if (is.null(object.name)) {
        if (length(loaded) != 1L) {
          stop(
            basename(file), " contains ", length(loaded),
            " objects. Set the corresponding '*.object' option to one of: ",
            paste(loaded, collapse = ", ")
          )
        }
        object.name <- loaded[[1L]]
      }

      if (!object.name %in% loaded) {
        stop(
          "Object '", object.name, "' was not found in ", basename(file),
          ". Available objects: ", paste(loaded, collapse = ", ")
        )
      }

      env[[object.name]]
    },
    stop(
      "Unsupported input format for ", basename(file),
      ". Supported extensions: csv, tsv, txt, rds, RData."
    )
  )

  if (!is.data.frame(x)) {
    stop(basename(file), " must contain a data.frame-like object.")
  }

  as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
}

# Confirm that all columns needed for a calculation are present.
require_columns <- function(x, columns, object.name) {
  columns <- columns[!is.na(columns) & nzchar(columns)]
  missing.columns <- setdiff(columns, names(x))

  if (length(missing.columns) > 0L) {
    stop(
      object.name, " is missing required column(s): ",
      paste(missing.columns, collapse = ", ")
    )
  }

  invisible(TRUE)
}

# Convert one input column to integer values and stop if conversion fails.
as_integer_column <- function(x, column, object.name) {
  value <- suppressWarnings(as.integer(x[[column]]))

  if (anyNA(value)) {
    stop(object.name, "$", column, " must contain nonmissing integers.")
  }

  value
}

# Convert one input column to finite numeric values and stop on invalid values.
as_numeric_column <- function(x, column, object.name) {
  value <- suppressWarnings(as.numeric(x[[column]]))

  if (anyNA(value) || any(!is.finite(value))) {
    stop(object.name, "$", column, " must contain finite numeric values.")
  }

  value
}

# Collapse unique values into a stable semicolon-separated label.
collapse_values <- function(x) {
  paste(sort(unique(as.character(x))), collapse = ";")
}


# ==============================================================================
# 2. PREPARE ALPHAFOLD COORDINATES AND MUTATION EVENTS
# ==============================================================================

# Validate C-alpha coordinates, apply optional pLDDT filtering and prepare
# residue-level opportunity weights for the positional null.
prepare_coordinates <- function(x, columns, min.confidence = 70) {
  require_columns(
    x,
    unlist(columns[c("residue", "x", "y", "z")], use.names = FALSE),
    "coordinate.data"
  )

  coordinates <- data.frame(
    residue = as_integer_column(x, columns$residue, "coordinate.data"),
    x = as_numeric_column(x, columns$x, "coordinate.data"),
    y = as_numeric_column(x, columns$y, "coordinate.data"),
    z = as_numeric_column(x, columns$z, "coordinate.data"),
    stringsAsFactors = FALSE
  )

  if (any(coordinates$residue < 1L)) {
    stop("Amino-acid residue numbers must be positive integers.")
  }

  if (anyDuplicated(coordinates$residue)) {
    duplicated.residues <- unique(
      coordinates$residue[duplicated(coordinates$residue)]
    )
    stop(
      "coordinate.data contains duplicate residue coordinates for: ",
      paste(head(duplicated.residues, 20L), collapse = ", ")
    )
  }

  if (!is.null(columns$confidence)) {
    require_columns(x, columns$confidence, "coordinate.data")
    coordinates$confidence <- as_numeric_column(
      x,
      columns$confidence,
      "coordinate.data"
    )
    coordinates <- coordinates[
      coordinates$confidence >= min.confidence,
      ,
      drop = FALSE
    ]
  } else {
    coordinates$confidence <- NA_real_
  }

  if (!is.null(columns$opportunity)) {
    require_columns(x, columns$opportunity, "coordinate.data")
    opportunity <- as_numeric_column(
      x,
      columns$opportunity,
      "coordinate.data"
    )

    if (any(opportunity < 0)) {
      stop("Residue opportunity weights must be nonnegative.")
    }

    coordinates$opportunity <- opportunity[
      match(coordinates$residue, as_integer_column(
        x,
        columns$residue,
        "coordinate.data"
      ))
    ]
  } else {
    coordinates$opportunity <- 1
  }

  if (nrow(coordinates) == 0L) {
    stop("No structurally eligible residues remain after coordinate filtering.")
  }

  if (sum(coordinates$opportunity) <= 0) {
    stop("At least one structurally eligible residue must have positive opportunity.")
  }

  coordinates <- coordinates[order(coordinates$residue), , drop = FALSE]
  rownames(coordinates) <- NULL
  coordinates
}

# Standardize mutation records and create event IDs when none are supplied.
prepare_events <- function(x, columns) {
  required <- unlist(columns[c("subject", "start", "end")], use.names = FALSE)
  if (!is.null(columns$event)) {
    required <- c(required, columns$event)
  }
  require_columns(x, required, "lesion.data")

  subject.id <- as.character(x[[columns$subject]])
  if (anyNA(subject.id) || any(!nzchar(subject.id))) {
    stop("lesion.data$", columns$subject, " contains a missing or empty ID.")
  }

  event.id <- if (is.null(columns$event)) {
    sprintf("event_%06d", seq_len(nrow(x)))
  } else {
    as.character(x[[columns$event]])
  }

  if (anyNA(event.id) || any(!nzchar(event.id))) {
    stop("The mutation event ID contains a missing or empty value.")
  }
  if (anyDuplicated(event.id)) {
    stop("Mutation event IDs must be unique. Duplicate event IDs were found.")
  }

  loc.start <- as_integer_column(x, columns$start, "lesion.data")
  loc.end <- as_integer_column(x, columns$end, "lesion.data")

  if (any(loc.start < 1L) || any(loc.end < loc.start)) {
    stop(
      "Mutation amino-acid intervals must satisfy 1 <= loc.start <= loc.end."
    )
  }

  data.frame(
    event_id = event.id,
    subject_id = subject.id,
    loc_start = loc.start,
    loc_end = loc.end,
    span_aa = loc.end - loc.start + 1L,
    input_row = seq_len(nrow(x)),
    stringsAsFactors = FALSE
  )
}

# Expand each mutation interval to one row per affected amino-acid residue.
expand_events <- function(events) {
  residue.list <- Map(seq.int, events$loc_start, events$loc_end)
  event.index <- rep(seq_len(nrow(events)), lengths(residue.list))

  data.frame(
    event_id = events$event_id[event.index],
    subject_id = events$subject_id[event.index],
    loc_start = events$loc_start[event.index],
    loc_end = events$loc_end[event.index],
    span_aa = events$span_aa[event.index],
    residue = unlist(residue.list, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

# Match mutation residues to eligible structural coordinates and separate
# retained events from incompletely mapped events.
map_events_to_structure <- function(
    events,
    coordinates,
    require.complete.mapping = TRUE) {
  expanded <- expand_events(events)
  expanded$coordinate_row <- match(expanded$residue, coordinates$residue)

  mapped.summary <- stats::aggregate(
    !is.na(expanded$coordinate_row),
    by = list(event_id = expanded$event_id),
    FUN = sum
  )
  names(mapped.summary)[[2L]] <- "mapped_residues"

  expected.summary <- stats::aggregate(
    expanded$residue,
    by = list(event_id = expanded$event_id),
    FUN = length
  )
  names(expected.summary)[[2L]] <- "expected_residues"

  mapping.summary <- merge(
    mapped.summary,
    expected.summary,
    by = "event_id",
    all = TRUE,
    sort = FALSE
  )
  mapping.summary$complete_mapping <- (
    mapping.summary$mapped_residues == mapping.summary$expected_residues
  )

  event.summary <- merge(
    events,
    mapping.summary,
    by = "event_id",
    all.x = TRUE,
    sort = FALSE
  )

  keep.event <- if (require.complete.mapping) {
    event.summary$complete_mapping
  } else {
    event.summary$mapped_residues > 0L
  }

  excluded.events <- event.summary[!keep.event, , drop = FALSE]
  retained.events <- event.summary[keep.event, names(events), drop = FALSE]

  mapping <- expanded[
    expanded$event_id %in% retained.events$event_id &
      !is.na(expanded$coordinate_row),
    ,
    drop = FALSE
  ]
  mapping$coordinate_row <- NULL

  mapping <- merge(
    mapping,
    coordinates,
    by = "residue",
    all.x = TRUE,
    sort = FALSE
  )
  mapping <- mapping[
    order(mapping$event_id, mapping$residue),
    ,
    drop = FALSE
  ]
  rownames(mapping) <- NULL

  list(
    events = retained.events,
    mapping = mapping,
    event_mapping_summary = event.summary,
    excluded_events = excluded.events
  )
}


# ==============================================================================
# 3. COMPLETE-LINKAGE CANDIDATE CLUSTERS
# ==============================================================================

# Recover the residue sets represented by every internal branch of an hclust tree.
extract_hclust_sets <- function(hc) {
  n.leaves <- length(hc$order)
  if (n.leaves < 2L) {
    return(list())
  }

  node.members <- vector("list", n.leaves - 1L)

  # Translate an hclust leaf or internal-node reference into leaf indices.
  resolve_node <- function(node) {
    if (node < 0L) {
      return(-node)
    }
    node.members[[node]]
  }

  for (i in seq_len(nrow(hc$merge))) {
    left <- resolve_node(hc$merge[i, 1L])
    right <- resolve_node(hc$merge[i, 2L])
    node.members[[i]] <- sort(unique(c(left, right)))
  }

  node.members
}

# Construct one complete-linkage tree and return its leaves and merged branches.
build_tree_sets <- function(residues, distance.matrix) {
  residues <- sort(unique(as.integer(residues)))

  leaf.sets <- lapply(residues, function(z) z)
  if (length(residues) < 2L) {
    return(leaf.sets)
  }

  hc <- stats::hclust(stats::as.dist(distance.matrix), method = "complete")
  index.sets <- extract_hclust_sets(hc)
  internal.sets <- lapply(index.sets, function(i) residues[i])

  c(leaf.sets, internal.sets)
}

# Calculate sequence span and maximum pairwise C-alpha distance for a residue set.
cluster_diameters <- function(residue.set, coordinates) {
  residue.set <- sort(unique(as.integer(residue.set)))

  diameter.1d <- if (length(residue.set) == 1L) {
    0
  } else {
    max(residue.set) - min(residue.set)
  }

  coordinate.rows <- match(residue.set, coordinates$residue)
  xyz <- as.matrix(coordinates[coordinate.rows, c("x", "y", "z")])

  diameter.3d <- if (nrow(xyz) == 1L) {
    0
  } else {
    max(as.matrix(stats::dist(xyz)))
  }

  c(diameter_1d = diameter.1d, diameter_3d = diameter.3d)
}

# Generate and deduplicate candidate clusters from the 1D and 3D trees, then
# retain clusters meeting the subject and event thresholds.
build_candidate_clusters <- function(
    mapping,
    coordinates,
    min.subjects = 2L,
    min.events = 2L) {
  altered.residues <- sort(unique(mapping$residue))
  coordinate.rows <- match(altered.residues, coordinates$residue)
  xyz <- as.matrix(coordinates[coordinate.rows, c("x", "y", "z")])

  distance.1d <- abs(outer(altered.residues, altered.residues, "-"))
  distance.3d <- if (length(altered.residues) == 1L) {
    matrix(0, nrow = 1L, ncol = 1L)
  } else {
    as.matrix(stats::dist(xyz))
  }

  sets.1d <- build_tree_sets(altered.residues, distance.1d)
  sets.3d <- build_tree_sets(altered.residues, distance.3d)

  source.map <- new.env(parent = emptyenv())
  residue.map <- new.env(parent = emptyenv())

  # Register candidate sets while tracking whether they came from 1D, 3D or both.
  register_sets <- function(sets, source) {
    for (residue.set in sets) {
      key <- paste(sort(unique(residue.set)), collapse = ",")
      previous <- if (exists(key, envir = source.map, inherits = FALSE)) {
        get(key, envir = source.map, inherits = FALSE)
      } else {
        character()
      }
      assign(key, sort(unique(c(previous, source))), envir = source.map)
      assign(key, sort(unique(residue.set)), envir = residue.map)
    }
  }

  register_sets(sets.1d, "1D")
  register_sets(sets.3d, "3D")

  keys <- ls(source.map, all.names = TRUE)
  rows <- vector("list", length(keys))

  for (i in seq_along(keys)) {
    key <- keys[[i]]
    residue.set <- get(key, envir = residue.map, inherits = FALSE)
    cluster.mapping <- mapping[
      mapping$residue %in% residue.set,
      ,
      drop = FALSE
    ]
    diameter <- cluster_diameters(residue.set, coordinates)

    rows[[i]] <- data.frame(
      candidate_key = paste(
        residue.set,
        collapse = ";"
      ),
      tree_sources = paste(
        get(key, envir = source.map, inherits = FALSE),
        collapse = ";"
      ),
      residues = paste(residue.set, collapse = ";"),
      residue_min = min(residue.set),
      residue_max = max(residue.set),
      n_residues = length(residue.set),
      n_events = length(unique(cluster.mapping$event_id)),
      n_subjects = length(unique(cluster.mapping$subject_id)),
      diameter_1d = unname(diameter[["diameter_1d"]]),
      diameter_3d = unname(diameter[["diameter_3d"]]),
      stringsAsFactors = FALSE
    )
  }

  candidates <- do.call(rbind, rows)
  candidates <- candidates[
    candidates$n_subjects >= min.subjects &
      candidates$n_events >= min.events,
    ,
    drop = FALSE
  ]

  if (nrow(candidates) == 0L) {
    return(candidates)
  }

  candidates <- candidates[
    order(
      candidates$n_subjects,
      candidates$diameter_3d,
      candidates$diameter_1d,
      candidates$candidate_key
    ),
    ,
    drop = FALSE
  ]
  rownames(candidates) <- NULL
  candidates
}


# ==============================================================================
# 4. SUBJECT- AND EVENT-PRESERVING MUTATION NULL
# ==============================================================================

# Enumerate every structurally eligible placement for each observed event span
# and assign optional opportunity-based sampling weights.
build_legal_windows <- function(events, coordinates) {
  eligible.residues <- coordinates$residue
  eligible.lookup <- setNames(rep(TRUE, length(eligible.residues)), eligible.residues)
  spans <- sort(unique(events$span_aa))
  windows <- vector("list", length(spans))
  names(windows) <- as.character(spans)

  for (span in spans) {
    possible.starts <- eligible.residues[
      eligible.residues + span - 1L <= max(eligible.residues)
    ]

    valid <- vapply(
      possible.starts,
      function(start) {
        residues <- seq.int(start, length.out = span)
        all(as.character(residues) %in% names(eligible.lookup))
      },
      logical(1L)
    )

    possible.starts <- possible.starts[valid]
    if (length(possible.starts) == 0L) {
      stop(
        "No structurally eligible amino-acid window can accommodate an ",
        "observed mutation span of ", span, " residue(s)."
      )
    }

    window.weights <- vapply(
      possible.starts,
      function(start) {
        residues <- seq.int(start, length.out = span)
        mean(coordinates$opportunity[match(residues, coordinates$residue)])
      },
      numeric(1L)
    )

    if (sum(window.weights) <= 0) {
      stop(
        "All opportunity weights are zero for mutation span ", span, "."
      )
    }

    windows[[as.character(span)]] <- list(
      starts = possible.starts,
      weights = window.weights
    )
  }

  windows
}

# Sample one value while handling a single available choice safely.
sample_one <- function(x, probability = NULL) {
  if (length(x) == 1L) {
    return(x[[1L]])
  }
  sample(x, size = 1L, prob = probability)
}

# Relocate all events while preserving event IDs, subjects, per-subject burden
# and amino-acid interval lengths.
simulate_event_mapping <- function(
    events,
    coordinates,
    legal.windows,
    avoid.within.subject.overlap = TRUE) {
  simulated <- vector("list", nrow(events))
  used.by.subject <- new.env(parent = emptyenv())

  for (i in seq_len(nrow(events))) {
    event <- events[i, , drop = FALSE]
    window.data <- legal.windows[[as.character(event$span_aa)]]
    starts <- window.data$starts
    weights <- window.data$weights

    subject <- event$subject_id
    used <- if (exists(subject, envir = used.by.subject, inherits = FALSE)) {
      get(subject, envir = used.by.subject, inherits = FALSE)
    } else {
      integer()
    }

    if (avoid.within.subject.overlap && length(used) > 0L) {
      nonoverlap <- vapply(
        starts,
        function(start) {
          !any(seq.int(start, length.out = event$span_aa) %in% used)
        },
        logical(1L)
      )

      if (any(nonoverlap)) {
        starts <- starts[nonoverlap]
        weights <- weights[nonoverlap]
      }
    }

    start <- sample_one(starts, weights)
    residues <- seq.int(start, length.out = event$span_aa)
    assign(
      subject,
      sort(unique(c(used, residues))),
      envir = used.by.subject
    )

    coordinate.rows <- match(residues, coordinates$residue)
    simulated[[i]] <- data.frame(
      event_id = event$event_id,
      subject_id = subject,
      loc_start = start,
      loc_end = start + event$span_aa - 1L,
      span_aa = event$span_aa,
      residue = residues,
      x = coordinates$x[coordinate.rows],
      y = coordinates$y[coordinate.rows],
      z = coordinates$z[coordinate.rows],
      confidence = coordinates$confidence[coordinate.rows],
      opportunity = coordinates$opportunity[coordinate.rows],
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, simulated)
}

# Record the smallest candidate diameter found for each affected-subject count.
minimum_diameter_by_size <- function(candidates, sizes, metric) {
  result <- rep(Inf, length(sizes))
  names(result) <- as.character(sizes)

  if (nrow(candidates) == 0L) {
    return(result)
  }

  for (i in seq_along(sizes)) {
    value <- candidates[candidates$n_subjects == sizes[[i]], metric]
    if (length(value) > 0L) {
      result[[i]] <- min(value)
    }
  }

  result
}

# Repeat complete event relocation and rebuild both trees in every simulation.
run_null_simulations <- function(
    events,
    coordinates,
    sizes,
    n.sim,
    min.subjects,
    min.events,
    avoid.within.subject.overlap,
    random.seed,
    progress.every = 100L) {
  legal.windows <- build_legal_windows(events, coordinates)
  null.1d <- matrix(
    Inf,
    nrow = n.sim,
    ncol = length(sizes),
    dimnames = list(NULL, as.character(sizes))
  )
  null.3d <- null.1d

  set.seed(random.seed)

  for (b in seq_len(n.sim)) {
    simulated.mapping <- simulate_event_mapping(
      events = events,
      coordinates = coordinates,
      legal.windows = legal.windows,
      avoid.within.subject.overlap = avoid.within.subject.overlap
    )

    simulated.candidates <- build_candidate_clusters(
      mapping = simulated.mapping,
      coordinates = coordinates,
      min.subjects = min.subjects,
      min.events = min.events
    )

    null.1d[b, ] <- minimum_diameter_by_size(
      simulated.candidates,
      sizes,
      "diameter_1d"
    )
    null.3d[b, ] <- minimum_diameter_by_size(
      simulated.candidates,
      sizes,
      "diameter_3d"
    )

    if (
      progress.every > 0L &&
        (b %% progress.every == 0L || b == n.sim)
    ) {
      message("Completed ", b, " of ", n.sim, " null simulations.")
    }
  }

  list(
    diameter_1d = null.1d,
    diameter_3d = null.3d
  )
}


# ==============================================================================
# 5. SIZE-SPECIFIC AND JOINT NULL CALIBRATION
# ==============================================================================

# Calculate a lower-tail empirical p-value with the standard plus-one correction.
empirical_lower_p <- function(observed, null.values) {
  (1 + sum(null.values <= observed)) / (length(null.values) + 1)
}

# Calculate lower-tail empirical p-values against a sorted null distribution.
empirical_p_from_reference <- function(values, sorted.reference) {
  counts <- findInterval(values, sorted.reference)
  (1 + counts) / (length(sorted.reference) + 1)
}

# Split simulations into calibration and evaluation sets, then calculate
# size-specific, protein-wide, joint and omnibus p-values for every candidate.
calibrate_candidates <- function(
    candidates,
    null.results,
    calibration.fraction = 0.50,
    alpha = 0.05) {
  n.sim <- nrow(null.results$diameter_1d)
  n.calibration <- floor(n.sim * calibration.fraction)

  if (n.calibration < 50L || n.sim - n.calibration < 50L) {
    stop(
      "Use enough simulations to retain at least 50 calibration and 50 ",
      "joint-evaluation replicates."
    )
  }

  calibration.rows <- seq_len(n.calibration)
  evaluation.rows <- seq.int(n.calibration + 1L, n.sim)
  sizes <- as.integer(colnames(null.results$diameter_1d))

  calibration.1d <- null.results$diameter_1d[calibration.rows, , drop = FALSE]
  calibration.3d <- null.results$diameter_3d[calibration.rows, , drop = FALSE]
  evaluation.1d <- null.results$diameter_1d[evaluation.rows, , drop = FALSE]
  evaluation.3d <- null.results$diameter_3d[evaluation.rows, , drop = FALSE]

  evaluation.p.1d <- matrix(
    NA_real_,
    nrow = nrow(evaluation.1d),
    ncol = ncol(evaluation.1d)
  )
  evaluation.p.3d <- evaluation.p.1d

  sorted.1d <- vector("list", length(sizes))
  sorted.3d <- vector("list", length(sizes))

  for (j in seq_along(sizes)) {
    sorted.1d[[j]] <- sort(calibration.1d[, j])
    sorted.3d[[j]] <- sort(calibration.3d[, j])

    evaluation.p.1d[, j] <- empirical_p_from_reference(
      evaluation.1d[, j],
      sorted.1d[[j]]
    )
    evaluation.p.3d[, j] <- empirical_p_from_reference(
      evaluation.3d[, j],
      sorted.3d[[j]]
    )
  }

  evaluation.min.1d <- apply(evaluation.p.1d, 1L, min)
  evaluation.min.3d <- apply(evaluation.p.3d, 1L, min)
  evaluation.min.joint <- pmin(evaluation.min.1d, evaluation.min.3d)

  candidates$p_1d_size <- NA_real_
  candidates$p_3d_size <- NA_real_
  candidates$p_1d_protein <- NA_real_
  candidates$p_3d_protein <- NA_real_
  candidates$p_1d_joint <- NA_real_
  candidates$p_3d_joint <- NA_real_
  candidates$p_any_joint <- NA_real_

  for (i in seq_len(nrow(candidates))) {
    j <- match(candidates$n_subjects[[i]], sizes)

    p.1d.size <- empirical_lower_p(
      candidates$diameter_1d[[i]],
      calibration.1d[, j]
    )
    p.3d.size <- empirical_lower_p(
      candidates$diameter_3d[[i]],
      calibration.3d[, j]
    )

    candidates$p_1d_size[[i]] <- p.1d.size
    candidates$p_3d_size[[i]] <- p.3d.size
    candidates$p_1d_protein[[i]] <- (
      1 + sum(evaluation.min.1d <= p.1d.size)
    ) / (length(evaluation.min.1d) + 1)
    candidates$p_3d_protein[[i]] <- (
      1 + sum(evaluation.min.3d <= p.3d.size)
    ) / (length(evaluation.min.3d) + 1)
    candidates$p_1d_joint[[i]] <- (
      1 + sum(evaluation.min.joint <= p.1d.size)
    ) / (length(evaluation.min.joint) + 1)
    candidates$p_3d_joint[[i]] <- (
      1 + sum(evaluation.min.joint <= p.3d.size)
    ) / (length(evaluation.min.joint) + 1)
    candidates$p_any_joint[[i]] <- (
      1 + sum(evaluation.min.joint <= min(p.1d.size, p.3d.size))
    ) / (length(evaluation.min.joint) + 1)
  }

  significant.1d <- candidates$p_1d_joint <= alpha
  significant.3d <- candidates$p_3d_joint <= alpha

  candidates$significant_1d <- significant.1d
  candidates$significant_3d <- significant.3d
  candidates$significant_any <- candidates$p_any_joint <= alpha

  candidates$best_1d_for_size <- vapply(
    seq_len(nrow(candidates)),
    function(i) {
      same.size <- candidates$n_subjects == candidates$n_subjects[[i]]
      candidates$diameter_1d[[i]] == min(candidates$diameter_1d[same.size])
    },
    logical(1L)
  )
  candidates$best_3d_for_size <- vapply(
    seq_len(nrow(candidates)),
    function(i) {
      same.size <- candidates$n_subjects == candidates$n_subjects[[i]]
      candidates$diameter_3d[[i]] == min(candidates$diameter_3d[same.size])
    },
    logical(1L)
  )

  candidates$hotspot_class <- ifelse(
    significant.1d & significant.3d,
    "Sequence-local hotspot also compact in structure",
    ifelse(
      !significant.1d & significant.3d,
      "Nonlocal, 3D-specific hotspot",
      ifelse(
        significant.1d & !significant.3d,
        "Linear hotspot without independent 3D enrichment",
        "Weak clustering evidence"
      )
    )
  )

  candidates <- candidates[
    order(
      candidates$p_any_joint,
      candidates$p_3d_joint,
      candidates$p_1d_joint,
      candidates$diameter_3d,
      candidates$candidate_key
    ),
    ,
    drop = FALSE
  ]
  candidates$cluster_id <- sprintf("cluster_%04d", seq_len(nrow(candidates)))
  rownames(candidates) <- NULL

  preferred.columns <- c(
    "cluster_id", "candidate_key", "tree_sources", "hotspot_class",
    "significant_1d", "significant_3d", "significant_any",
    "best_1d_for_size", "best_3d_for_size",
    "n_subjects", "n_events", "n_residues", "residues",
    "residue_min", "residue_max", "diameter_1d", "diameter_3d",
    "p_1d_size", "p_3d_size", "p_1d_protein", "p_3d_protein",
    "p_1d_joint", "p_3d_joint", "p_any_joint"
  )

  list(
    clusters = candidates[, preferred.columns, drop = FALSE],
    evaluation_min_p = data.frame(
      min_p_1d = evaluation.min.1d,
      min_p_3d = evaluation.min.3d,
      min_p_joint = evaluation.min.joint
    ),
    n_calibration = n.calibration,
    n_evaluation = length(evaluation.rows)
  )
}


# ==============================================================================
# 6. OUTPUT TABLES
# ==============================================================================

# Expand the cluster table into one row per mutation event in each cluster.
build_cluster_member_table <- function(clusters, observed.mapping) {
  rows <- vector("list", nrow(clusters))

  for (i in seq_len(nrow(clusters))) {
    residues <- as.integer(strsplit(
      clusters$residues[[i]],
      ";",
      fixed = TRUE
    )[[1L]])

    member.mapping <- observed.mapping[
      observed.mapping$residue %in% residues,
      c(
        "event_id", "subject_id", "loc_start", "loc_end",
        "span_aa", "residue"
      ),
      drop = FALSE
    ]

    event.rows <- split(member.mapping, member.mapping$event_id)
    rows[[i]] <- do.call(
      rbind,
      lapply(event.rows, function(z) {
        data.frame(
          cluster_id = clusters$cluster_id[[i]],
          event_id = z$event_id[[1L]],
          subject_id = z$subject_id[[1L]],
          loc_start = z$loc_start[[1L]],
          loc_end = z$loc_end[[1L]],
          span_aa = z$span_aa[[1L]],
          affected_residues_in_cluster = collapse_values(z$residue),
          stringsAsFactors = FALSE
        )
      })
    )
  }

  members <- do.call(rbind, rows)
  members <- members[
    order(members$cluster_id, members$subject_id, members$event_id),
    ,
    drop = FALSE
  ]
  rownames(members) <- NULL
  members
}

# Summarize the most compact observed 1D and 3D cluster at each subject count.
build_size_summary <- function(clusters) {
  sizes <- sort(unique(clusters$n_subjects))
  rows <- lapply(sizes, function(size) {
    x <- clusters[clusters$n_subjects == size, , drop = FALSE]
    best.1d <- x[which.min(x$diameter_1d), , drop = FALSE]
    best.3d <- x[which.min(x$diameter_3d), , drop = FALSE]

    data.frame(
      n_subjects = size,
      best_1d_cluster = best.1d$cluster_id,
      minimum_1d_diameter = best.1d$diameter_1d,
      best_1d_joint_p = best.1d$p_1d_joint,
      best_3d_cluster = best.3d$cluster_id,
      minimum_3d_diameter = best.3d$diameter_3d,
      best_3d_joint_p = best.3d$p_3d_joint,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}


# ==============================================================================
# 7. USER-FACING ANALYSIS RUNNER
# ==============================================================================

# Run one complete protein analysis, write the result files and return the full
# result object invisibly for interactive use in R.
run_grin3d_mutation_hotspots <- function(
    lesion.file,
    coordinate.file,
    results.dir = "results_mutation_hotspots",
    protein = "protein_of_interest",
    lesion.columns = list(
      subject = "ID",
      event = NULL,
      start = "loc.start",
      end = "loc.end"
    ),
    coordinate.columns = list(
      residue = "residue",
      x = "x",
      y = "y",
      z = "z",
      confidence = NULL,
      opportunity = NULL
    ),
    lesion.object = NULL,
    coordinate.object = NULL,
    min.subjects = 2L,
    min.events = 2L,
    min.confidence = 70,
    require.complete.event.mapping = TRUE,
    avoid.within.subject.overlap = TRUE,
    alpha = 0.05,
    n.sim = 1000L,
    calibration.fraction = 0.50,
    random.seed = 20260828L,
    progress.every = 100L) {

  validate_mutation_column_map(
    lesion.columns,
    required.keys = c("subject", "start", "end"),
    optional.keys = "event",
    map.name = "lesion.columns"
  )
  validate_mutation_column_map(
    coordinate.columns,
    required.keys = c("residue", "x", "y", "z"),
    optional.keys = c("confidence", "opportunity"),
    map.name = "coordinate.columns"
  )
  validate_mutation_analysis_settings(
    min.subjects = min.subjects,
    min.events = min.events,
    alpha = alpha,
    n.sim = n.sim,
    calibration.fraction = calibration.fraction,
    progress.every = progress.every
  )
  if (
    !is.null(coordinate.columns$confidence) &&
      (length(min.confidence) != 1L || !is.finite(min.confidence))
  ) {
    stop("min.confidence must be one finite number when filtering is enabled.")
  }

  analysis.options <- list(
    protein = protein,
    min.subjects = min.subjects,
    min.events = min.events,
    min.confidence = min.confidence,
    require.complete.event.mapping = require.complete.event.mapping,
    avoid.within.subject.overlap = avoid.within.subject.overlap,
    alpha = alpha,
    n.sim = n.sim,
    calibration.fraction = calibration.fraction,
    random.seed = random.seed,
    progress.every = progress.every
  )

  dir.create(results.dir, recursive = TRUE, showWarnings = FALSE)

  message("Reading mutation and AlphaFold coordinate data.")
  lesion.data <- read_input_file(lesion.file, lesion.object)
  coordinate.data <- read_input_file(coordinate.file, coordinate.object)

  if (is.null(coordinate.columns$confidence)) {
    message(
      "No confidence column was selected; min.confidence will not be applied."
    )
  } else {
    message(
      "Filtering structural coordinates using ",
      coordinate.columns$confidence, " >= ", min.confidence, "."
    )
  }
  if (is.null(coordinate.columns$opportunity)) {
    message("Using a uniform positional mutation null.")
  } else {
    message(
      "Using residue weights from the ", coordinate.columns$opportunity,
      " column in the positional null."
    )
  }

  coordinates <- prepare_coordinates(
    coordinate.data,
    coordinate.columns,
    min.confidence = analysis.options$min.confidence
  )
  events <- prepare_events(lesion.data, lesion.columns)
  mapped <- map_events_to_structure(
    events,
    coordinates,
    require.complete.mapping = analysis.options$require.complete.event.mapping
  )

  if (nrow(mapped$events) == 0L) {
    stop("No mutation events remain after structural-coordinate mapping.")
  }

  if (nrow(mapped$excluded_events) > 0L) {
    warning(
      nrow(mapped$excluded_events),
      " mutation event(s) were excluded because their amino-acid interval was ",
      "not completely represented among structurally eligible residues."
    )
  }

  message("Building observed 1D and 3D lesion-residue trees.")
  observed.candidates <- build_candidate_clusters(
    mapping = mapped$mapping,
    coordinates = coordinates,
    min.subjects = analysis.options$min.subjects,
    min.events = analysis.options$min.events
  )

  if (nrow(observed.candidates) == 0L) {
    stop(
      "No candidate cluster satisfies min.subjects = ",
      analysis.options$min.subjects,
      " and min.events = ", analysis.options$min.events, "."
    )
  }

  cluster.sizes <- sort(unique(observed.candidates$n_subjects))

  message(
    "Running ", analysis.options$n.sim,
    " subject- and event-preserving null simulations."
  )
  null.results <- run_null_simulations(
    events = mapped$events,
    coordinates = coordinates,
    sizes = cluster.sizes,
    n.sim = analysis.options$n.sim,
    min.subjects = analysis.options$min.subjects,
    min.events = analysis.options$min.events,
    avoid.within.subject.overlap = analysis.options$avoid.within.subject.overlap,
    random.seed = analysis.options$random.seed,
    progress.every = analysis.options$progress.every
  )

  message("Calibrating size-specific and joint 1D/3D p-values.")
  calibrated <- calibrate_candidates(
    candidates = observed.candidates,
    null.results = null.results,
    calibration.fraction = analysis.options$calibration.fraction,
    alpha = analysis.options$alpha
  )

  cluster.table <- calibrated$clusters
  cluster.table$protein <- analysis.options$protein
  cluster.table <- cluster.table[, c(
    "protein",
    setdiff(names(cluster.table), "protein")
  )]

  cluster.member.table <- build_cluster_member_table(
    cluster.table,
    mapped$mapping
  )
  cluster.size.summary <- build_size_summary(cluster.table)

  result.object <- list(
    settings = list(
      lesion.file = lesion.file,
      coordinate.file = coordinate.file,
      results.dir = results.dir,
      lesion.columns = lesion.columns,
      coordinate.columns = coordinate.columns,
      confidence.filter.applied = !is.null(coordinate.columns$confidence),
      opportunity.weights.applied = !is.null(coordinate.columns$opportunity),
      analysis.options = analysis.options
    ),
    coordinates = coordinates,
    mapped_events = mapped$events,
    excluded_events = mapped$excluded_events,
    event_mapping_summary = mapped$event_mapping_summary,
    event_residue_mapping = mapped$mapping,
    clusters = cluster.table,
    cluster_members = cluster.member.table,
    size_summary = cluster.size.summary,
    null_diameters = null.results,
    null_min_p = calibrated$evaluation_min_p,
    n_calibration = calibrated$n_calibration,
    n_evaluation = calibrated$n_evaluation
  )

  utils::write.csv(
    result.object$clusters,
    file.path(results.dir, "mutation_hotspot_clusters.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    result.object$cluster_members,
    file.path(results.dir, "mutation_hotspot_cluster_members.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    result.object$size_summary,
    file.path(results.dir, "mutation_hotspot_size_summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    result.object$event_mapping_summary,
    file.path(results.dir, "mutation_event_mapping_summary.csv"),
    row.names = FALSE
  )
  saveRDS(
    result.object,
    file.path(results.dir, "GRIN3D_mutation_hotspot_results.rds")
  )

  message("Analysis complete. Results were written to: ", results.dir)
  message(
    "Use n.sim = 100000L for final inference after validating this prototype ",
    "on known positive and null examples."
  )

  invisible(result.object)
}
