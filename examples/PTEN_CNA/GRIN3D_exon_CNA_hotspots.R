# GRIN3D genomic CNA hotspot prototype scratch code
#
# Required input columns
# ----------------------
# CNA lesion file:
#   ID, chrom, loc.start, loc.end, lsn.type
#   An optional event-ID column may also be supplied. If it is omitted, the
#   code creates a unique event ID for each input row.
#
# Coding-exon file for one selected transcript:
#   exon_id, exon_rank, chrom, strand, cds_genomic_start, cds_genomic_end,
#   aa_start, aa_end
#   gene_id, gene_name, transcript_id and protein_id are recommended metadata.
#   this file can be generated for selected proteins using GRIN3D_exon_to_protein_download_alphafold.R script
#
# AlphaFold C-alpha coordinate file:
#   residue, x, y, z
#   plddt is optional. To apply a confidence filter, set
#   coordinate.columns$confidence = "plddt". If confidence is NULL, the
#   min.confidence setting is not applied.
#   This file ca be generated for selected proteins using GRIN3D_exon_to_protein_download_alphafold.R script
#
# Column names can differ from these defaults. Supply the corresponding names
# through cna.columns, exon.columns and coordinate.columns in the main runner.
#
# Main analyses
# -------------
# 1. Map every intact CNA segment to the coding exons it overlaps.
#
# 2. Count events and unique subjects per exon and CNA type, separately for
#    all coding-overlap events (cover all exons) and the structurally localizable subset.
#
# 3. Evaluate CNA breadth at 50%, 75%, 90% and 100% coding-sequence
#    coverage. Conditional simulations test whether CNAs of the observed
#    genomic sizes cover unusually large fractions of the coding sequence.
#
# 4. Test whether individual coding exons are affected by more subjects or
#    events than expected under two separate recurrence models:
#    ALL_OVERLAP includes every CNA overlapping the exon, whereas
#    PARTIAL_LOCALIZABLE includes only partial CNAs that can localize evidence
#    within the protein. Each scope is calibrated with its corresponding null.
#
# 5. Build two complete-linkage exon trees:
#    - 1D tree: groups exons that are close in transcript order.
#    - 3D tree: calculate the distance between each 2 exons in the 3D space based on
#      using the mean of the five smallest reliable C-alpha distances between their
#      encoded residues. The number of residue pairs can be varied, for example from 5
#      to 10 or 15, in sensitivity analyses.
#
# 6. Test each candidate cluster against simulations that preserve each CNA as
#    one intact event, along with its subject, CNA type, exon span and approximate
#    coding-residue span. Run lesion-type-specific and combined CNA analyses.
#
# 7. For each merged branch, report how many new subjects and events it adds
#    beyond its most strongly supported smaller cluster. This helps distinguish
#    an informative larger hotspot from a strong core plus weakly related exons.
#
# 8. Check that each CNA has enough alternative exon placements for a reliable
#    simulation test. Flag CNAs with too few possible placements or little
#    variation across simulations.
#
# Coverage p-values condition on a CNA already intersecting this gene's coding
# sequence. They test coverage localization, not genome-wide gene recurrence tested by GRIN2.
# Whole-coding-sequence events are retained in coverage summaries but are not
# assigned an artificial 3D hotspot position.
#
# ---------------------------------------------
# Expected output
# ---------------------------------------------
# 1. CNA_event_mapping_summary.csv
#    One row per CNA, showing its coding coverage, affected exons, localization
#    class and eligibility for structural analysis.
#
# 2. CNA_event_exon_membership.csv
#    One row per CNA-exon overlap, including the number and fraction of coding
#    bases affected within that exon.
#
# 3. CNA_exon_counts_by_type.csv
#    Event and subject counts for each exon and CNA type. Counts are reported
#    for all overlapping CNAs and for the partial-localizable subset.
#
# 4. CNA_gene_coverage_results.csv
#    Subject and event counts, with conditional p-values, for CNAs covering at
#    least 50%, 75%, 90% or 100% of the coding sequence. These are gene-coverage
#    results, not exon-level structural-hotspot results.
#
# 5. CNA_exon_recurrence_results.csv
#    Exon-level recurrence results under two scopes:
#    ALL_OVERLAP includes every CNA affecting the exon.
#    PARTIAL_LOCALIZABLE includes only partial CNAs capable of localizing
#    evidence within the protein.
#
# 6. CNA_exon_structural_distances.csv
#    Pairwise 3D distances between coding exons. The primary distance is the
#    mean of the five smallest reliable cross-exon C-alpha distances.
#
# 7. CNA_structural_testability.csv
#    Indicates whether each CNA type has enough localizable subjects, events
#    and structurally eligible exons for 1D/3D hotspot testing.
#
# 8. CNA_exon_hotspot_clusters.csv
#    The main structural-hotspot table. It reports cluster exons, CNA type,
#    affected subjects and events, 1D and 3D diameters, empirical p-values and
#    hotspot classification.
#
#    p_1d_size:
#      Tests whether any simulated cluster with the same number of affected subjects
#      is as compact in exon order as the observed cluster.
#
#    p_3d_size:
#      Tests whether any simulated cluster with the same number of affected
#      subjects has a 3D diameter as small as the observed cluster.
#
#    p_1d_protein and p_3d_protein:
#      Correct the corresponding 1D or 3D result for searching across candidate
#      locations, tree branches and supported subject counts within the protein.
#
#    p_1d_joint and p_3d_joint:
#      Evaluate the observed 1D or 3D evidence against a null search that
#      considers both trees.
#
#    p_any_joint:
#      Tests whether the cluster is significant through either its 1D or 3D
#      evidence after calibr complete 1D-3D search correction.
#
#    p_hotspot_omnibus:
#      Primary cluster-level result. It summarizes the strongest 1D or 3D
#      evidence using the jointly calibrated null.
#
#    significant_1d, significant_3d and significant_any:
#      Indicate whether the corresponding corrected p-value is at or below the
#      selected alpha threshold.
#
# 9. CNA_exon_hotspot_cluster_members.csv
#    Lists the individual exons and CNA events contributing to each cluster.
#
# 10. CNA_partial_recurrence_permutation_diagnostics.csv
#     Reports whether each partial CNA has enough alternative exon placements
#     to support an informative simulation-based p-value.
#
# 11. GRIN3D_exon_CNA_results.rds
#     Contains the complete analysis object, including output tables, null
#     distributions, distance matrices, mappings and analysis settings.
#
# ----------------------------------
# Interpretation
# ---------------------------------
# 1D-only hotspot:
#   Recurrent CNAs affect consecutive or nearby exons, but the exons do not
#   show independent evidence of unusual structural proximity.
#
# 3D-specific hotspot:
#   CNA-affected exons may be separated in linear exon order but converge in
#   the folded protein.
#
# Jointly supported hotspot:
#   The cluster is compact in both exon order and 3D structure.
#
# Weak clustering evidence:
#   The cluster is not significant after correction for the complete tree search.
#
# Exon recurrence and structural clustering answer different questions.
# A recurrent individual exon is not necessarily part of a significant
# multi-exon 3D cluster.
#
# p_hotspot_omnibus is the primary corrected cluster-level p-value because it
# accounts for searching both the 1D and 3D trees. Size-specific p-values are
# useful supporting results but should not be interpreted as final significance.
#
# Whole-gene CNAs support gene-level disruption but do not localize a structural
# hotspot. Their absence from the cluster table is therefore intentional.
#
# ==============================================================================
# BIOHACKATHON PRIORITIES
# ==============================================================================
#
# 1. Define CNA localization
#    - Compare coding-coverage classes: <=50%, 50-75%, 75-<100% and 100%.
#    - Determine which partial CNAs qualify for exon-level 3D analysis.
#    - Retain broad and whole-gene CNAs as gene-level evidence.
#
# 2. Validate the statistical framework
#    - Confirm that simulations preserve subject, CNA type and intact event span.
#    - Verify exon-recurrence, 1D, 3D and omnibus p-value calculations.
#    - Check that each CNA has enough alternative placement windows.
#    - Use subject counts as primary evidence and event counts as supporting evidence.
#
# 3. Evaluate exon proximity
#    - Assess whether coding exons are the appropriate structural unit.
#    - Compare the mean of the 5 versus 10 closest cross-exon C-alpha distances.
#    - Select the primary distance metric before comparing statistical results.
#    - Use pLDDT to restrict 3D inference to reliable protein regions.
#
# 4. Test robustness
#    - Can use PTEN to evaluate multi-exon structural convergence.
#    - Use CDKN2A to examine genes with few coding exons.
#    - Test additional genes with different CNA and exon patterns.
#    - Compare cluster stability across different analysis settings (EX: distance
#      metrics and confidence levels).
#
# Goal
# ----
# Demonstrate feasibility, identify limitations and define the settings required
# for a validated exon-level CNA hotspot model.

# ==============================================================================
# 1. GENERAL HELPERS
# ==============================================================================

# Display the default required and optional columns for the three input files.
grin3d_cna_input_spec <- function() {
  data.frame(
    input_file = c(
      "CNA lesions",
      "Coding-exon map",
      "AlphaFold C-alpha coordinates"
    ),
    required_columns = c(
      "ID, chrom, loc.start, loc.end, lsn.type",
      paste(
        "exon_id, exon_rank, chrom, strand, cds_genomic_start,",
        "cds_genomic_end, aa_start, aa_end"
      ),
      "residue, x, y, z"
    ),
    optional_columns = c(
      "unique event ID",
      "gene_id, gene_name, transcript_id, protein_id",
      "plddt"
    ),
    stringsAsFactors = FALSE
  )
}

# Validate the names supplied in a column-mapping list before reading values.
validate_grin3d_column_map <- function(
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
  allowed.keys <- c(required.keys, optional.keys)
  unknown.keys <- setdiff(names(columns), allowed.keys)
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

# Validate hierarchical-clustering linkage method.
validate_linkage_method <- function(linkage.method) {
  allowed.methods <- c("complete", "average", "single", "ward.D2")
  if (
    length(linkage.method) != 1L ||
      !is.character(linkage.method) ||
      is.na(linkage.method) ||
      !(linkage.method %in% allowed.methods)
  ) {
    stop(
      "Invalid linkage.method '", paste(linkage.method, collapse = ", "),
      "'. Allowed methods are: ", paste(allowed.methods, collapse = ", "), "."
    )
  }
  invisible(TRUE)
}

# Read a CSV, TSV, TXT, RDS or RData input and return it as a data frame.
read_grin3d_input <- function(file, object.name = NULL) {
  if (!file.exists(file)) stop("Input file not found: ", file)

  extension <- tolower(tools::file_ext(file))
  x <- switch(
    extension,
    csv = utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE),
    tsv = utils::read.delim(file, stringsAsFactors = FALSE, check.names = FALSE),
    txt = utils::read.delim(file, stringsAsFactors = FALSE, check.names = FALSE),
    rds = readRDS(file),
    rdata = {
      env <- new.env(parent = emptyenv())
      loaded <- load(file, envir = env)
      if (is.null(object.name)) {
        if (length(loaded) != 1L) {
          stop(
            basename(file), " contains ", length(loaded),
            " objects. Supply object.name. Available objects: ",
            paste(loaded, collapse = ", ")
          )
        }
        object.name <- loaded[[1L]]
      }
      if (!object.name %in% loaded) {
        stop("Object '", object.name, "' was not found in ", basename(file))
      }
      env[[object.name]]
    },
    stop("Unsupported input extension: ", extension)
  )

  if (!is.data.frame(x)) stop(basename(file), " must contain a data.frame.")
  as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
}

# Confirm that all columns needed for a calculation are present.
require_grin3d_columns <- function(x, columns, object.name) {
  columns <- unique(columns[!is.na(columns) & nzchar(columns)])
  missing <- setdiff(columns, names(x))
  if (length(missing)) {
    stop(object.name, " is missing column(s): ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

# Convert one input column to integer values and stop if conversion fails.
as_grin3d_integer <- function(x, column, object.name) {
  out <- suppressWarnings(as.integer(x[[column]]))
  if (anyNA(out)) stop(object.name, "$", column, " must contain integers.")
  out
}

# Convert one input column to finite numeric values and stop on invalid values.
as_grin3d_numeric <- function(x, column, object.name) {
  out <- suppressWarnings(as.numeric(x[[column]]))
  if (anyNA(out) || any(!is.finite(out))) {
    stop(object.name, "$", column, " must contain finite numeric values.")
  }
  out
}

# Standardize chromosome labels so that, for example, chr10 and 10 match.
normalize_grin3d_chromosome <- function(x) {
  toupper(sub("^CHR", "", trimws(as.character(x)), ignore.case = TRUE))
}

# Collapse unique values into a stable semicolon-separated label.
collapse_grin3d_values <- function(x) {
  paste(sort(unique(x)), collapse = ";")
}

# Convert a lesion-type label into a safe suffix for output column names.
safe_grin3d_name <- function(x) {
  gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", toupper(x)))
}

# Combine a list of data-frame rows while preserving a defined empty schema.
rbind_grin3d_rows <- function(rows, empty.columns = character()) {
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (!length(rows)) {
    return(as.data.frame(
      setNames(replicate(length(empty.columns), logical(), simplify = FALSE),
               empty.columns),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Calculate a lower-tail empirical p-value with the standard plus-one correction.
cna_empirical_lower_p <- function(observed, null.values) {
  (1 + sum(null.values <= observed)) / (length(null.values) + 1)
}

# Calculate an upper-tail empirical p-value with the standard plus-one correction.
cna_empirical_upper_p <- function(observed, null.values) {
  (1 + sum(null.values >= observed)) / (length(null.values) + 1)
}

# Calculate lower-tail empirical p-values against a sorted null distribution.
cna_empirical_p_from_reference <- function(values, sorted.reference) {
  counts <- findInterval(values, sorted.reference)
  (1 + counts) / (length(sorted.reference) + 1)
}

# Calculate upper-tail empirical p-values for integer recurrence counts.
cna_empirical_upper_p_from_reference <- function(values, sorted.reference) {
  n <- length(sorted.reference)
  # Recurrence null statistics are integer counts. Subtracting 0.5 makes
  # findInterval count values strictly below each observed count, preserving
  # conservative treatment of ties in the upper tail.
  counts.ge <- n - findInterval(values - 0.5, sorted.reference)
  (1 + counts.ge) / (n + 1)
}


# ==============================================================================
# 2. COORDINATES AND CODING-EXON MAP
# ==============================================================================

# Validate C-alpha coordinates and optionally remove residues below the
# requested structural-confidence threshold.
prepare_cna_coordinates <- function(x, columns, min.confidence = 70) {
  required <- unlist(columns[c("residue", "x", "y", "z")], use.names = FALSE)
  require_grin3d_columns(x, required, "coordinate.data")

  coordinates <- data.frame(
    residue = as_grin3d_integer(x, columns$residue, "coordinate.data"),
    x = as_grin3d_numeric(x, columns$x, "coordinate.data"),
    y = as_grin3d_numeric(x, columns$y, "coordinate.data"),
    z = as_grin3d_numeric(x, columns$z, "coordinate.data"),
    stringsAsFactors = FALSE
  )
  if (any(coordinates$residue < 1L)) stop("Residue numbers must be positive.")
  if (anyDuplicated(coordinates$residue)) {
    stop("coordinate.data must have one row per residue.")
  }

  if (!is.null(columns$confidence)) {
    require_grin3d_columns(x, columns$confidence, "coordinate.data")
    coordinates$confidence <- as_grin3d_numeric(
      x, columns$confidence, "coordinate.data"
    )
    coordinates <- coordinates[
      coordinates$confidence >= min.confidence, , drop = FALSE
    ]
  } else {
    coordinates$confidence <- NA_real_
  }
  if (!nrow(coordinates)) stop("No structurally eligible residues remain.")
  coordinates <- coordinates[order(coordinates$residue), , drop = FALSE]
  rownames(coordinates) <- NULL
  coordinates
}

# Prepare the selected transcript's coding exons, map their residue ranges to
# available coordinates and determine which exons are structurally eligible.
prepare_coding_exons <- function(
    x,
    columns,
    coordinates,
    min.mapped.fraction = 0.50,
    min.mapped.residues = 1L) {
  required <- unlist(
    columns[c(
      "exon_id", "exon_rank", "chrom", "strand", "genomic_start",
      "genomic_end", "residue_start", "residue_end"
    )],
    use.names = FALSE
  )
  require_grin3d_columns(x, required, "exon.data")

  optional_value <- function(name, default = NA_character_) {
    column <- columns[[name]]
    if (is.null(column) || !column %in% names(x)) {
      rep(default, nrow(x))
    } else {
      x[[column]]
    }
  }

  exons <- data.frame(
    gene_id = as.character(optional_value("gene_id")),
    gene_name = as.character(optional_value("gene_name")),
    transcript_id = as.character(optional_value("transcript_id")),
    protein_id = as.character(optional_value("protein_id")),
    exon_id = as.character(x[[columns$exon_id]]),
    transcript_exon_rank = as_grin3d_integer(
      x, columns$exon_rank, "exon.data"
    ),
    chrom = normalize_grin3d_chromosome(x[[columns$chrom]]),
    strand = as.character(x[[columns$strand]]),
    cds_genomic_start = as_grin3d_integer(
      x, columns$genomic_start, "exon.data"
    ),
    cds_genomic_end = as_grin3d_integer(
      x, columns$genomic_end, "exon.data"
    ),
    residue_start = as_grin3d_integer(
      x, columns$residue_start, "exon.data"
    ),
    residue_end = as_grin3d_integer(
      x, columns$residue_end, "exon.data"
    ),
    stringsAsFactors = FALSE
  )

  if (anyNA(exons$exon_id) || any(!nzchar(exons$exon_id)) ||
      anyDuplicated(exons$exon_id)) {
    stop("Coding exon IDs must be nonempty and unique.")
  }
  if (length(unique(exons$chrom)) != 1L ||
      length(unique(exons$strand)) != 1L) {
    stop("exon.data must describe one transcript on one chromosome and strand.")
  }
  if (!all(exons$strand %in% c("+", "-"))) {
    stop("The transcript strand must be '+' or '-'.")
  }
  if (any(exons$cds_genomic_start > exons$cds_genomic_end) ||
      any(exons$residue_start > exons$residue_end) ||
      any(exons$residue_start < 1L)) {
    stop("Exon genomic and residue intervals must satisfy start <= end.")
  }

  exons <- exons[order(exons$transcript_exon_rank), , drop = FALSE]
  exons$exon_order <- seq_len(nrow(exons))
  exons$cds_length_nt <- exons$cds_genomic_end - exons$cds_genomic_start + 1L

  residue.rows <- lapply(seq_len(nrow(exons)), function(i) {
    residues <- seq.int(exons$residue_start[[i]], exons$residue_end[[i]])
    coordinate.rows <- match(residues, coordinates$residue)
    data.frame(
      exon_order = exons$exon_order[[i]],
      transcript_exon_rank = exons$transcript_exon_rank[[i]],
      exon_id = exons$exon_id[[i]],
      residue = residues,
      coordinate_row = coordinate.rows,
      structurally_eligible = !is.na(coordinate.rows),
      stringsAsFactors = FALSE
    )
  })
  exon.residue.map <- do.call(rbind, residue.rows)
  rownames(exon.residue.map) <- NULL

  exons$n_coding_residues <- vapply(
    exons$exon_order,
    function(i) sum(exon.residue.map$exon_order == i),
    integer(1L)
  )
  exons$n_structural_residues <- vapply(
    exons$exon_order,
    function(i) sum(
      exon.residue.map$exon_order == i &
        exon.residue.map$structurally_eligible
    ),
    integer(1L)
  )
  exons$mapped_fraction <- exons$n_structural_residues / exons$n_coding_residues
  exons$structurally_eligible <-
    exons$n_structural_residues >= min.mapped.residues &
    exons$mapped_fraction >= min.mapped.fraction
  exons$structural_eligibility_reason <- ifelse(
    exons$structurally_eligible,
    "eligible",
    ifelse(
      exons$n_structural_residues < min.mapped.residues,
      "fewer_than_minimum_mapped_residues",
      "mapped_fraction_below_threshold"
    )
  )

  list(exons = exons, exon_residue_map = exon.residue.map)
}


# ==============================================================================
# 3. GENOMIC CNA-TO-EXON MAPPING
# ==============================================================================

# Standardize raw CNA records and create event IDs when none are supplied.
prepare_genomic_cna_events <- function(x, columns) {
  required <- unlist(
    columns[c("subject", "chrom", "start", "end", "type")],
    use.names = FALSE
  )
  if (!is.null(columns$event)) required <- c(required, columns$event)
  require_grin3d_columns(x, required, "cna.data")

  subject <- as.character(x[[columns$subject]])
  if (anyNA(subject) || any(!nzchar(subject))) {
    stop("The CNA subject column contains missing or empty values.")
  }
  event <- if (is.null(columns$event)) {
    sprintf("CNA_event_%06d", seq_len(nrow(x)))
  } else {
    as.character(x[[columns$event]])
  }
  if (anyNA(event) || any(!nzchar(event)) || anyDuplicated(event)) {
    stop("CNA event IDs must be nonmissing and unique.")
  }

  out <- data.frame(
    event_id = event,
    subject_id = subject,
    cna_type = toupper(trimws(as.character(x[[columns$type]]))),
    chrom = normalize_grin3d_chromosome(x[[columns$chrom]]),
    genomic_start = as_grin3d_integer(x, columns$start, "cna.data"),
    genomic_end = as_grin3d_integer(x, columns$end, "cna.data"),
    input_row = seq_len(nrow(x)),
    stringsAsFactors = FALSE
  )
  if (anyNA(out$cna_type) || any(!nzchar(out$cna_type))) {
    stop("CNA types must be nonmissing and nonempty.")
  }
  if (any(out$genomic_start < 1L) ||
      any(out$genomic_end < out$genomic_start)) {
    stop("CNA intervals must satisfy 1 <= loc.start <= loc.end.")
  }
  out$width_bp <- out$genomic_end - out$genomic_start + 1L
  out
}

# Assign a coding-coverage resolution class to one CNA event.
classify_cna_localization <- function(fraction, tolerance = 1e-12) {
  if (!is.finite(fraction) || fraction <= 0) return("no_coding_overlap")
  if (fraction >= 1 - tolerance) return("whole_coding_sequence")
  if (fraction <= 0.50) return("high_resolution_le_50pct")
  if (fraction <= 0.75) return("moderate_resolution_50_75pct")
  if (fraction <= 0.90) return("low_resolution_75_90pct")
  "very_low_resolution_90_lt_100pct"
}

# Intersect every intact CNA with the selected transcript's coding exons and
# summarize its coding coverage, affected exons and structural eligibility.
map_genomic_cnas_to_exons <- function(
    events,
    exons,
    exon.residue.map,
    min.exon.coding.overlap = 0) {
  target.chrom <- unique(exons$chrom)
  total.cds.nt <- sum(exons$cds_length_nt)
  total.residues <- length(unique(exon.residue.map$residue))
  total.exons <- nrow(exons)
  mapping.rows <- list()
  map.index <- 0L

  for (i in seq_len(nrow(events))) {
    if (events$chrom[[i]] != target.chrom) next
    overlap.start <- pmax(events$genomic_start[[i]], exons$cds_genomic_start)
    overlap.end <- pmin(events$genomic_end[[i]], exons$cds_genomic_end)
    overlap.bp <- pmax(0L, overlap.end - overlap.start + 1L)
    hit <- which(overlap.bp > 0L)
    if (!length(hit)) next

    for (j in hit) {
      map.index <- map.index + 1L
      fraction <- overlap.bp[[j]] / exons$cds_length_nt[[j]]
      mapping.rows[[map.index]] <- data.frame(
        event_id = events$event_id[[i]],
        subject_id = events$subject_id[[i]],
        cna_type = events$cna_type[[i]],
        input_row = events$input_row[[i]],
        chrom = events$chrom[[i]],
        genomic_start = events$genomic_start[[i]],
        genomic_end = events$genomic_end[[i]],
        width_bp = events$width_bp[[i]],
        exon_order = exons$exon_order[[j]],
        transcript_exon_rank = exons$transcript_exon_rank[[j]],
        exon_id = exons$exon_id[[j]],
        exon_cds_start = exons$cds_genomic_start[[j]],
        exon_cds_end = exons$cds_genomic_end[[j]],
        overlap_start = overlap.start[[j]],
        overlap_end = overlap.end[[j]],
        coding_overlap_bp = overlap.bp[[j]],
        exon_coding_fraction = fraction,
        passes_exon_overlap = fraction > min.exon.coding.overlap,
        stringsAsFactors = FALSE
      )
    }
  }

  mapping.columns <- c(
    "event_id", "subject_id", "cna_type", "input_row", "chrom",
    "genomic_start", "genomic_end", "width_bp", "exon_order",
    "transcript_exon_rank", "exon_id", "exon_cds_start", "exon_cds_end",
    "overlap_start", "overlap_end", "coding_overlap_bp",
    "exon_coding_fraction", "passes_exon_overlap"
  )
  mapping <- rbind_grin3d_rows(mapping.rows, mapping.columns)
  mapped.ids <- unique(mapping$event_id)
  summary.rows <- vector("list", nrow(events))

  for (i in seq_len(nrow(events))) {
    event <- events[i, , drop = FALSE]
    z <- mapping[mapping$event_id == event$event_id, , drop = FALSE]
    if (!nrow(z)) {
      summary.rows[[i]] <- data.frame(
        event,
        coding_overlap_bp = 0L,
        coding_fraction = 0,
        n_affected_exons = 0L,
        affected_exon_orders = "",
        affected_transcript_exon_ranks = "",
        affected_exon_ids = "",
        exon_start = NA_integer_,
        exon_end = NA_integer_,
        n_coding_residues = 0L,
        residue_fraction_exon_proxy = 0,
        complete_structural_mapping = FALSE,
        structurally_localizable = FALSE,
        localization_class = "no_coding_overlap",
        covers_ge_50pct = FALSE,
        covers_ge_75pct = FALSE,
        covers_ge_90pct = FALSE,
        covers_100pct = FALSE,
        stringsAsFactors = FALSE
      )
      next
    }

    affected.orders <- sort(unique(z$exon_order[z$passes_exon_overlap]))
    if (!length(affected.orders)) affected.orders <- sort(unique(z$exon_order))
    affected.residues <- unique(exon.residue.map$residue[
      exon.residue.map$exon_order %in% affected.orders
    ])
    coding.bp <- sum(z$coding_overlap_bp)
    coding.fraction <- min(1, coding.bp / total.cds.nt)
    complete.structure <- all(
      exons$structurally_eligible[match(affected.orders, exons$exon_order)]
    )
    whole.coding <- coding.fraction >= 1 - 1e-12 ||
      length(affected.orders) == total.exons

    summary.rows[[i]] <- data.frame(
      event,
      coding_overlap_bp = coding.bp,
      coding_fraction = coding.fraction,
      n_affected_exons = length(affected.orders),
      affected_exon_orders = collapse_grin3d_values(affected.orders),
      affected_transcript_exon_ranks = collapse_grin3d_values(
        exons$transcript_exon_rank[match(affected.orders, exons$exon_order)]
      ),
      affected_exon_ids = collapse_grin3d_values(
        exons$exon_id[match(affected.orders, exons$exon_order)]
      ),
      exon_start = min(affected.orders),
      exon_end = max(affected.orders),
      n_coding_residues = length(affected.residues),
      residue_fraction_exon_proxy = length(affected.residues) / total.residues,
      complete_structural_mapping = complete.structure,
      structurally_localizable = !whole.coding && complete.structure,
      localization_class = classify_cna_localization(coding.fraction),
      covers_ge_50pct = coding.fraction >= 0.50,
      covers_ge_75pct = coding.fraction >= 0.75,
      covers_ge_90pct = coding.fraction >= 0.90,
      covers_100pct = coding.fraction >= 1 - 1e-12,
      stringsAsFactors = FALSE
    )
  }

  event.summary <- do.call(rbind, summary.rows)
  mapping <- merge(
    mapping,
    event.summary[, c(
      "event_id", "coding_fraction", "n_affected_exons",
      "localization_class", "complete_structural_mapping",
      "structurally_localizable"
    )],
    by = "event_id", all.x = TRUE, sort = FALSE
  )
  mapping <- mapping[
    order(mapping$input_row, mapping$exon_order), , drop = FALSE
  ]
  rownames(mapping) <- NULL
  rownames(event.summary) <- NULL

  list(
    events = event.summary[
      event.summary$event_id %in% mapped.ids, , drop = FALSE
    ],
    all_events = event.summary,
    unmapped_events = event.summary[
      !event.summary$event_id %in% mapped.ids, , drop = FALSE
    ],
    mapping = mapping,
    total_cds_nt = total.cds.nt,
    total_unique_residues = total.residues
  )
}

# Count unique CNA events and subjects per exon, lesion type and recurrence
# scope (all overlap versus structurally localizable events).
build_exon_type_counts <- function(
    exons,
    mapping,
    cna.types,
    structural.mapping = NULL) {
  out <- exons[, c(
    "gene_id", "gene_name", "transcript_id", "protein_id", "exon_id",
    "transcript_exon_rank", "exon_order", "chrom", "strand",
    "cds_genomic_start", "cds_genomic_end", "residue_start", "residue_end",
    "n_coding_residues", "n_structural_residues", "mapped_fraction",
    "structurally_eligible", "structural_eligibility_reason"
  ), drop = FALSE]

  for (type in c(cna.types, "ALL_CNA")) {
    safe <- safe_grin3d_name(type)
    z.all <- if (type == "ALL_CNA") {
      mapping
    } else {
      mapping[mapping$cna_type == type, , drop = FALSE]
    }
    z.structural <- if (is.null(structural.mapping)) {
      mapping[FALSE, , drop = FALSE]
    } else if (type == "ALL_CNA") {
      structural.mapping
    } else {
      structural.mapping[
        structural.mapping$cna_type == type, , drop = FALSE
      ]
    }

    count_one <- function(z) {
      event.values <- integer(nrow(out))
      subject.values <- integer(nrow(out))
      if (!nrow(z)) {
        return(list(events = event.values, subjects = subject.values))
      }
      event.count <- tapply(
        z$event_id, z$exon_order, function(v) length(unique(v))
      )
      subject.count <- tapply(
        z$subject_id, z$exon_order, function(v) length(unique(v))
      )
      event.values <- as.integer(
        event.count[match(as.character(out$exon_order), names(event.count))]
      )
      subject.values <- as.integer(
        subject.count[match(as.character(out$exon_order), names(subject.count))]
      )
      event.values[is.na(event.values)] <- 0L
      subject.values[is.na(subject.values)] <- 0L
      list(events = event.values, subjects = subject.values)
    }

    all.counts <- count_one(z.all)
    structural.counts <- count_one(z.structural)

    # Retain the original column names for backward compatibility.
    out[[paste0("n_events_1d_", safe)]] <- all.counts$events
    out[[paste0("n_subjects_1d_", safe)]] <- all.counts$subjects

    # Explicit names distinguish all coding-overlap events from the subset
    # that is eligible for the exon-structure analysis.
    out[[paste0("n_events_all_overlap_", safe)]] <- all.counts$events
    out[[paste0("n_subjects_all_overlap_", safe)]] <- all.counts$subjects
    out[[paste0("n_events_structural_localizable_", safe)]] <-
      structural.counts$events
    out[[paste0("n_subjects_structural_localizable_", safe)]] <-
      structural.counts$subjects
  }
  out
}


# ==============================================================================
# 4. CONDITIONAL GENE-COVERAGE NULL
# ==============================================================================

# Merge overlapping or adjacent integer intervals into nonoverlapping ranges.
merge_integer_intervals <- function(starts, ends) {
  order.index <- order(starts, ends)
  starts <- as.numeric(starts[order.index])
  ends <- as.numeric(ends[order.index])
  merged.start <- starts[[1L]]
  merged.end <- ends[[1L]]
  out.start <- numeric()
  out.end <- numeric()

  if (length(starts) > 1L) {
    for (i in 2:length(starts)) {
      if (starts[[i]] <= merged.end + 1) {
        merged.end <- max(merged.end, ends[[i]])
      } else {
        out.start <- c(out.start, merged.start)
        out.end <- c(out.end, merged.end)
        merged.start <- starts[[i]]
        merged.end <- ends[[i]]
      }
    }
  }
  data.frame(
    start = c(out.start, merged.start),
    end = c(out.end, merged.end)
  )
}

# Find every genomic start position where an interval of the observed width
# would overlap at least one coding exon of the selected transcript.
legal_cna_starts_conditioned_on_cds <- function(
    width.bp,
    exons,
    chromosome.length = NULL) {
  starts <- exons$cds_genomic_start - width.bp + 1
  ends <- exons$cds_genomic_end

  lower <- 1
  upper <- Inf
  if (!is.null(chromosome.length)) {
    if (length(chromosome.length) != 1L || !is.finite(chromosome.length) ||
        chromosome.length < 1) {
      stop("chromosome.length must be NULL or one positive finite number.")
    }
    if (width.bp > chromosome.length) {
      stop("A CNA event is longer than the supplied chromosome length.")
    }
    upper <- chromosome.length - width.bp + 1
  }

  starts <- pmax(starts, lower)
  ends <- pmin(ends, upper)
  keep <- starts <= ends
  if (!any(keep)) {
    stop("No legal interval placement can overlap the coding sequence.")
  }
  merge_integer_intervals(starts[keep], ends[keep])
}

# Draw one genomic coordinate uniformly from a set of integer intervals.
sample_integer_from_intervals <- function(intervals) {
  sizes <- intervals$end - intervals$start + 1
  total <- sum(sizes)
  offset <- floor(stats::runif(1L, min = 0, max = total))
  cumulative <- cumsum(sizes)
  selected <- which(offset < cumulative)[[1L]]
  previous <- if (selected == 1L) 0 else cumulative[[selected - 1L]]
  intervals$start[[selected]] + offset - previous
}

# Calculate the fraction of the transcript's coding bases covered by an interval.
coding_fraction_for_interval <- function(start, end, exons) {
  overlap <- pmax(
    0,
    pmin(end, exons$cds_genomic_end) -
      pmax(start, exons$cds_genomic_start) + 1
  )
  min(1, sum(overlap) / sum(exons$cds_length_nt))
}

# Summarize subjects and events reaching each coding-coverage threshold.
coverage_count_table <- function(
    events,
    fractions,
    cna.types,
    thresholds,
    include.combined.analysis = TRUE) {
  rows <- list()
  index <- 0L
  labels <- cna.types
  if (include.combined.analysis && length(cna.types) > 1L) {
    labels <- c(labels, "ALL_CNA")
  }
  for (type in labels) {
    type.keep <- if (type == "ALL_CNA") {
      rep(TRUE, nrow(events))
    } else {
      events$cna_type == type
    }
    for (threshold in thresholds) {
      index <- index + 1L
      keep <- type.keep & fractions + 1e-12 >= threshold
      rows[[index]] <- data.frame(
        analysis_type = type,
        coverage_threshold = threshold,
        n_events = length(unique(events$event_id[keep])),
        n_subjects = length(unique(events$subject_id[keep])),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out$row_key <- paste(
    out$analysis_type,
    sprintf("%.6f", out$coverage_threshold),
    sep = "|"
  )
  out
}

# Convert event-to-exon mappings into a binary event-by-exon membership matrix.
cna_exon_membership_matrix <- function(
    starts,
    widths,
    exons,
    min.exon.coding.overlap = 0) {
  if (length(starts) != length(widths)) {
    stop("starts and widths must have the same length.")
  }
  ends <- starts + widths - 1
  membership <- vapply(seq_len(nrow(exons)), function(j) {
    overlap <- pmax(
      0,
      pmin(ends, exons$cds_genomic_end[[j]]) -
        pmax(starts, exons$cds_genomic_start[[j]]) + 1
    )
    overlap / exons$cds_length_nt[[j]] > min.exon.coding.overlap
  }, logical(length(starts)))
  membership <- matrix(
    membership,
    nrow = length(starts),
    ncol = nrow(exons),
    dimnames = list(NULL, as.character(exons$exon_order))
  )
  membership
}

# Count affected subjects and events for each exon in one recurrence scope.
exon_recurrence_count_table <- function(
    events,
    membership,
    exons,
    cna.types,
    include.combined.analysis = TRUE,
    recurrence.scope = "ALL_OVERLAP") {
  labels <- cna.types
  if (include.combined.analysis && length(cna.types) > 1L) {
    labels <- c(labels, "ALL_CNA")
  }
  rows <- vector("list", length(labels) * nrow(exons))
  index <- 0L
  for (type in labels) {
    type.keep <- if (type == "ALL_CNA") {
      rep(TRUE, nrow(events))
    } else {
      events$cna_type == type
    }
    for (j in seq_len(nrow(exons))) {
      index <- index + 1L
      keep <- type.keep & membership[, j]
      rows[[index]] <- data.frame(
        recurrence_scope = recurrence.scope,
        analysis_type = type,
        exon_order = exons$exon_order[[j]],
        transcript_exon_rank = exons$transcript_exon_rank[[j]],
        exon_id = exons$exon_id[[j]],
        n_coding_residues = exons$n_coding_residues[[j]],
        n_structural_residues = exons$n_structural_residues[[j]],
        mapped_fraction = exons$mapped_fraction[[j]],
        structurally_eligible = exons$structurally_eligible[[j]],
        structural_eligibility_reason =
          exons$structural_eligibility_reason[[j]],
        n_events = length(unique(events$event_id[keep])),
        n_subjects = length(unique(events$subject_id[keep])),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out$test_key <- paste(out$analysis_type, out$exon_order, sep = "|")
  out$row_key <- paste(
    out$recurrence_scope, out$analysis_type, out$exon_order, sep = "|"
  )
  rownames(out) <- NULL
  out
}

# Return exon recurrence counts as fixed-length vectors for null simulations.
exon_recurrence_count_vectors <- function(
    events,
    membership,
    analysis.labels) {
  n.tests <- length(analysis.labels) * ncol(membership)
  event.counts <- integer(n.tests)
  subject.counts <- integer(n.tests)
  index <- 0L
  for (type in analysis.labels) {
    type.keep <- if (type == "ALL_CNA") {
      rep(TRUE, nrow(events))
    } else {
      events$cna_type == type
    }
    for (j in seq_len(ncol(membership))) {
      index <- index + 1L
      keep <- type.keep & membership[, j]
      event.counts[[index]] <- sum(keep)
      subject.counts[[index]] <- length(unique(events$subject_id[keep]))
    }
  }
  list(events = event.counts, subjects = subject.counts)
}

# Compare observed exon recurrence with simulated counts and calculate exon-level
# and search-corrected empirical p-values.
calibrate_cna_exon_recurrence <- function(
    observed,
    null.events,
    null.subjects,
    alpha = 0.05) {
  n.sim <- nrow(null.events)
  if (!identical(dim(null.events), dim(null.subjects))) {
    stop("Event and subject recurrence null matrices must have equal dimensions.")
  }

  observed$n_unique_null_event_counts <- vapply(
    seq_len(ncol(null.events)),
    function(j) length(unique(null.events[, j])),
    integer(1L)
  )
  observed$n_unique_null_subject_counts <- vapply(
    seq_len(ncol(null.subjects)),
    function(j) length(unique(null.subjects[, j])),
    integer(1L)
  )
  observed$null_event_count_sd <- vapply(
    seq_len(ncol(null.events)),
    function(j) if (n.sim > 1L) stats::sd(null.events[, j]) else 0,
    numeric(1L)
  )
  observed$null_subject_count_sd <- vapply(
    seq_len(ncol(null.subjects)),
    function(j) if (n.sim > 1L) stats::sd(null.subjects[, j]) else 0,
    numeric(1L)
  )
  observed$event_recurrence_testable <-
    observed$n_unique_null_event_counts > 1L
  observed$subject_recurrence_testable <-
    observed$n_unique_null_subject_counts > 1L
  observed$event_recurrence_testability_reason <- ifelse(
    observed$event_recurrence_testable,
    "",
    "invariant_null_event_count"
  )
  observed$subject_recurrence_testability_reason <- ifelse(
    observed$subject_recurrence_testable,
    "",
    "invariant_null_subject_count"
  )

  observed$p_events_exon <- vapply(seq_len(nrow(observed)), function(j) {
    cna_empirical_upper_p(observed$n_events[[j]], null.events[, j])
  }, numeric(1L))
  observed$p_subjects_exon <- vapply(seq_len(nrow(observed)), function(j) {
    cna_empirical_upper_p(observed$n_subjects[[j]], null.subjects[, j])
  }, numeric(1L))

  null.p.events <- matrix(
    NA_real_, nrow = n.sim, ncol = ncol(null.events),
    dimnames = dimnames(null.events)
  )
  null.p.subjects <- null.p.events
  for (j in seq_len(ncol(null.events))) {
    null.p.events[, j] <- cna_empirical_upper_p_from_reference(
      null.events[, j], sort(null.events[, j])
    )
    null.p.subjects[, j] <- cna_empirical_upper_p_from_reference(
      null.subjects[, j], sort(null.subjects[, j])
    )
  }

  observed$p_events_typewide <- NA_real_
  observed$p_subjects_typewide <- NA_real_
  for (type in unique(observed$analysis_type)) {
    columns <- which(observed$analysis_type == type)
    min.event.p <- apply(null.p.events[, columns, drop = FALSE], 1L, min)
    min.subject.p <- apply(
      null.p.subjects[, columns, drop = FALSE], 1L, min
    )
    observed$p_events_typewide[columns] <- vapply(columns, function(j) {
      (1 + sum(min.event.p <= observed$p_events_exon[[j]])) / (n.sim + 1)
    }, numeric(1L))
    observed$p_subjects_typewide[columns] <- vapply(columns, function(j) {
      (1 + sum(min.subject.p <= observed$p_subjects_exon[[j]])) /
        (n.sim + 1)
    }, numeric(1L))
  }

  min.event.p.all <- apply(null.p.events, 1L, min)
  min.subject.p.all <- apply(null.p.subjects, 1L, min)
  observed$p_events_omnibus <- vapply(seq_len(nrow(observed)), function(j) {
    (1 + sum(min.event.p.all <= observed$p_events_exon[[j]])) / (n.sim + 1)
  }, numeric(1L))
  observed$p_subjects_omnibus <- vapply(seq_len(nrow(observed)), function(j) {
    (1 + sum(min.subject.p.all <= observed$p_subjects_exon[[j]])) /
      (n.sim + 1)
  }, numeric(1L))
  observed$p_events_bonferroni <- pmin(
    1, observed$p_events_exon * nrow(observed)
  )
  observed$p_subjects_bonferroni <- pmin(
    1, observed$p_subjects_exon * nrow(observed)
  )

  event.p.columns <- c(
    "p_events_exon", "p_events_typewide", "p_events_omnibus",
    "p_events_bonferroni"
  )
  subject.p.columns <- c(
    "p_subjects_exon", "p_subjects_typewide", "p_subjects_omnibus",
    "p_subjects_bonferroni"
  )
  for (column in event.p.columns) {
    observed[[column]][!observed$event_recurrence_testable] <- NA_real_
  }
  for (column in subject.p.columns) {
    observed[[column]][!observed$subject_recurrence_testable] <- NA_real_
  }
  observed$significant_subject_recurrence <-
    observed$subject_recurrence_testable &
      !is.na(observed$p_subjects_omnibus) &
      observed$p_subjects_omnibus <= alpha
  observed$significant_event_recurrence <-
    observed$event_recurrence_testable &
      !is.na(observed$p_events_omnibus) &
      observed$p_events_omnibus <= alpha
  observed$recurrence_class <- ifelse(
    !observed$subject_recurrence_testable &
      !observed$event_recurrence_testable,
    "Recurrence not testable: invariant null",
    ifelse(
      observed$significant_subject_recurrence &
        observed$significant_event_recurrence,
      "Subject and event recurrence",
      ifelse(
        observed$significant_subject_recurrence,
        "Subject recurrence",
        ifelse(
          observed$significant_event_recurrence,
          "Event recurrence",
          ifelse(
            !observed$subject_recurrence_testable,
            "Subject recurrence not testable",
            ifelse(
              !observed$event_recurrence_testable,
              "Event recurrence not testable",
              "No significant exon recurrence"
            )
          )
        )
      )
    )
  )

  list(
    results = observed,
    null_event_counts = null.events,
    null_subject_counts = null.subjects,
    null_event_p = null.p.events,
    null_subject_p = null.p.subjects,
    null_min_event_p = min.event.p.all,
    null_min_subject_p = min.subject.p.all
  )
}

# Simulate intact genomic intervals conditional on coding-sequence overlap to
# calibrate coding-coverage and ALL_OVERLAP exon-recurrence statistics.
run_cna_coverage_null <- function(
    events,
    exons,
    cna.types,
    thresholds = c(0.50, 0.75, 0.90, 1.00),
    n.sim = 10000L,
    chromosome.length = NULL,
    min.exon.coding.overlap = 0,
    include.combined.analysis = TRUE,
    alpha = 0.05,
    random.seed = 20260907L,
    progress.every = 100L) {
  if (n.sim < 1L) stop("n.sim.coverage must be at least 1.")
  thresholds <- sort(unique(as.numeric(thresholds)))
  if (anyNA(thresholds) || any(thresholds <= 0 | thresholds > 1)) {
    stop("coverage.thresholds must be values in (0, 1].")
  }

  observed <- coverage_count_table(
    events, events$coding_fraction, cna.types, thresholds,
    include.combined.analysis
  )
  observed.membership <- cna_exon_membership_matrix(
    events$genomic_start,
    events$width_bp,
    exons,
    min.exon.coding.overlap
  )
  observed.recurrence <- exon_recurrence_count_table(
    events,
    observed.membership,
    exons,
    cna.types,
    include.combined.analysis,
    recurrence.scope = "ALL_OVERLAP"
  )
  recurrence.labels <- unique(observed.recurrence$analysis_type)
  null.events <- matrix(
    0L, nrow = n.sim, ncol = nrow(observed),
    dimnames = list(NULL, observed$row_key)
  )
  null.subjects <- null.events
  null.recurrence.events <- matrix(
    0L, nrow = n.sim, ncol = nrow(observed.recurrence),
    dimnames = list(NULL, observed.recurrence$row_key)
  )
  null.recurrence.subjects <- null.recurrence.events
  legal.starts <- lapply(
    events$width_bp,
    legal_cna_starts_conditioned_on_cds,
    exons = exons,
    chromosome.length = chromosome.length
  )

  set.seed(random.seed)
  for (b in seq_len(n.sim)) {
    simulated.starts <- vapply(seq_len(nrow(events)), function(i) {
      sample_integer_from_intervals(legal.starts[[i]])
    }, numeric(1L))
    simulated.fractions <- vapply(seq_len(nrow(events)), function(i) {
      coding_fraction_for_interval(
        simulated.starts[[i]],
        simulated.starts[[i]] + events$width_bp[[i]] - 1,
        exons
      )
    }, numeric(1L))
    one <- coverage_count_table(
      events, simulated.fractions, cna.types, thresholds,
      include.combined.analysis
    )
    one <- one[match(observed$row_key, one$row_key), , drop = FALSE]
    null.events[b, ] <- one$n_events
    null.subjects[b, ] <- one$n_subjects

    simulated.membership <- cna_exon_membership_matrix(
      simulated.starts,
      events$width_bp,
      exons,
      min.exon.coding.overlap
    )
    one.recurrence <- exon_recurrence_count_vectors(
      events,
      simulated.membership,
      recurrence.labels
    )
    null.recurrence.events[b, ] <- one.recurrence$events
    null.recurrence.subjects[b, ] <- one.recurrence$subjects

    if (progress.every > 0L &&
        (b %% progress.every == 0L || b == n.sim)) {
      message("Completed ", b, " of ", n.sim, " coverage simulations.")
    }
  }

  observed$p_events_conditional <- vapply(seq_len(nrow(observed)), function(i) {
    cna_empirical_upper_p(observed$n_events[[i]], null.events[, i])
  }, numeric(1L))
  observed$p_subjects_conditional <- vapply(
    seq_len(nrow(observed)),
    function(i) {
      cna_empirical_upper_p(observed$n_subjects[[i]], null.subjects[, i])
    },
    numeric(1L)
  )
  observed$n_unique_null_event_counts <- vapply(
    seq_len(ncol(null.events)),
    function(j) length(unique(null.events[, j])),
    integer(1L)
  )
  observed$n_unique_null_subject_counts <- vapply(
    seq_len(ncol(null.subjects)),
    function(j) length(unique(null.subjects[, j])),
    integer(1L)
  )
  observed$event_coverage_testable <-
    observed$n_unique_null_event_counts > 1L
  observed$subject_coverage_testable <-
    observed$n_unique_null_subject_counts > 1L
  observed$event_coverage_testability_reason <- ifelse(
    observed$event_coverage_testable,
    "",
    "invariant_null_event_count"
  )
  observed$subject_coverage_testability_reason <- ifelse(
    observed$subject_coverage_testable,
    "",
    "invariant_null_subject_count"
  )
  n.tests <- nrow(observed)
  observed$p_events_bonferroni <- pmin(
    1, observed$p_events_conditional * n.tests
  )
  observed$p_subjects_bonferroni <- pmin(
    1, observed$p_subjects_conditional * n.tests
  )
  observed$p_events_conditional[!observed$event_coverage_testable] <-
    NA_real_
  observed$p_events_bonferroni[!observed$event_coverage_testable] <-
    NA_real_
  observed$p_subjects_conditional[!observed$subject_coverage_testable] <-
    NA_real_
  observed$p_subjects_bonferroni[!observed$subject_coverage_testable] <-
    NA_real_
  observed$significant_subject_coverage <-
    observed$subject_coverage_testable &
      !is.na(observed$p_subjects_bonferroni) &
      observed$p_subjects_bonferroni <= alpha
  observed$significant_event_coverage <-
    observed$event_coverage_testable &
      !is.na(observed$p_events_bonferroni) &
      observed$p_events_bonferroni <= alpha
  observed$interpretation_scope <-
    "Within-gene coding coverage conditional on a CNA overlapping the CDS"
  observed$coverage_threshold_percent <- 100 * observed$coverage_threshold
  observed <- observed[, c(
    "analysis_type", "coverage_threshold", "coverage_threshold_percent",
    "n_subjects", "n_events", "p_subjects_conditional",
    "p_events_conditional", "p_subjects_bonferroni",
    "p_events_bonferroni", "significant_subject_coverage",
    "significant_event_coverage", "subject_coverage_testable",
    "event_coverage_testable", "n_unique_null_subject_counts",
    "n_unique_null_event_counts", "subject_coverage_testability_reason",
    "event_coverage_testability_reason", "interpretation_scope", "row_key"
  )]

  recurrence <- calibrate_cna_exon_recurrence(
    observed.recurrence,
    null.recurrence.events,
    null.recurrence.subjects,
    alpha
  )

  list(
    results = observed,
    null_subject_counts = null.subjects,
    null_event_counts = null.events,
    exon_recurrence = recurrence
  )
}

# Build a binary event-by-exon matrix from an existing mapping table.
membership_matrix_from_exon_mapping <- function(events, mapping, exons) {
  membership <- matrix(
    FALSE,
    nrow = nrow(events),
    ncol = nrow(exons),
    dimnames = list(events$event_id, as.character(exons$exon_order))
  )
  if (!nrow(events) || !nrow(mapping)) return(membership)

  z <- mapping[mapping$event_id %in% events$event_id, , drop = FALSE]
  if ("passes_exon_overlap" %in% names(z)) {
    z <- z[z$passes_exon_overlap, , drop = FALSE]
  }
  if (!nrow(z)) return(membership)

  row.index <- match(z$event_id, events$event_id)
  column.index <- match(z$exon_order, exons$exon_order)
  keep <- !is.na(row.index) & !is.na(column.index)
  membership[cbind(row.index[keep], column.index[keep])] <- TRUE
  membership
}

# Combine recurrence results from ALL_OVERLAP and PARTIAL_LOCALIZABLE analyses
# into one consistently formatted output table.
combine_exon_recurrence_scopes <- function(scope.results, alpha = 0.05) {
  scope.results <- scope.results[
    !vapply(scope.results, is.null, logical(1L))
  ]
  if (!length(scope.results)) return(data.frame())

  scope.tables <- lapply(scope.results, function(x) x$results)
  # V4 adds permutation-opportunity summaries to PARTIAL_LOCALIZABLE rows.
  # Align the two schemas before binding and leave those fields unavailable
  # for ALL_OVERLAP, whose null operates on genomic intervals instead.
  all.columns <- unique(unlist(lapply(scope.tables, names), use.names = FALSE))
  scope.tables <- lapply(scope.tables, function(x) {
    missing.columns <- setdiff(all.columns, names(x))
    for (column in missing.columns) x[[column]] <- NA
    x[, all.columns, drop = FALSE]
  })
  out <- do.call(rbind, scope.tables)
  rownames(out) <- NULL
  n.scopes <- length(unique(out$recurrence_scope))
  out$recurrence_scope_description <- ifelse(
    out$recurrence_scope == "ALL_OVERLAP",
    paste(
      "All analyzed CNAs overlapping the coding sequence; genomic width",
      "is preserved in the conditional interval-placement null"
    ),
    paste(
      "Only partial structurally localizable CNAs; subject, type, exon span",
      "and approximate coding-residue span are preserved in the null"
    )
  )

  # p_*_omnibus already corrects across all exons and CNA types within one
  # recurrence scope. The additional Bonferroni factor protects an analysis
  # that examines both prespecified scopes and reports the more favorable one.
  out$p_events_scope_omnibus <- out$p_events_omnibus
  out$p_subjects_scope_omnibus <- out$p_subjects_omnibus
  out$p_events_across_scopes <- pmin(
    1, out$p_events_scope_omnibus * n.scopes
  )
  out$p_subjects_across_scopes <- pmin(
    1, out$p_subjects_scope_omnibus * n.scopes
  )
  out$significant_subject_recurrence_within_scope <-
    out$significant_subject_recurrence
  out$significant_event_recurrence_within_scope <-
    out$significant_event_recurrence
  out$significant_subject_recurrence <-
    out$subject_recurrence_testable &
      !is.na(out$p_subjects_across_scopes) &
      out$p_subjects_across_scopes <= alpha
  out$significant_event_recurrence <-
    out$event_recurrence_testable &
      !is.na(out$p_events_across_scopes) &
      out$p_events_across_scopes <= alpha
  out$recurrence_class <- ifelse(
    !out$subject_recurrence_testable &
      !out$event_recurrence_testable,
    "Recurrence not testable: invariant null",
    ifelse(
      out$significant_subject_recurrence &
        out$significant_event_recurrence,
      "Subject and event recurrence",
      ifelse(
        out$significant_subject_recurrence,
        "Subject recurrence",
        ifelse(
          out$significant_event_recurrence,
          "Event recurrence",
          ifelse(
            !out$subject_recurrence_testable,
            "Subject recurrence not testable",
            ifelse(
              !out$event_recurrence_testable,
              "Event recurrence not testable",
              "No significant exon recurrence"
            )
          )
        )
      )
    )
  )

  leading <- c(
    "recurrence_scope", "recurrence_scope_description",
    "analysis_type", "exon_order",
    "transcript_exon_rank", "exon_id", "n_coding_residues",
    "n_structural_residues", "mapped_fraction",
    "structurally_eligible", "structural_eligibility_reason",
    "n_events", "n_subjects", "test_key", "row_key"
  )
  out[, c(leading, setdiff(names(out), leading)), drop = FALSE]
}

# Calibrate exon recurrence among partial, structurally localizable CNAs by
# relocating each intact exon span among its eligible transcript windows.
run_partial_localizable_recurrence_null <- function(
    events,
    observed.mapping,
    exons,
    exon.residue.map,
    cna.types,
    include.combined.analysis = TRUE,
    n.sim = 10000L,
    max.residue.fraction = 0.999999,
    residue.span.temperature = 0.50,
    uniform.window.mixture = 0.25,
    min.effective.windows = 1.25,
    min.alternative.probability = 0.05,
    avoid.within.subject.overlap = TRUE,
    alpha = 0.05,
    random.seed = 20260907L,
    progress.every = 100L) {
  if (n.sim < 1L) stop("n.sim.coverage must be at least 1.")

  if (!"span_exons" %in% names(events)) {
    events$span_exons <- events$n_affected_exons
  }
  if (nrow(events)) {
    events$span_exons <- normalize_cna_exon_spans(
      events$span_exons,
      nrow(exons),
      context = "PARTIAL_LOCALIZABLE CNA events"
    )
  }
  observed.membership <- membership_matrix_from_exon_mapping(
    events, observed.mapping, exons
  )
  observed <- exon_recurrence_count_table(
    events,
    observed.membership,
    exons,
    cna.types,
    include.combined.analysis,
    recurrence.scope = "PARTIAL_LOCALIZABLE"
  )
  labels <- unique(observed$analysis_type)
  null.events <- matrix(
    0L, nrow = n.sim, ncol = nrow(observed),
    dimnames = list(NULL, observed$row_key)
  )
  null.subjects <- null.events
  # Initialize the return object with a zero-row event table. Passing the full
  # event table here would attempt to retrieve windows before legal.windows is
  # constructed below.
  permutation.diagnostics <- build_cna_permutation_diagnostics(
    events = events[0, , drop = FALSE],
    legal.windows = list(),
    residue.span.temperature = residue.span.temperature,
    uniform.window.mixture = uniform.window.mixture,
    min.effective.windows = min.effective.windows,
    min.alternative.probability = min.alternative.probability
  )

  if (nrow(events)) {
    legal.windows <- build_legal_cna_windows(
      events,
      exons,
      exon.residue.map,
      max.residue.fraction
    )
    permutation.diagnostics <- build_cna_permutation_diagnostics(
      events,
      legal.windows,
      residue.span.temperature,
      uniform.window.mixture,
      min.effective.windows,
      min.alternative.probability
    )
    set.seed(random.seed)
    for (b in seq_len(n.sim)) {
      simulated.mapping <- simulate_cna_event_mapping(
        events,
        legal.windows,
        exons,
        residue.span.temperature,
        uniform.window.mixture,
        avoid.within.subject.overlap
      )
      simulated.membership <- membership_matrix_from_exon_mapping(
        events, simulated.mapping, exons
      )
      one <- exon_recurrence_count_vectors(
        events, simulated.membership, labels
      )
      null.events[b, ] <- one$events
      null.subjects[b, ] <- one$subjects

      if (progress.every > 0L &&
          (b %% progress.every == 0L || b == n.sim)) {
        message(
          "Completed ", b, " of ", n.sim,
          " partial-localizable recurrence simulations."
        )
      }
    }
  }

  calibrated <- calibrate_cna_exon_recurrence(
    observed, null.events, null.subjects, alpha
  )
  calibrated$results$n_scope_events_total <- 0L
  calibrated$results$n_events_with_alternative_windows <- 0L
  calibrated$results$n_events_permutation_testable <- 0L
  calibrated$results$fraction_events_permutation_testable <- NA_real_
  calibrated$results$median_effective_number_of_windows <- NA_real_
  calibrated$results$minimum_effective_number_of_windows <- NA_real_
  calibrated$results$maximum_observed_window_probability <- NA_real_

  for (type in unique(calibrated$results$analysis_type)) {
    keep.result <- calibrated$results$analysis_type == type
    z <- if (type == "ALL_CNA") {
      permutation.diagnostics
    } else {
      permutation.diagnostics[
        permutation.diagnostics$cna_type == type, , drop = FALSE
      ]
    }
    calibrated$results$n_scope_events_total[keep.result] <- nrow(z)
    if (!nrow(z)) next
    calibrated$results$n_events_with_alternative_windows[keep.result] <-
      sum(z$n_alternative_windows > 0L)
    calibrated$results$n_events_permutation_testable[keep.result] <-
      sum(z$permutation_testable)
    calibrated$results$fraction_events_permutation_testable[keep.result] <-
      mean(z$permutation_testable)
    calibrated$results$median_effective_number_of_windows[keep.result] <-
      stats::median(z$effective_number_of_windows)
    calibrated$results$minimum_effective_number_of_windows[keep.result] <-
      min(z$effective_number_of_windows)
    observed.probabilities <- z$probability_observed_window[
      is.finite(z$probability_observed_window)
    ]
    if (length(observed.probabilities)) {
      calibrated$results$maximum_observed_window_probability[keep.result] <-
        max(observed.probabilities)
    }
  }
  calibrated$permutation_diagnostics <- permutation.diagnostics
  calibrated
}


# ==============================================================================
# 5. EXON-TO-EXON STRUCTURAL DISTANCES
# ==============================================================================

# Calculate eligible C-alpha distances between residues encoded by two exons.
cross_exon_distances <- function(
    exon.a,
    exon.b,
    exon.residue.map,
    coordinates,
    min.sequence.separation = 10L) {
  a <- exon.residue.map[
    exon.residue.map$exon_order == exon.a &
      exon.residue.map$structurally_eligible,
    , drop = FALSE
  ]
  b <- exon.residue.map[
    exon.residue.map$exon_order == exon.b &
      exon.residue.map$structurally_eligible,
    , drop = FALSE
  ]
  xyz.a <- as.matrix(
    coordinates[a$coordinate_row, c("x", "y", "z"), drop = FALSE]
  )
  xyz.b <- as.matrix(
    coordinates[b$coordinate_row, c("x", "y", "z"), drop = FALSE]
  )

  dx <- outer(xyz.a[, 1L], xyz.b[, 1L], "-")
  dy <- outer(xyz.a[, 2L], xyz.b[, 2L], "-")
  dz <- outer(xyz.a[, 3L], xyz.b[, 3L], "-")
  distances <- sqrt(dx^2 + dy^2 + dz^2)
  separation <- abs(outer(a$residue, b$residue, "-"))
  keep <- separation >= min.sequence.separation
  fallback <- FALSE
  if (!any(keep)) {
    keep[,] <- TRUE
    fallback <- TRUE
  }

  list(
    distances = as.numeric(distances[keep]),
    used_sequence_fallback = fallback
  )
}

# Build a robust exon-to-exon distance matrix using the mean of the closest
# eligible cross-exon residue pairs.
build_exon_distance_matrix <- function(
    exons,
    exon.residue.map,
    coordinates,
    n.closest = 5L,
    min.sequence.separation = 10L) {
  eligible <- exons[exons$structurally_eligible, , drop = FALSE]
  if (n.closest < 1L) stop("n.closest must be at least 1.")

  orders <- eligible$exon_order
  robust.matrix <- matrix(
    0, nrow = length(orders), ncol = length(orders),
    dimnames = list(as.character(orders), as.character(orders))
  )
  if (nrow(eligible) < 2L) {
    pair.table <- data.frame(
      exon_order_1 = integer(), transcript_exon_rank_1 = integer(),
      exon_id_1 = character(), exon_order_2 = integer(),
      transcript_exon_rank_2 = integer(), exon_id_2 = character(),
      adjacent_in_transcript = logical(), n_cross_residue_pairs = integer(),
      n_smallest_distances_used = integer(), minimum_distance = numeric(),
      robust_distance = numeric(), median_cross_distance = numeric(),
      sequence_filter_fallback = logical(), stringsAsFactors = FALSE
    )
    return(list(
      distance_matrix = robust.matrix,
      pair_table = pair.table,
      eligible_exons = eligible
    ))
  }
  pair.rows <- vector("list", choose(length(orders), 2L))
  row.index <- 0L

  for (i in seq_len(length(orders) - 1L)) {
    for (j in seq.int(i + 1L, length(orders))) {
      row.index <- row.index + 1L
      cross <- cross_exon_distances(
        orders[[i]], orders[[j]], exon.residue.map, coordinates,
        min.sequence.separation
      )
      sorted <- sort(cross$distances)
      n.used <- min(n.closest, length(sorted))
      robust <- mean(sorted[seq_len(n.used)])
      robust.matrix[i, j] <- robust
      robust.matrix[j, i] <- robust

      pair.rows[[row.index]] <- data.frame(
        exon_order_1 = orders[[i]],
        transcript_exon_rank_1 = eligible$transcript_exon_rank[[i]],
        exon_id_1 = eligible$exon_id[[i]],
        exon_order_2 = orders[[j]],
        transcript_exon_rank_2 = eligible$transcript_exon_rank[[j]],
        exon_id_2 = eligible$exon_id[[j]],
        adjacent_in_transcript = abs(orders[[i]] - orders[[j]]) == 1L,
        n_cross_residue_pairs = length(sorted),
        n_smallest_distances_used = n.used,
        minimum_distance = min(sorted),
        robust_distance = robust,
        median_cross_distance = stats::median(sorted),
        sequence_filter_fallback = cross$used_sequence_fallback,
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    distance_matrix = robust.matrix,
    pair_table = do.call(rbind, pair.rows),
    eligible_exons = eligible
  )
}


# ==============================================================================
# 6. ADAPTIVE 1D AND 3D COMPLETE-LINKAGE EXON TREES
# ==============================================================================

# Recover the exon sets represented by every internal branch of an hclust tree.
extract_cna_hclust_sets <- function(hc) {
  n <- length(hc$order)
  if (n <= 1L) return(list())
  sets <- vector("list", n - 1L)
  for (i in seq_len(n - 1L)) {
    children <- hc$merge[i, ]
    members <- integer()
    for (child in children) {
      members <- c(members, if (child < 0L) -child else sets[[child]])
    }
    sets[[i]] <- sort(unique(members))
  }
  sets
}

# Construct a hierarchical-clustering tree and return its single-exon and merged sets.
build_cna_tree_sets <- function(items, distance.matrix, linkage.method = "complete") {
  validate_linkage_method(linkage.method)
  leaf.sets <- lapply(items, function(x) x)
  if (length(items) == 1L) return(leaf.sets)
  hc <- stats::hclust(stats::as.dist(distance.matrix), method = linkage.method)
  internal <- lapply(extract_cna_hclust_sets(hc), function(i) items[i])
  c(leaf.sets, internal)
}

# Calculate the 1D exon-order span and 3D complete-linkage diameter of a set.
cna_cluster_diameters <- function(exon.set, exon.distance.matrix) {
  exon.set <- sort(unique(as.integer(exon.set)))
  diameter.1d <- if (length(exon.set) == 1L) {
    0
  } else {
    max(exon.set) - min(exon.set)
  }
  diameter.3d <- if (length(exon.set) == 1L) {
    0
  } else {
    max(exon.distance.matrix[
      as.character(exon.set), as.character(exon.set), drop = FALSE
    ])
  }
  c(diameter_1d = diameter.1d, diameter_3d = diameter.3d)
}

# Generate and deduplicate candidate CNA clusters from the 1D and 3D exon trees,
# retaining only clusters that meet the subject and event thresholds.
build_cna_candidate_clusters <- function(
    mapping,
    exons,
    exon.distance.matrix,
    min.subjects = 2L,
    min.events = 2L,
    linkage.method = "complete") {
  validate_linkage_method(linkage.method)
  altered <- sort(unique(mapping$exon_order))
  if (!length(altered)) return(data.frame())
  if (!all(as.character(altered) %in% rownames(exon.distance.matrix))) {
    stop("An altered exon is absent from the structural distance matrix.")
  }

  distance.1d <- abs(outer(altered, altered, "-"))
  distance.3d <- exon.distance.matrix[
    as.character(altered), as.character(altered), drop = FALSE
  ]
  sets.1d <- build_cna_tree_sets(altered, distance.1d, linkage.method = linkage.method)
  sets.3d <- build_cna_tree_sets(altered, distance.3d, linkage.method = linkage.method)
  source.map <- new.env(parent = emptyenv())
  exon.map <- new.env(parent = emptyenv())

  register <- function(sets, source) {
    for (exon.set in sets) {
      key <- paste(sort(unique(exon.set)), collapse = ";")
      previous <- if (exists(key, source.map, inherits = FALSE)) {
        get(key, source.map, inherits = FALSE)
      } else {
        character()
      }
      assign(key, sort(unique(c(previous, source))), source.map)
      assign(key, sort(unique(exon.set)), exon.map)
    }
  }
  register(sets.1d, "1D")
  register(sets.3d, "3D")

  keys <- ls(source.map, all.names = TRUE)
  rows <- lapply(keys, function(key) {
    exon.set <- get(key, exon.map, inherits = FALSE)
    z <- mapping[mapping$exon_order %in% exon.set, , drop = FALSE]
    diameter <- cna_cluster_diameters(exon.set, exon.distance.matrix)
    data.frame(
      candidate_key = key,
      tree_sources = paste(get(key, source.map, inherits = FALSE), collapse = ";"),
      exon_orders = collapse_grin3d_values(exon.set),
      transcript_exon_ranks = collapse_grin3d_values(
        exons$transcript_exon_rank[match(exon.set, exons$exon_order)]
      ),
      exon_ids = collapse_grin3d_values(
        exons$exon_id[match(exon.set, exons$exon_order)]
      ),
      exon_order_min = min(exon.set),
      exon_order_max = max(exon.set),
      n_exons = length(exon.set),
      n_events = length(unique(z$event_id)),
      n_subjects = length(unique(z$subject_id)),
      diameter_1d = unname(diameter[["diameter_1d"]]),
      diameter_3d = unname(diameter[["diameter_3d"]]),
      stringsAsFactors = FALSE
    )
  })
  candidates <- do.call(rbind, rows)
  candidates <- candidates[
    candidates$n_subjects >= min.subjects &
      candidates$n_events >= min.events,
    , drop = FALSE
  ]
  rownames(candidates) <- NULL
  candidates
}

# Add lesion-type-specific subject and event counts to each candidate cluster.
add_cluster_type_counts <- function(clusters, mapping, cna.types) {
  for (type in c(cna.types, "ALL_CNA")) {
    safe <- safe_grin3d_name(type)
    event.column <- paste0("n_events_", safe)
    subject.column <- paste0("n_subjects_", safe)
    clusters[[event.column]] <- 0L
    clusters[[subject.column]] <- 0L
    for (i in seq_len(nrow(clusters))) {
      exons <- as.integer(strsplit(
        clusters$exon_orders[[i]], ";", fixed = TRUE
      )[[1L]])
      z <- mapping[mapping$exon_order %in% exons, , drop = FALSE]
      if (type != "ALL_CNA") z <- z[z$cna_type == type, , drop = FALSE]
      clusters[[event.column]][[i]] <- length(unique(z$event_id))
      clusters[[subject.column]][[i]] <- length(unique(z$subject_id))
    }
  }
  clusters
}

# Quantify the subjects and events added beyond each cluster's best-supported
# proper subcluster.
add_cluster_support_metrics <- function(clusters, mapping) {
  if (!nrow(clusters)) return(clusters)
  parse_set <- function(x) {
    as.integer(strsplit(x, ";", fixed = TRUE)[[1L]])
  }
  cluster.sets <- lapply(clusters$exon_orders, parse_set)

  clusters$events_per_exon <- ""
  clusters$subjects_per_exon <- ""
  clusters$min_events_per_exon <- 0L
  clusters$min_subjects_per_exon <- 0L
  clusters$n_events_spanning_multiple_cluster_exons <- 0L
  clusters$n_subjects_spanning_multiple_cluster_exons <- 0L
  clusters$n_subjects_with_multiple_events_in_cluster <- 0L
  clusters$best_supported_subcluster_key <- ""
  clusters$n_events_added_vs_best_subcluster <- NA_integer_
  clusters$n_subjects_added_vs_best_subcluster <- NA_integer_
  clusters$branch_adds_independent_subjects <- NA

  for (i in seq_len(nrow(clusters))) {
    exon.set <- cluster.sets[[i]]
    z <- mapping[mapping$exon_order %in% exon.set, , drop = FALSE]

    event.counts <- vapply(exon.set, function(exon.order) {
      length(unique(z$event_id[z$exon_order == exon.order]))
    }, integer(1L))
    subject.counts <- vapply(exon.set, function(exon.order) {
      length(unique(z$subject_id[z$exon_order == exon.order]))
    }, integer(1L))
    clusters$events_per_exon[[i]] <- paste(
      paste(exon.set, event.counts, sep = ":"), collapse = ";"
    )
    clusters$subjects_per_exon[[i]] <- paste(
      paste(exon.set, subject.counts, sep = ":"), collapse = ";"
    )
    clusters$min_events_per_exon[[i]] <- min(event.counts)
    clusters$min_subjects_per_exon[[i]] <- min(subject.counts)

    event.exon.count <- tapply(
      z$exon_order,
      z$event_id,
      function(v) length(unique(v))
    )
    multi.event.ids <- names(event.exon.count)[event.exon.count >= 2L]
    clusters$n_events_spanning_multiple_cluster_exons[[i]] <-
      length(multi.event.ids)
    clusters$n_subjects_spanning_multiple_cluster_exons[[i]] <-
      length(unique(z$subject_id[z$event_id %in% multi.event.ids]))
    subject.event.count <- tapply(
      z$event_id,
      z$subject_id,
      function(v) length(unique(v))
    )
    clusters$n_subjects_with_multiple_events_in_cluster[[i]] <-
      sum(subject.event.count >= 2L)

    if (length(exon.set) > 1L) {
      proper <- which(vapply(cluster.sets, function(candidate) {
        length(candidate) < length(exon.set) && all(candidate %in% exon.set)
      }, logical(1L)))
      if (length(proper)) {
        ranked <- proper[order(
          -clusters$n_subjects[proper],
          -clusters$n_events[proper],
          -clusters$n_exons[proper],
          clusters$candidate_key[proper]
        )]
        best <- ranked[[1L]]
        clusters$best_supported_subcluster_key[[i]] <-
          clusters$candidate_key[[best]]
        clusters$n_events_added_vs_best_subcluster[[i]] <-
          clusters$n_events[[i]] - clusters$n_events[[best]]
        clusters$n_subjects_added_vs_best_subcluster[[i]] <-
          clusters$n_subjects[[i]] - clusters$n_subjects[[best]]
        clusters$branch_adds_independent_subjects[[i]] <-
          clusters$n_subjects_added_vs_best_subcluster[[i]] > 0L
      }
    }
  }
  clusters
}

# Link each structural cluster to the strongest exon-recurrence evidence among
# its member exons.
annotate_clusters_with_exon_recurrence <- function(
    clusters,
    recurrence.results,
    alpha = 0.05) {
  if (!nrow(clusters)) return(clusters)
  scopes <- unique(recurrence.results$recurrence_scope)
  clusters$recurrence_scope_used <- if (length(scopes) == 1L) {
    scopes
  } else {
    paste(scopes, collapse = ";")
  }
  clusters$best_recurrent_exon_order <- NA_integer_
  clusters$best_recurrent_transcript_exon_rank <- NA_integer_
  clusters$best_exon_p_subjects_raw <- NA_real_
  clusters$best_exon_p_subjects_scope_omnibus <- NA_real_
  clusters$best_exon_p_subjects_across_scopes <- NA_real_
  clusters$best_exon_p_subjects_omnibus <- NA_real_
  clusters$n_testable_recurrence_exons <- 0L
  clusters$n_significant_recurrent_exons <- 0L
  clusters$integrated_evidence_class <- ""

  for (i in seq_len(nrow(clusters))) {
    exon.set <- as.integer(strsplit(
      clusters$exon_orders[[i]], ";", fixed = TRUE
    )[[1L]])
    z <- recurrence.results[
      recurrence.results$analysis_type == clusters$analysis_type[[i]] &
        recurrence.results$exon_order %in% exon.set,
      , drop = FALSE
    ]
    if (!nrow(z)) {
      clusters$integrated_evidence_class[[i]] <-
        "No matching exon-recurrence result"
      next
    }
    best <- order(
      is.na(z$p_subjects_across_scopes),
      z$p_subjects_across_scopes,
      is.na(z$p_subjects_scope_omnibus),
      z$p_subjects_scope_omnibus,
      z$p_subjects_exon
    )[[1L]]
    clusters$best_recurrent_exon_order[[i]] <- z$exon_order[[best]]
    clusters$best_recurrent_transcript_exon_rank[[i]] <-
      z$transcript_exon_rank[[best]]
    clusters$best_exon_p_subjects_raw[[i]] <- z$p_subjects_exon[[best]]
    clusters$best_exon_p_subjects_scope_omnibus[[i]] <-
      z$p_subjects_scope_omnibus[[best]]
    clusters$best_exon_p_subjects_across_scopes[[i]] <-
      z$p_subjects_across_scopes[[best]]
    # Backward-compatible alias. In v4 this is the fully adjusted p-value
    # across exons, CNA types and the two recurrence scopes.
    clusters$best_exon_p_subjects_omnibus[[i]] <-
      z$p_subjects_across_scopes[[best]]
    clusters$n_testable_recurrence_exons[[i]] <-
      sum(z$subject_recurrence_testable)
    clusters$n_significant_recurrent_exons[[i]] <-
      sum(z$significant_subject_recurrence)

    recurrence.testable <- clusters$n_testable_recurrence_exons[[i]] > 0L
    recurrent <- clusters$n_significant_recurrent_exons[[i]] > 0L
    compactness.testable <- isTRUE(clusters$compactness_testable[[i]])
    compact <- isTRUE(clusters$p_hotspot_omnibus[[i]] <= alpha)
    clusters$integrated_evidence_class[[i]] <- if (
      !recurrence.testable && !compactness.testable
    ) {
      "Exon recurrence and multi-exon compactness not testable"
    } else if (!recurrence.testable && compact) {
      "Significant compact cluster; exon recurrence not testable"
    } else if (!recurrence.testable) {
      "Exon recurrence not testable; no significant compactness"
    } else if (recurrent && !compactness.testable) {
      "Recurrent exon(s); multi-exon compactness not testable"
    } else if (!compactness.testable) {
      "No significant recurrence; multi-exon compactness not testable"
    } else if (recurrent && compact) {
      "Recurrent exon(s) in a significant 1D/3D compact cluster"
    } else if (recurrent) {
      "Recurrent exon(s) without significant multi-exon compactness"
    } else if (compact) {
      "Significant compact cluster without exon-level recurrence"
    } else {
      "No significant recurrence or compactness"
    }
  }
  clusters
}


# ==============================================================================
# 7. SUBJECT-, TYPE- AND SPAN-PRESERVING STRUCTURAL NULL
# ==============================================================================

# Count unique protein residues represented by a consecutive exon window.
count_window_residues <- function(exon.orders, exon.residue.map) {
  length(unique(exon.residue.map$residue[
    exon.residue.map$exon_order %in% exon.orders
  ]))
}

# Validate event exon spans and cap them at the number of available exons.
normalize_cna_exon_spans <- function(
    spans,
    n.exons,
    context = "CNA events") {
  normalized <- suppressWarnings(as.integer(as.character(spans)))
  invalid <- is.na(normalized) | normalized < 1L | normalized >= n.exons
  if (any(invalid)) {
    bad <- unique(as.character(spans[invalid]))
    stop(
      context, " contains invalid partial-CNA exon span(s): ",
      paste(bad, collapse = "; "), ". Valid spans for a transcript with ",
      n.exons, " coding exons are 1 through ", n.exons - 1L, "."
    )
  }
  normalized
}

# Retrieve the eligible relocation windows for one observed CNA exon span.
get_cna_legal_windows <- function(
    legal.windows,
    span.exons,
    event.id = NA_character_) {
  span <- suppressWarnings(as.integer(as.character(span.exons[[1L]])))
  event.label <- if (!is.na(event.id) && nzchar(event.id)) {
    paste0(" for event '", event.id, "'")
  } else {
    ""
  }
  if (is.na(span) || span < 1L) {
    stop("Invalid affected-exon span", event.label, ": ", span.exons[[1L]], ".")
  }

  key <- as.character(span)
  if (is.null(names(legal.windows)) || !key %in% names(legal.windows)) {
    stop(
      "No legal-window table was constructed", event.label,
      " for an event spanning ", span, " exon(s). Available spans: ",
      paste(names(legal.windows), collapse = "; "), "."
    )
  }
  candidates <- legal.windows[[key]]
  if (is.null(candidates)) {
    stop("The legal-window table is NULL", event.label, " for span ", span, ".")
  }
  candidates <- as.data.frame(candidates, stringsAsFactors = FALSE)
  if (base::nrow(candidates) == 0L) {
    stop("The legal-window table has no rows", event.label, " for span ", span, ".")
  }
  required <- c(
    "exon_start", "exon_end", "span_exons", "n_coding_residues",
    "structurally_eligible", "residue_fraction"
  )
  missing <- setdiff(required, names(candidates))
  if (length(missing)) {
    stop(
      "The legal-window table", event.label, " is missing column(s): ",
      paste(missing, collapse = "; "), "."
    )
  }
  candidates
}

# Enumerate consecutive exon windows that can receive each partial CNA while
# satisfying coding-residue and structural-localization constraints.
build_legal_cna_windows <- function(
    events,
    exons,
    exon.residue.map,
    max.residue.fraction = 0.999999) {
  eligible.orders <- exons$exon_order[exons$structurally_eligible]
  total.exons <- nrow(exons)
  total.residues <- length(unique(exon.residue.map$residue))
  spans <- sort(unique(normalize_cna_exon_spans(
    events$span_exons,
    total.exons,
    context = "Partial structurally localizable CNA events"
  )))
  windows <- vector("list", length(spans))
  names(windows) <- as.character(spans)

  for (span in spans) {
    starts <- seq_len(total.exons - span + 1L)
    rows <- lapply(starts, function(start) {
      orders <- seq.int(start, length.out = span)
      n.residues <- count_window_residues(orders, exon.residue.map)
      data.frame(
        exon_start = start,
        exon_end = start + span - 1L,
        span_exons = span,
        n_coding_residues = n.residues,
        structurally_eligible = all(orders %in% eligible.orders),
        residue_fraction = n.residues / total.residues,
        stringsAsFactors = FALSE
      )
    })
    x <- do.call(rbind, rows)
    x <- x[
      x$structurally_eligible &
        x$residue_fraction <= max.residue.fraction,
      , drop = FALSE
    ]
    if (!nrow(x)) {
      stop(
        "No legal structurally eligible window can accommodate an event ",
        "spanning ", span, " coding exon(s)."
      )
    }
    windows[[as.character(span)]] <- x
  }
  windows
}

# Assign relocation probabilities using residue-span similarity plus a uniform
# mixture so that alternative windows retain nonzero probability.
cna_window_sampling_probabilities <- function(
    candidates,
    observed.residue.span,
    residue.span.temperature = 0.50,
    uniform.window.mixture = 0.25) {
  if (is.null(candidates)) stop("The legal CNA-window object is NULL.")
  candidates <- as.data.frame(candidates, stringsAsFactors = FALSE)
  if (base::nrow(candidates) == 0L) {
    stop("At least one legal CNA window is required.")
  }
  if (!is.finite(residue.span.temperature) ||
      residue.span.temperature <= 0) {
    stop("residue.span.temperature must be positive.")
  }
  if (!is.finite(uniform.window.mixture) ||
      uniform.window.mixture < 0 || uniform.window.mixture > 1) {
    stop("uniform.window.mixture must be between 0 and 1.")
  }

  # A symmetric log-ratio measures how different each candidate coding span is
  # from the observed span. Unlike v3, no otherwise legal window is discarded.
  log.span.error <- abs(log(
    pmax(candidates$n_coding_residues, 1) /
      max(observed.residue.span, 1)
  ))
  similarity <- exp(-log.span.error / residue.span.temperature)
  if (!all(is.finite(similarity)) || sum(similarity) <= 0) {
    similarity <- rep(1, nrow(candidates))
  }
  similarity <- similarity / sum(similarity)
  uniform <- rep(1 / nrow(candidates), nrow(candidates))
  probability <-
    (1 - uniform.window.mixture) * similarity +
    uniform.window.mixture * uniform
  probability <- probability / sum(probability)

  candidates$log_residue_span_error <- log.span.error
  candidates$sampling_probability <- probability
  candidates
}

# Summarize the number, probability and effective number of legal placements
# available to every event in the partial-CNA null.
build_cna_permutation_diagnostics <- function(
    events,
    legal.windows,
    residue.span.temperature = 0.50,
    uniform.window.mixture = 0.25,
    min.effective.windows = 1.25,
    min.alternative.probability = 0.05) {
  if (!nrow(events)) {
    return(data.frame(
      event_id = character(), subject_id = character(), cna_type = character(),
      observed_exon_start = integer(), observed_exon_end = integer(),
      span_exons = integer(), observed_coding_residues = integer(),
      n_legal_windows = integer(), n_alternative_windows = integer(),
      probability_observed_window = numeric(),
      total_alternative_probability = numeric(),
      maximum_window_probability = numeric(),
      effective_number_of_windows = numeric(), normalized_entropy = numeric(),
      permutation_testable = logical(), testability_reason = character(),
      stringsAsFactors = FALSE
    ))
  }

  rows <- vector("list", nrow(events))
  for (i in seq_len(nrow(events))) {
    event <- events[i, , drop = FALSE]
    candidates <- get_cna_legal_windows(
      legal.windows,
      event$span_exons,
      event.id = as.character(event$event_id[[1L]])
    )
    candidates <- cna_window_sampling_probabilities(
      candidates,
      event$n_coding_residues[[1L]],
      residue.span.temperature,
      uniform.window.mixture
    )
    observed <-
      candidates$exon_start == event$exon_start &
      candidates$exon_end == event$exon_end
    observed.probability <- if (any(observed)) {
      sum(candidates$sampling_probability[observed])
    } else {
      NA_real_
    }
    alternative.probability <- if (is.na(observed.probability)) {
      1
    } else {
      1 - observed.probability
    }
    effective.windows <- 1 / sum(candidates$sampling_probability^2)
    entropy <- -sum(
      candidates$sampling_probability *
        log(candidates$sampling_probability)
    )
    normalized.entropy <- if (nrow(candidates) > 1L) {
      entropy / log(nrow(candidates))
    } else {
      0
    }
    testable <-
      nrow(candidates) > 1L &&
      effective.windows >= min.effective.windows &&
      alternative.probability >= min.alternative.probability
    reason <- if (nrow(candidates) <= 1L) {
      "only_one_legal_window"
    } else if (effective.windows < min.effective.windows) {
      "effective_number_of_windows_too_small"
    } else if (alternative.probability < min.alternative.probability) {
      "insufficient_probability_of_alternative_placement"
    } else {
      ""
    }

    rows[[i]] <- data.frame(
      event_id = event$event_id,
      subject_id = event$subject_id,
      cna_type = event$cna_type,
      observed_exon_start = event$exon_start,
      observed_exon_end = event$exon_end,
      span_exons = event$span_exons,
      observed_coding_residues = event$n_coding_residues,
      n_legal_windows = nrow(candidates),
      n_alternative_windows = nrow(candidates) - as.integer(any(observed)),
      probability_observed_window = observed.probability,
      total_alternative_probability = alternative.probability,
      maximum_window_probability = max(candidates$sampling_probability),
      effective_number_of_windows = effective.windows,
      normalized_entropy = normalized.entropy,
      permutation_testable = testable,
      testability_reason = reason,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# Sample one legal exon window for an event under the partial-CNA null.
sample_cna_window <- function(
    candidates,
    observed.residue.span,
    residue.span.temperature = 0.50,
    uniform.window.mixture = 0.25) {
  candidates <- cna_window_sampling_probabilities(
    candidates,
    observed.residue.span,
    residue.span.temperature,
    uniform.window.mixture
  )
  candidates[
    sample(
      seq_len(nrow(candidates)),
      1L,
      prob = candidates$sampling_probability
    ),
    , drop = FALSE
  ]
}

# Relocate all selected CNA events as intact exon intervals while preserving
# their subjects, lesion types and observed exon-span sizes.
simulate_cna_event_mapping <- function(
    events,
    legal.windows,
    exons,
    residue.span.temperature = 0.50,
    uniform.window.mixture = 0.25,
    avoid.within.subject.overlap = TRUE) {
  simulated <- vector("list", nrow(events))
  used.by.subject <- new.env(parent = emptyenv())

  for (i in seq_len(nrow(events))) {
    event <- events[i, , drop = FALSE]
    candidates <- get_cna_legal_windows(
      legal.windows,
      event$span_exons,
      event.id = as.character(event$event_id[[1L]])
    )
    subject <- as.character(event$subject_id[[1L]])
    used <- if (exists(subject, used.by.subject, inherits = FALSE)) {
      get(subject, used.by.subject, inherits = FALSE)
    } else {
      integer()
    }

    if (avoid.within.subject.overlap && length(used)) {
      nonoverlap <- vapply(seq_len(base::nrow(candidates)), function(j) {
        orders <- seq.int(
          candidates$exon_start[[j]], candidates$exon_end[[j]]
        )
        !any(orders %in% used)
      }, logical(1L))
      if (any(nonoverlap)) candidates <- candidates[nonoverlap, , drop = FALSE]
    }

    selected <- sample_cna_window(
      candidates,
      observed.residue.span = event$n_coding_residues[[1L]],
      residue.span.temperature = residue.span.temperature,
      uniform.window.mixture = uniform.window.mixture
    )
    affected <- seq.int(selected$exon_start, selected$exon_end)
    assign(subject, sort(unique(c(used, affected))), used.by.subject)
    simulated[[i]] <- data.frame(
      event_id = event$event_id,
      subject_id = subject,
      cna_type = event$cna_type,
      exon_start = selected$exon_start,
      exon_end = selected$exon_end,
      span_exons = as.integer(event$span_exons[[1L]]),
      n_coding_residues = selected$n_coding_residues,
      exon_order = affected,
      exon_id = exons$exon_id[match(affected, exons$exon_order)],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, simulated)
}

# Record the smallest candidate diameter found for each affected-subject count.
minimum_cna_diameter_by_size <- function(candidates, sizes, metric) {
  result <- rep(Inf, length(sizes))
  names(result) <- as.character(sizes)
  if (is.null(candidates) || base::nrow(as.data.frame(candidates)) == 0L) {
    return(result)
  }
  for (i in seq_along(sizes)) {
    value <- candidates[candidates$n_subjects == sizes[[i]], metric]
    if (length(value)) result[[i]] <- min(value)
  }
  result
}

# Repeat the complete partial-CNA relocation and 1D/3D tree search to build
# size-specific null distributions of the best cluster diameters.
run_cna_null_simulations <- function(
    events,
    exons,
    exon.residue.map,
    exon.distance.matrix,
    sizes,
    n.sim = 10000L,
    min.subjects,
    min.events,
    max.residue.fraction,
    residue.span.temperature,
    uniform.window.mixture,
    avoid.within.subject.overlap,
    random.seed,
    progress.every = 100L,
    linkage.method = "complete",
    mc.cores = NULL) {
  validate_linkage_method(linkage.method)
  legal.windows <- build_legal_cna_windows(
    events, exons, exon.residue.map, max.residue.fraction
  )

  n.cores <- if (!is.null(mc.cores)) {
    as.integer(mc.cores)
  } else {
    getOption("mc.cores", parallel::detectCores())
  }
  if (is.na(n.cores) || n.cores < 1L) {
    n.cores <- 1L
  }

  can_fork <- (.Platform$OS.type != "windows") &&
    (n.cores > 1L) &&
    exists("mclapply", where = asNamespace("parallel"), mode = "function")

  old.rng <- RNGkind("L'Ecuyer-CMRG")
  on.exit(RNGkind(old.rng[1L], old.rng[2L], old.rng[3L]), add = TRUE)
  set.seed(random.seed)

  run_one_sim <- function(b) {
    simulated.mapping <- simulate_cna_event_mapping(
      events,
      legal.windows,
      exons,
      residue.span.temperature,
      uniform.window.mixture,
      avoid.within.subject.overlap
    )
    simulated.candidates <- build_cna_candidate_clusters(
      simulated.mapping, exons, exon.distance.matrix,
      min.subjects, min.events,
      linkage.method = linkage.method
    )
    list(
      d1 = minimum_cna_diameter_by_size(
        simulated.candidates, sizes, "diameter_1d"
      ),
      d3 = minimum_cna_diameter_by_size(
        simulated.candidates, sizes, "diameter_3d"
      )
    )
  }

  if (can_fork) {
    effective.cores <- min(as.integer(n.cores), as.integer(n.sim))
    sim_results <- parallel::mclapply(
      seq_len(n.sim),
      run_one_sim,
      mc.cores = effective.cores,
      mc.preschedule = TRUE,
      mc.set.seed = TRUE
    )
    failed <- vapply(sim_results, inherits, logical(1L), what = "try-error")
    if (any(failed)) {
      err_msg <- as.character(sim_results[[which(failed)[1L]]])
      stop("Structural simulation worker failed: ", err_msg)
    }
    null.1d <- matrix(
      unlist(lapply(sim_results, `[[`, "d1"), use.names = FALSE),
      nrow = n.sim, byrow = TRUE,
      dimnames = list(NULL, as.character(sizes))
    )
    null.3d <- matrix(
      unlist(lapply(sim_results, `[[`, "d3"), use.names = FALSE),
      nrow = n.sim, byrow = TRUE,
      dimnames = list(NULL, as.character(sizes))
    )
    if (progress.every > 0L) {
      message("Completed ", n.sim, " of ", n.sim, " structural simulations.")
    }
  } else {
    null.1d <- matrix(
      Inf, nrow = n.sim, ncol = length(sizes),
      dimnames = list(NULL, as.character(sizes))
    )
    null.3d <- null.1d
    for (b in seq_len(n.sim)) {
      res <- run_one_sim(b)
      null.1d[b, ] <- res$d1
      null.3d[b, ] <- res$d3

      if (progress.every > 0L &&
          (b %% progress.every == 0L || b == n.sim)) {
        message("Completed ", b, " of ", n.sim, " structural simulations.")
      }
    }
  }

  list(diameter_1d = null.1d, diameter_3d = null.3d)
}


# ==============================================================================
# 8. SIZE-SPECIFIC AND PROTEIN-WIDE CALIBRATION
# ==============================================================================

# Calculate size-specific, protein-wide, joint and omnibus empirical p-values
# for the observed CNA candidate clusters.
calibrate_cna_candidates <- function(
    candidates,
    null.results,
    analysis.type,
    calibration.fraction = 0.50,
    alpha = 0.05) {
  n.sim <- nrow(null.results$diameter_1d)
  n.calibration <- floor(n.sim * calibration.fraction)
  if (n.calibration < 20L || n.sim - n.calibration < 20L) {
    stop(
      "Use at least 40 structural simulations so both calibration and ",
      "protein-wide evaluation contain at least 20 replicates."
    )
  }
  calibration.rows <- seq_len(n.calibration)
  evaluation.rows <- seq.int(n.calibration + 1L, n.sim)
  if (length(calibration.rows) != n.calibration ||
      length(evaluation.rows) != (n.sim - n.calibration)) {
    stop("Calibration and evaluation split partition mismatch.")
  }
  sizes <- as.integer(colnames(null.results$diameter_1d))
  calibration.1d <- null.results$diameter_1d[
    calibration.rows, , drop = FALSE
  ]
  calibration.3d <- null.results$diameter_3d[
    calibration.rows, , drop = FALSE
  ]
  evaluation.1d <- null.results$diameter_1d[
    evaluation.rows, , drop = FALSE
  ]
  evaluation.3d <- null.results$diameter_3d[
    evaluation.rows, , drop = FALSE
  ]

  evaluation.p.1d <- matrix(
    NA_real_, nrow = nrow(evaluation.1d), ncol = ncol(evaluation.1d)
  )
  evaluation.p.3d <- evaluation.p.1d
  sorted.1d <- vector("list", length(sizes))
  sorted.3d <- vector("list", length(sizes))
  for (j in seq_along(sizes)) {
    sorted.1d[[j]] <- sort(calibration.1d[, j])
    sorted.3d[[j]] <- sort(calibration.3d[, j])
    evaluation.p.1d[, j] <- cna_empirical_p_from_reference(
      evaluation.1d[, j], sorted.1d[[j]]
    )
    evaluation.p.3d[, j] <- cna_empirical_p_from_reference(
      evaluation.3d[, j], sorted.3d[[j]]
    )
  }
  evaluation.min.1d <- apply(evaluation.p.1d, 1L, min)
  evaluation.min.3d <- apply(evaluation.p.3d, 1L, min)
  evaluation.min.joint <- pmin(evaluation.min.1d, evaluation.min.3d)

  p.columns <- c(
    "p_1d_size", "p_3d_size", "p_1d_protein", "p_3d_protein",
    "p_1d_joint", "p_3d_joint", "p_any_joint", "p_hotspot_omnibus"
  )
  candidates[p.columns] <- NA_real_
  candidates$n_unique_null_1d_diameters <- 0L
  candidates$n_unique_null_3d_diameters <- 0L
  candidates$null_1d_variable <- FALSE
  candidates$null_3d_variable <- FALSE
  for (i in seq_len(nrow(candidates))) {
    j <- match(candidates$n_subjects[[i]], sizes)
    candidates$n_unique_null_1d_diameters[[i]] <-
      length(unique(calibration.1d[, j]))
    candidates$n_unique_null_3d_diameters[[i]] <-
      length(unique(calibration.3d[, j]))
    candidates$null_1d_variable[[i]] <-
      candidates$n_unique_null_1d_diameters[[i]] > 1L
    candidates$null_3d_variable[[i]] <-
      candidates$n_unique_null_3d_diameters[[i]] > 1L
    p.1d <- cna_empirical_lower_p(
      candidates$diameter_1d[[i]], calibration.1d[, j]
    )
    p.3d <- cna_empirical_lower_p(
      candidates$diameter_3d[[i]], calibration.3d[, j]
    )
    candidates$p_1d_size[[i]] <- p.1d
    candidates$p_3d_size[[i]] <- p.3d
    candidates$p_1d_protein[[i]] <-
      (1 + sum(evaluation.min.1d <= p.1d)) /
      (length(evaluation.min.1d) + 1)
    candidates$p_3d_protein[[i]] <-
      (1 + sum(evaluation.min.3d <= p.3d)) /
      (length(evaluation.min.3d) + 1)
    candidates$p_1d_joint[[i]] <-
      (1 + sum(evaluation.min.joint <= p.1d)) /
      (length(evaluation.min.joint) + 1)
    candidates$p_3d_joint[[i]] <-
      (1 + sum(evaluation.min.joint <= p.3d)) /
      (length(evaluation.min.joint) + 1)
    candidates$p_any_joint[[i]] <-
      (1 + sum(evaluation.min.joint <= min(p.1d, p.3d))) /
      (length(evaluation.min.joint) + 1)
    candidates$p_hotspot_omnibus[[i]] <- candidates$p_any_joint[[i]]
  }

  single.exon <- candidates$n_exons < 2L
  candidates$compactness_testable <-
    !single.exon &
    (candidates$null_1d_variable | candidates$null_3d_variable)
  candidates$compactness_testability_reason <- ifelse(
    single.exon,
    "single_exon_cluster",
    ifelse(
      !candidates$null_1d_variable & !candidates$null_3d_variable,
      "invariant_1d_and_3d_null",
      ifelse(
        !candidates$null_1d_variable,
        "invariant_1d_null",
        ifelse(
          !candidates$null_3d_variable,
          "invariant_3d_null",
          ""
        )
      )
    )
  )
  candidates[single.exon, p.columns] <- NA_real_
  candidates[
    !single.exon & !candidates$null_1d_variable,
    c("p_1d_size", "p_1d_protein", "p_1d_joint")
  ] <- NA_real_
  candidates[
    !single.exon & !candidates$null_3d_variable,
    c("p_3d_size", "p_3d_protein", "p_3d_joint")
  ] <- NA_real_
  candidates[
    !single.exon & !candidates$null_1d_variable &
      !candidates$null_3d_variable,
    c("p_any_joint", "p_hotspot_omnibus")
  ] <- NA_real_
  candidates$significant_1d <-
    !single.exon & candidates$null_1d_variable &
      !is.na(candidates$p_1d_joint) & candidates$p_1d_joint <= alpha
  candidates$significant_3d <-
    !single.exon & candidates$null_3d_variable &
      !is.na(candidates$p_3d_joint) & candidates$p_3d_joint <= alpha
  candidates$significant_any <-
    candidates$compactness_testable &
      !is.na(candidates$p_hotspot_omnibus) &
      candidates$p_hotspot_omnibus <= alpha
  candidates$hotspot_class <- ifelse(
    single.exon,
    "Single-exon recurrence assessed separately",
    ifelse(
      !candidates$compactness_testable,
      "Multi-exon compactness not testable: invariant null",
      ifelse(
        candidates$significant_1d & candidates$significant_3d,
        "Exon-local CNA hotspot also compact in structure",
        ifelse(
          !candidates$significant_1d & candidates$significant_3d,
          "Nonlocal 3D-specific exon CNA hotspot",
          ifelse(
            candidates$significant_1d & !candidates$significant_3d,
            "Exon-order hotspot without independent 3D enrichment",
            "Weak clustering evidence"
          )
        )
      )
    )
  )
  candidates <- candidates[
    order(
      candidates$p_hotspot_omnibus, candidates$p_3d_joint,
      candidates$p_1d_joint, candidates$diameter_3d
    ), , drop = FALSE
  ]
  safe.type <- safe_grin3d_name(analysis.type)
  candidates$cluster_id <- sprintf(
    "%s_cluster_%04d", safe.type, seq_len(nrow(candidates))
  )
  candidates$analysis_type <- analysis.type
  rownames(candidates) <- NULL

  list(
    clusters = candidates,
    evaluation_min_p = data.frame(
      analysis_type = analysis.type,
      min_p_1d = evaluation.min.1d,
      min_p_3d = evaluation.min.3d,
      min_p_joint = evaluation.min.joint
    ),
    n_calibration = n.calibration,
    n_evaluation = length(evaluation.rows)
  )
}

# Expand the cluster table into one row per cluster-exon member.
build_cna_cluster_member_table <- function(clusters, mapping) {
  if (!nrow(clusters)) return(data.frame())
  rows <- vector("list", nrow(clusters))
  for (i in seq_len(nrow(clusters))) {
    exon.orders <- as.integer(strsplit(
      clusters$exon_orders[[i]], ";", fixed = TRUE
    )[[1L]])
    z <- mapping[mapping$exon_order %in% exon.orders, , drop = FALSE]
    event.rows <- split(z, z$event_id)
    rows[[i]] <- do.call(rbind, lapply(event.rows, function(one) {
      n.exons.hit <- length(unique(one$exon_order))
      data.frame(
        cluster_id = clusters$cluster_id[[i]],
        analysis_type = clusters$analysis_type[[i]],
        cna_type = one$cna_type[[1L]],
        supports_analysis_type =
          clusters$analysis_type[[i]] == "ALL_CNA" ||
          one$cna_type[[1L]] == clusters$analysis_type[[i]],
        event_id = one$event_id[[1L]],
        subject_id = one$subject_id[[1L]],
        event_coding_fraction = one$coding_fraction[[1L]],
        event_localization_class = one$localization_class[[1L]],
        n_cluster_exons_hit_by_event = n.exons.hit,
        spans_multiple_cluster_exons = n.exons.hit >= 2L,
        affected_exon_orders_in_cluster = collapse_grin3d_values(
          one$exon_order
        ),
        affected_transcript_exon_ranks_in_cluster = collapse_grin3d_values(
          one$transcript_exon_rank
        ),
        affected_exon_ids_in_cluster = collapse_grin3d_values(one$exon_id),
        stringsAsFactors = FALSE
      )
    }))
  }
  out <- do.call(rbind, rows)
  out <- out[order(out$cluster_id, out$subject_id, out$event_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Create an empty cluster table with the same columns as a completed analysis.
empty_cna_cluster_table <- function(cna.types) {
  out <- data.frame(
    protein = character(), transcript = character(),
    cluster_id = character(), analysis_type = character(),
    candidate_key = character(), tree_sources = character(),
    exon_orders = character(), transcript_exon_ranks = character(),
    exon_ids = character(), exon_order_min = integer(),
    exon_order_max = integer(), n_exons = integer(),
    n_events = integer(), n_subjects = integer(),
    diameter_1d = numeric(), diameter_3d = numeric(),
    events_per_exon = character(), subjects_per_exon = character(),
    min_events_per_exon = integer(), min_subjects_per_exon = integer(),
    n_events_spanning_multiple_cluster_exons = integer(),
    n_subjects_spanning_multiple_cluster_exons = integer(),
    n_subjects_with_multiple_events_in_cluster = integer(),
    best_supported_subcluster_key = character(),
    n_events_added_vs_best_subcluster = integer(),
    n_subjects_added_vs_best_subcluster = integer(),
    branch_adds_independent_subjects = logical(),
    best_recurrent_exon_order = integer(),
    best_recurrent_transcript_exon_rank = integer(),
    best_exon_p_subjects_raw = numeric(),
    recurrence_scope_used = character(),
    best_exon_p_subjects_scope_omnibus = numeric(),
    best_exon_p_subjects_across_scopes = numeric(),
    best_exon_p_subjects_omnibus = numeric(),
    n_testable_recurrence_exons = integer(),
    n_significant_recurrent_exons = integer(),
    integrated_evidence_class = character(),
    p_1d_size = numeric(), p_3d_size = numeric(),
    p_1d_protein = numeric(), p_3d_protein = numeric(),
    p_1d_joint = numeric(), p_3d_joint = numeric(),
    p_any_joint = numeric(), p_hotspot_omnibus = numeric(),
    n_unique_null_1d_diameters = integer(),
    n_unique_null_3d_diameters = integer(),
    null_1d_variable = logical(), null_3d_variable = logical(),
    compactness_testable = logical(),
    compactness_testability_reason = character(),
    significant_1d = logical(), significant_3d = logical(),
    significant_any = logical(), hotspot_class = character(),
    stringsAsFactors = FALSE
  )
  for (type in c(cna.types, "ALL_CNA")) {
    safe <- safe_grin3d_name(type)
    out[[paste0("n_events_", safe)]] <- integer()
    out[[paste0("n_subjects_", safe)]] <- integer()
  }
  out
}

# Create an empty cluster-member table with the expected output columns.
empty_cna_member_table <- function() {
  data.frame(
    cluster_id = character(), analysis_type = character(),
    cna_type = character(), supports_analysis_type = logical(),
    event_id = character(), subject_id = character(),
    event_coding_fraction = numeric(),
    event_localization_class = character(),
    n_cluster_exons_hit_by_event = integer(),
    spans_multiple_cluster_exons = logical(),
    affected_exon_orders_in_cluster = character(),
    affected_transcript_exon_ranks_in_cluster = character(),
    affected_exon_ids_in_cluster = character(),
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 9. USER-FACING RUNNER
# ==============================================================================

# Run the complete exon-level CNA analysis, write all result files and return
# the full result object invisibly for interactive use in R.
run_grin3d_exon_cna_hotspots <- function(
    cna.file,
    exon.file,
    coordinate.file,
    results.dir = "results_exon_CNA_hotspots",
    protein = "protein_of_interest",
    transcript = "selected_transcript",
    analysis.types = NULL,
    include.combined.analysis = TRUE,
    cna.columns = list(
      subject = "ID",
      event = NULL,
      chrom = "chrom",
      start = "loc.start",
      end = "loc.end",
      type = "lsn.type"
    ),
    exon.columns = list(
      gene_id = "gene_id",
      gene_name = "gene_name",
      transcript_id = "transcript_id",
      protein_id = "protein_id",
      exon_id = "exon_id",
      exon_rank = "exon_rank",
      chrom = "chrom",
      strand = "strand",
      genomic_start = "cds_genomic_start",
      genomic_end = "cds_genomic_end",
      residue_start = "aa_start",
      residue_end = "aa_end"
    ),
    coordinate.columns = list(
      residue = "residue",
      x = "x",
      y = "y",
      z = "z",
      confidence = NULL
    ),
    cna.object = NULL,
    exon.object = NULL,
    coordinate.object = NULL,
    coverage.thresholds = c(0.50, 0.75, 0.90, 1.00),
    chromosome.length = NULL,
    min.exon.coding.overlap = 0,
    structural.max.coding.fraction = 0.999999,
    min.confidence = 70,
    min.exon.mapped.fraction = 0.50,
    min.exon.mapped.residues = 1L,
    n.closest.cross.exon.distances = 5L,
    min.sequence.separation = 10L,
    residue.span.temperature = 0.50,
    uniform.window.mixture = 0.25,
    min.effective.windows = 1.25,
    min.alternative.probability = 0.05,
    avoid.within.subject.overlap = TRUE,
    min.subjects = 2L,
    min.events = 2L,
    alpha = 0.05,
    n.sim.coverage = 10000L,
    n.sim.structural = 10000L,
    calibration.fraction = 0.50,
    random.seed = 20260907L,
    progress.every = 100L,
    linkage.method = "complete",
    mc.cores = NULL) {
  validate_linkage_method(linkage.method)
  validate_grin3d_column_map(
    cna.columns,
    required.keys = c("subject", "chrom", "start", "end", "type"),
    optional.keys = "event",
    map.name = "cna.columns"
  )
  validate_grin3d_column_map(
    exon.columns,
    required.keys = c(
      "exon_id", "exon_rank", "chrom", "strand", "genomic_start",
      "genomic_end", "residue_start", "residue_end"
    ),
    optional.keys = c("gene_id", "gene_name", "transcript_id", "protein_id"),
    map.name = "exon.columns"
  )
  validate_grin3d_column_map(
    coordinate.columns,
    required.keys = c("residue", "x", "y", "z"),
    optional.keys = "confidence",
    map.name = "coordinate.columns"
  )

  dir.create(results.dir, recursive = TRUE, showWarnings = FALSE)
  message("Reading CNA, coding-exon and AlphaFold coordinate files.")
  cna.data <- read_grin3d_input(cna.file, cna.object)
  exon.data <- read_grin3d_input(exon.file, exon.object)
  coordinate.data <- read_grin3d_input(coordinate.file, coordinate.object)

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

  coordinates <- prepare_cna_coordinates(
    coordinate.data, coordinate.columns, min.confidence
  )
  prepared.exons <- prepare_coding_exons(
    exon.data, exon.columns, coordinates,
    min.exon.mapped.fraction, min.exon.mapped.residues
  )
  exons <- prepared.exons$exons
  exon.residue.map <- prepared.exons$exon_residue_map
  if (any(!exons$structurally_eligible)) {
    failed <- exons[!exons$structurally_eligible, , drop = FALSE]
    warning(
      "Structurally ineligible coding exon(s): ",
      paste(
        paste0(
          "rank ", failed$transcript_exon_rank,
          " [", failed$n_structural_residues, "/",
          failed$n_coding_residues, " mapped; ",
          failed$structural_eligibility_reason, "]"
        ),
        collapse = "; "
      ),
      ". See CNA_exon_counts_by_type.csv for diagnostics."
    )
  }

  raw.events <- prepare_genomic_cna_events(cna.data, cna.columns)
  mapped <- map_genomic_cnas_to_exons(
    raw.events, exons, exon.residue.map, min.exon.coding.overlap
  )
  if (!nrow(mapped$events)) {
    stop("No CNA event overlaps the selected transcript's coding sequence.")
  }

  observed.types <- sort(unique(mapped$events$cna_type))
  if (is.null(analysis.types)) {
    analysis.types <- observed.types
  } else {
    analysis.types <- unique(toupper(as.character(analysis.types)))
    absent <- setdiff(analysis.types, observed.types)
    if (length(absent)) {
      warning(
        "Ignoring CNA type(s) absent from this gene: ",
        paste(absent, collapse = ", ")
      )
    }
    analysis.types <- intersect(analysis.types, observed.types)
  }
  if (!length(analysis.types)) stop("No requested CNA types overlap the gene.")

  mapped$events <- mapped$events[
    mapped$events$cna_type %in% analysis.types, , drop = FALSE
  ]
  mapped$mapping <- mapped$mapping[
    mapped$mapping$cna_type %in% analysis.types, , drop = FALSE
  ]

  # The partial-localizable scope is fixed before recurrence calibration and
  # is the same event set used by the 1D/3D exon-tree analysis.
  structural.events <- mapped$events[
    mapped$events$structurally_localizable &
      mapped$events$coding_fraction <= structural.max.coding.fraction,
    , drop = FALSE
  ]
  structural.events$span_exons <- structural.events$n_affected_exons
  structural.mapping <- mapped$mapping[
    mapped$mapping$event_id %in% structural.events$event_id &
      mapped$mapping$passes_exon_overlap,
    , drop = FALSE
  ]

  message(
    "Calculating conditional coding-coverage and ALL_OVERLAP exon-recurrence ",
    "p-values."
  )
  coverage <- run_cna_coverage_null(
    events = mapped$events,
    exons = exons,
    cna.types = analysis.types,
    thresholds = coverage.thresholds,
    n.sim = n.sim.coverage,
    chromosome.length = chromosome.length,
    min.exon.coding.overlap = min.exon.coding.overlap,
    include.combined.analysis = include.combined.analysis,
    alpha = alpha,
    random.seed = random.seed,
    progress.every = progress.every
  )

  message(
    "Calculating PARTIAL_LOCALIZABLE exon-recurrence p-values."
  )
  partial.recurrence <- run_partial_localizable_recurrence_null(
    events = structural.events,
    observed.mapping = structural.mapping,
    exons = exons,
    exon.residue.map = exon.residue.map,
    cna.types = analysis.types,
    include.combined.analysis = include.combined.analysis,
    n.sim = n.sim.coverage,
    max.residue.fraction = structural.max.coding.fraction,
    residue.span.temperature = residue.span.temperature,
    uniform.window.mixture = uniform.window.mixture,
    min.effective.windows = min.effective.windows,
    min.alternative.probability = min.alternative.probability,
    avoid.within.subject.overlap = avoid.within.subject.overlap,
    alpha = alpha,
    random.seed = random.seed + 1000L,
    progress.every = progress.every
  )
  exon.recurrence.by.scope <- list(
    ALL_OVERLAP = coverage$exon_recurrence,
    PARTIAL_LOCALIZABLE = partial.recurrence
  )
  exon.recurrence.results <- combine_exon_recurrence_scopes(
    exon.recurrence.by.scope,
    alpha
  )

  message("Calculating robust exon-to-exon C-alpha distances.")
  structural <- build_exon_distance_matrix(
    exons, exon.residue.map, coordinates,
    n.closest.cross.exon.distances, min.sequence.separation
  )
  structural.global.reason <- if (
    nrow(structural$eligible_exons) < 2L
  ) {
    "fewer_than_two_structurally_eligible_exons"
  } else {
    NULL
  }

  exon.counts <- build_exon_type_counts(
    exons,
    mapped$mapping,
    analysis.types,
    structural.mapping = structural.mapping
  )

  analysis.labels <- analysis.types
  if (include.combined.analysis && length(analysis.types) > 1L) {
    analysis.labels <- c(analysis.labels, "ALL_CNA")
  }
  cluster.results <- list()
  member.results <- list()
  null.results <- list()
  null.min.p <- list()
  testability <- list()

  for (a in seq_along(analysis.labels)) {
    label <- analysis.labels[[a]]
    selected.events <- if (label == "ALL_CNA") {
      structural.events
    } else {
      structural.events[
        structural.events$cna_type == label, , drop = FALSE
      ]
    }
    selected.mapping <- structural.mapping[
      structural.mapping$event_id %in% selected.events$event_id,
      , drop = FALSE
    ]

    mapped.events.for.label <- if (label == "ALL_CNA") {
      mapped$events
    } else {
      mapped$events[mapped$events$cna_type == label, , drop = FALSE]
    }
    n.coding.mapped.subjects <- length(unique(
      mapped.events.for.label$subject_id
    ))
    n.coding.mapped.events <- length(unique(mapped.events.for.label$event_id))
    n.subjects <- length(unique(selected.events$subject_id))
    n.events <- length(unique(selected.events$event_id))
    localizable.fraction <- if (n.coding.mapped.events) {
      n.events / n.coding.mapped.events
    } else {
      NA_real_
    }
    if (!is.null(structural.global.reason)) {
      testability[[a]] <- data.frame(
        analysis_type = label,
        n_coding_mapped_subjects = n.coding.mapped.subjects,
        n_coding_mapped_events = n.coding.mapped.events,
        n_structurally_localizable_subjects = n.subjects,
        n_structurally_localizable_events = n.events,
        structurally_localizable_event_fraction = localizable.fraction,
        testable = FALSE,
        reason = structural.global.reason,
        stringsAsFactors = FALSE
      )
      next
    }
    if (n.subjects < min.subjects || n.events < min.events) {
      testability[[a]] <- data.frame(
        analysis_type = label,
        n_coding_mapped_subjects = n.coding.mapped.subjects,
        n_coding_mapped_events = n.coding.mapped.events,
        n_structurally_localizable_subjects = n.subjects,
        n_structurally_localizable_events = n.events,
        structurally_localizable_event_fraction = localizable.fraction,
        testable = FALSE,
        reason = "fewer_than_minimum_subjects_or_events",
        stringsAsFactors = FALSE
      )
      next
    }

    candidates <- build_cna_candidate_clusters(
      selected.mapping, exons, structural$distance_matrix,
      min.subjects, min.events,
      linkage.method = linkage.method
    )
    if (is.null(candidates) || base::nrow(candidates) == 0L) {
      testability[[a]] <- data.frame(
        analysis_type = label,
        n_coding_mapped_subjects = n.coding.mapped.subjects,
        n_coding_mapped_events = n.coding.mapped.events,
        n_structurally_localizable_subjects = n.subjects,
        n_structurally_localizable_events = n.events,
        structurally_localizable_event_fraction = localizable.fraction,
        testable = FALSE,
        reason = "no_candidate_cluster_met_testing_criteria",
        stringsAsFactors = FALSE
      )
      next
    }

    message(
      "Running ", n.sim.structural, " structural simulations for ", label, "."
    )
    sizes <- sort(unique(candidates$n_subjects))
    one.null <- run_cna_null_simulations(
      events = selected.events,
      exons = exons,
      exon.residue.map = exon.residue.map,
      exon.distance.matrix = structural$distance_matrix,
      sizes = sizes,
      n.sim = n.sim.structural,
      min.subjects = min.subjects,
      min.events = min.events,
      max.residue.fraction = structural.max.coding.fraction,
      residue.span.temperature = residue.span.temperature,
      uniform.window.mixture = uniform.window.mixture,
      avoid.within.subject.overlap = avoid.within.subject.overlap,
      random.seed = random.seed + a,
      progress.every = progress.every,
      linkage.method = linkage.method,
      mc.cores = mc.cores
    )
    calibrated <- calibrate_cna_candidates(
      candidates, one.null, label, calibration.fraction, alpha
    )
    calibrated$clusters <- add_cluster_support_metrics(
      calibrated$clusters, selected.mapping
    )
    calibrated$clusters <- annotate_clusters_with_exon_recurrence(
      calibrated$clusters,
      exon.recurrence.results[
        exon.recurrence.results$recurrence_scope ==
          "PARTIAL_LOCALIZABLE",
        , drop = FALSE
      ],
      alpha
    )
    clusters <- add_cluster_type_counts(
      calibrated$clusters, structural.mapping, analysis.types
    )
    clusters$protein <- protein
    clusters$transcript <- transcript
    clusters <- clusters[, c(
      "protein", "transcript", "cluster_id", "analysis_type",
      setdiff(
        names(clusters),
        c("protein", "transcript", "cluster_id", "analysis_type")
      )
    )]

    cluster.results[[a]] <- clusters
    member.results[[a]] <- build_cna_cluster_member_table(
      clusters, structural.mapping
    )
    null.results[[label]] <- one.null
    null.min.p[[a]] <- calibrated$evaluation_min_p
    testability[[a]] <- data.frame(
      analysis_type = label,
      n_coding_mapped_subjects = n.coding.mapped.subjects,
      n_coding_mapped_events = n.coding.mapped.events,
      n_structurally_localizable_subjects = n.subjects,
      n_structurally_localizable_events = n.events,
      structurally_localizable_event_fraction = localizable.fraction,
      testable = TRUE,
      reason = "",
      stringsAsFactors = FALSE
    )
  }

  nonempty.clusters <- cluster.results[
    !vapply(cluster.results, is.null, logical(1L))
  ]
  nonempty.members <- member.results[
    !vapply(member.results, is.null, logical(1L))
  ]
  clusters <- if (length(nonempty.clusters)) {
    do.call(rbind, nonempty.clusters)
  } else {
    empty_cna_cluster_table(analysis.types)
  }
  members <- if (length(nonempty.members)) {
    do.call(rbind, nonempty.members)
  } else {
    empty_cna_member_table()
  }
  testability <- rbind_grin3d_rows(testability)
  nonempty.min.p <- null.min.p[!vapply(null.min.p, is.null, logical(1L))]
  null.min.p <- if (length(nonempty.min.p)) {
    do.call(rbind, nonempty.min.p)
  } else {
    data.frame(
      analysis_type = character(), min_p_1d = numeric(),
      min_p_3d = numeric(), min_p_joint = numeric()
    )
  }

  settings <- list(
    cna.file = cna.file,
    exon.file = exon.file,
    coordinate.file = coordinate.file,
    results.dir = results.dir,
    protein = protein,
    transcript = transcript,
    analysis.types = analysis.types,
    recurrence.scopes = c("ALL_OVERLAP", "PARTIAL_LOCALIZABLE"),
    coverage.thresholds = coverage.thresholds,
    chromosome.length = chromosome.length,
    min.exon.coding.overlap = min.exon.coding.overlap,
    structural.max.coding.fraction = structural.max.coding.fraction,
    confidence.column = coordinate.columns$confidence,
    confidence.filter.applied = !is.null(coordinate.columns$confidence),
    min.confidence = min.confidence,
    min.exon.mapped.fraction = min.exon.mapped.fraction,
    min.exon.mapped.residues = min.exon.mapped.residues,
    n.closest.cross.exon.distances = n.closest.cross.exon.distances,
    min.sequence.separation = min.sequence.separation,
    residue.span.temperature = residue.span.temperature,
    uniform.window.mixture = uniform.window.mixture,
    min.effective.windows = min.effective.windows,
    min.alternative.probability = min.alternative.probability,
    n.sim.coverage = n.sim.coverage,
    n.sim.structural = n.sim.structural,
    calibration.fraction = calibration.fraction,
    random.seed = random.seed,
    linkage.method = linkage.method,
    mc.cores = mc.cores
  )
  result <- list(
    settings = settings,
    coordinates = coordinates,
    exons = exons,
    exon_residue_mapping = exon.residue.map,
    all_input_events = mapped$all_events,
    mapped_events = mapped$events,
    unmapped_events = mapped$unmapped_events,
    event_exon_mapping = mapped$mapping,
    exon_counts_by_type = exon.counts,
    gene_coverage_results = coverage$results,
    coverage_null_subject_counts = coverage$null_subject_counts,
    coverage_null_event_counts = coverage$null_event_counts,
    exon_recurrence_results = exon.recurrence.results,
    exon_recurrence_results_all_overlap =
      coverage$exon_recurrence$results,
    exon_recurrence_results_partial_localizable =
      partial.recurrence$results,
    exon_recurrence_null_by_scope = exon.recurrence.by.scope,
    # The following four fields retain the v2 ALL_OVERLAP meaning for
    # backward compatibility.
    exon_recurrence_null_subject_counts =
      coverage$exon_recurrence$null_subject_counts,
    exon_recurrence_null_event_counts =
      coverage$exon_recurrence$null_event_counts,
    exon_recurrence_null_min_subject_p =
      coverage$exon_recurrence$null_min_subject_p,
    exon_recurrence_null_min_event_p =
      coverage$exon_recurrence$null_min_event_p,
    partial_exon_recurrence_null_subject_counts =
      partial.recurrence$null_subject_counts,
    partial_exon_recurrence_null_event_counts =
      partial.recurrence$null_event_counts,
    partial_exon_recurrence_null_min_subject_p =
      partial.recurrence$null_min_subject_p,
    partial_exon_recurrence_null_min_event_p =
      partial.recurrence$null_min_event_p,
    partial_recurrence_permutation_diagnostics =
      partial.recurrence$permutation_diagnostics,
    exon_distance_matrix = structural$distance_matrix,
    exon_pair_distances = structural$pair_table,
    structural_events = structural.events,
    structural_testability = testability,
    linkage_method = linkage.method,
    clusters = clusters,
    cluster_members = members,
    structural_null_diameters = null.results,
    structural_null_min_p = null.min.p
  )

  utils::write.csv(
    result$mapped_events,
    file.path(results.dir, "CNA_event_mapping_summary.csv"),
    row.names = FALSE, na = ""
  )
  utils::write.csv(
    result$event_exon_mapping,
    file.path(results.dir, "CNA_event_exon_membership.csv"),
    row.names = FALSE, na = ""
  )
  utils::write.csv(
    result$unmapped_events,
    file.path(results.dir, "CNA_unmapped_events.csv"),
    row.names = FALSE, na = ""
  )
  utils::write.csv(
    result$exon_counts_by_type,
    file.path(results.dir, "CNA_exon_counts_by_type.csv"),
    row.names = FALSE, na = ""
  )
  utils::write.csv(
    result$gene_coverage_results,
    file.path(results.dir, "CNA_gene_coverage_results.csv"),
    row.names = FALSE, na = ""
  )
  utils::write.csv(
    result$exon_recurrence_results,
    file.path(results.dir, "CNA_exon_recurrence_results.csv"),
    row.names = FALSE, na = ""
  )
  utils::write.csv(
    result$partial_recurrence_permutation_diagnostics,
    file.path(
      results.dir,
      "CNA_partial_recurrence_permutation_diagnostics.csv"
    ),
    row.names = FALSE, na = ""
  )
  utils::write.csv(
    result$exon_pair_distances,
    file.path(results.dir, "CNA_exon_structural_distances.csv"),
    row.names = FALSE, na = ""
  )
  utils::write.csv(
    result$structural_testability,
    file.path(results.dir, "CNA_structural_testability.csv"),
    row.names = FALSE, na = ""
  )
  utils::write.csv(
    result$clusters,
    file.path(results.dir, "CNA_exon_hotspot_clusters.csv"),
    row.names = FALSE, na = ""
  )
  utils::write.csv(
    result$cluster_members,
    file.path(results.dir, "CNA_exon_hotspot_cluster_members.csv"),
    row.names = FALSE, na = ""
  )
  saveRDS(result, file.path(results.dir, "GRIN3D_exon_CNA_results.rds"))

  message("GRIN3D exon/CNA analysis complete: ", results.dir)
  message(
    "Use 100 simulations only for a smoke test; increase to 10,000 or more ",
    "for stable development results and 100,000 for final inference."
  )
  invisible(result)
}
