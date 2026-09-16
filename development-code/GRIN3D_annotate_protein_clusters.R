# GRIN3D: biological annotation of mutation-hotspot clusters.
#
# This module runs after hotspot inference. It never changes cluster membership,
# p-values, or significance. It retrieves source-attributed protein features,
# maps cluster residues by exact UniProt sequence coordinates, writes normalized
# CSV files, and can generate a browsable HTML/SVG report for every cluster.
#
# Automated sources:
#   UniProtKB REST: https://rest.uniprot.org/
#   InterPro REST:  https://www.ebi.ac.uk/interpro/api/
#   PDBe-KB graph:   https://www.ebi.ac.uk/pdbe/graph-api/
#   ChEMBL REST:     https://www.ebi.ac.uk/chembl/api/data/
#
# Required R package: jsonlite. Plotting and HTML generation use base R only.

`%||%` <- function(x, fallback) {
  if (is.null(x) || length(x) == 0L) fallback else x
}

grin3d_annotation_urls <- list(
  uniprot.api = "https://rest.uniprot.org/uniprotkb/%s.json",
  uniprot.web = "https://www.uniprot.org/uniprotkb/%s/entry",
  interpro.matches = paste0(
    "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/",
    "%s/?page_size=200"
  ),
  interpro.entry = "https://www.ebi.ac.uk/interpro/api/entry/interpro/%s/",
  interpro.web = "https://www.ebi.ac.uk/interpro/entry/InterPro/%s/",
  pdbe.interfaces = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot/interface_residues/%s",
  pdbe.ligands = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot/ligand_sites/%s",
  pdbe.web = "https://www.ebi.ac.uk/pdbe/pdbe-kb/proteins/%s",
  chembl.targets = paste0(
    "https://www.ebi.ac.uk/chembl/api/data/target.json?",
    "target_components__accession=%s&limit=100"
  ),
  chembl.mechanisms = paste0(
    "https://www.ebi.ac.uk/chembl/api/data/mechanism.json?",
    "target_chembl_id__in=%s&limit=1000"
  ),
  chembl.molecules = paste0(
    "https://www.ebi.ac.uk/chembl/api/data/molecule.json?",
    "molecule_chembl_id__in=%s&limit=1000"
  ),
  chembl.target.web = "https://www.ebi.ac.uk/chembl/explore/target/%s",
  chembl.molecule.web = "https://www.ebi.ac.uk/chembl/explore/compound/%s"
)

grin3d_uniprot_feature_types <- c(
  "Active site", "Binding site", "Coiled coil", "Compositional bias",
  "Cross-link", "DNA binding", "Disulfide bond", "Domain",
  "Glycosylation", "Lipidation", "Modified residue", "Motif", "Region",
  "Repeat", "Site", "Topological domain", "Transmembrane", "Zinc finger"
)

grin3d_cluster_columns <- c(
  "protein", "cluster_id", "candidate_key", "tree_sources", "hotspot_class",
  "significant_1d", "significant_3d", "significant_any", "n_subjects",
  "n_events", "n_residues", "residues", "p_1d_joint", "p_3d_joint",
  "p_any_joint"
)

grin3d_feature_columns <- c(
  "feature_id", "protein", "uniprot_accession", "source",
  "source_record_id", "annotation_class", "feature_category",
  "feature_subcategory", "classification_rule", "feature_type", "feature_name",
  "description", "start", "end", "feature_segments", "start_modifier",
  "end_modifier", "evidence_codes", "evidence_sources", "evidence_ids",
  "cross_references", "go_terms", "source_url", "source_record_version",
  "source_annotation_date", "retrieved_at_utc"
)

grin3d_publication_columns <- c(
  "feature_id", "protein", "uniprot_accession", "source",
  "source_record_id", "evidence_scope", "evidence_code",
  "publication_id_type", "publication_id", "pubmed_id", "doi", "title",
  "journal", "year", "authors", "publication_url"
)

grin3d_source_columns <- c(
  "source", "source_record_id", "protein", "uniprot_accession", "api_url",
  "web_url", "record_version", "annotation_date", "retrieved_at_utc",
  "raw_response_file", "notes"
)

grin3d_pdbe_ligand_columns <- c(
  "protein", "uniprot_accession", "ligand_code", "ligand_name",
  "contact_residues", "pdb_ids", "chembl_id", "drugbank_id",
  "cofactor_id", "reaction_id", "is_solvent", "evidence_scope",
  "source_url", "retrieved_at_utc"
)

grin3d_chembl_target_columns <- c(
  "protein", "uniprot_accession", "target_chembl_id", "target_name",
  "target_type", "organism", "tax_id", "component_accessions",
  "component_relationships", "evidence_scope", "source_url",
  "retrieved_at_utc"
)

grin3d_chembl_mechanism_columns <- c(
  "protein", "uniprot_accession", "target_chembl_id", "target_name",
  "target_type", "molecule_chembl_id", "molecule_name", "molecule_type",
  "max_phase", "first_approval", "action_type", "mechanism_of_action",
  "mechanism_comment", "direct_interaction", "molecular_mechanism",
  "evidence_scope", "target_url", "molecule_url", "retrieved_at_utc"
)

grin3d_empty_table <- function(columns) {
  out <- as.data.frame(
    setNames(replicate(length(columns), character(), simplify = FALSE), columns),
    stringsAsFactors = FALSE
  )
  out
}

grin3d_as_text <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[[1L]])) "" else as.character(x[[1L]])
}

grin3d_join <- function(x, sort.values = TRUE) {
  if (is.null(x)) return("")
  x <- trimws(as.character(unlist(x, recursive = TRUE, use.names = FALSE)))
  x <- unique(x[nzchar(x) & !is.na(x)])
  if (sort.values) x <- sort(x)
  paste(x, collapse = ";")
}

grin3d_html_text <- function(x) {
  x <- grin3d_as_text(x)
  x <- gsub("<[^>]+>", " ", x)
  x <- gsub("&nbsp;", " ", x, fixed = TRUE)
  x <- gsub("&amp;", "&", x, fixed = TRUE)
  x <- gsub("&lt;", "<", x, fixed = TRUE)
  x <- gsub("&gt;", ">", x, fixed = TRUE)
  trimws(gsub("[[:space:]]+", " ", x))
}

grin3d_html_escape <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

grin3d_safe_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
}

grin3d_stable_id <- function(...) {
  bytes <- utf8ToInt(paste(..., collapse = "|"))
  hash <- 0
  for (byte in bytes) hash <- (hash * 131 + byte) %% 2147483647
  sprintf("%08x", as.integer(hash))
}

grin3d_normalize_doi <- function(x) {
  x <- grin3d_as_text(x)
  x <- sub("^https?://(dx\\.)?doi\\.org/", "", x, ignore.case = TRUE)
  if (grepl("^10\\.[0-9]{4,9}/[^[:space:]]+$", x, ignore.case = TRUE)) x else ""
}

grin3d_parse_residues <- function(x) {
  values <- trimws(unlist(strsplit(as.character(x), "[;,]")))
  values <- values[nzchar(values)]
  residues <- suppressWarnings(as.integer(values))
  if (!length(residues) || anyNA(residues) || any(residues < 1L)) {
    stop("Residues must be positive integers: ", x)
  }
  sort(unique(residues))
}

grin3d_feature_category <- function(type, description = "") {
  type.lower <- tolower(type)
  description.lower <- tolower(description)
  # Description keywords refine generic UniProt Region records only. A domain
  # remains a domain even when its narrative mentions disorder or interaction.
  if (type.lower == "region" && grepl("disorder", description.lower)) return("disorder")
  if (type.lower == "region" && grepl("interaction", description.lower)) return("interaction_region")
  if (type.lower == "active site") return("catalytic_site")
  if (type.lower %in% c("binding site", "site")) return("binding_site")
  if (type.lower %in% c(
    "modified residue", "cross-link", "glycosylation", "lipidation",
    "disulfide bond"
  )) return("ptm")
  if (type.lower == "motif") return("motif")
  if (type.lower %in% c(
    "domain", "dna binding", "homologous_superfamily", "family"
  )) return("domain")
  if (type.lower == "conserved_site") return("conserved_site")
  "functional_region"
}

# Three-level annotation hierarchy:
#   annotation_class -> feature_subcategory -> source-native type/name.
# The original source fields are always retained; these deterministic rules only
# provide consistent grouping across UniProtKB and InterPro.
grin3d_feature_classification <- function(type, name = "", description = "") {
  category <- grin3d_feature_category(type, paste(name, description))
  text <- tolower(paste(type, name, description))
  identity.text <- tolower(paste(type, name))
  type.lower <- tolower(type)

  annotation.class <- switch(
    category,
    domain = "domain",
    binding_site = "site",
    catalytic_site = "site",
    conserved_site = "site",
    ptm = "modification",
    motif = "motif",
    interaction_region = "interaction",
    disorder = "structural_region",
    functional_region = "functional_region",
    "functional_region"
  )

  subcategory <- switch(
    category,
    domain = if (type.lower == "family") {
      "protein_family"
    } else if (type.lower == "homologous_superfamily") {
      "homologous_superfamily"
    } else if (grepl("tetrameri[sz]|oligomeri[sz]", identity.text)) {
      "oligomerization_domain"
    } else if (grepl("transactivation", identity.text)) {
      "transactivation_domain"
    } else if (grepl("zinc finger", identity.text)) {
      "zinc_finger_domain"
    } else if (grepl("dna[- ]?binding|dna binding", identity.text)) {
      "dna_binding_domain"
    } else {
      "protein_domain"
    },
    binding_site = if (grepl("zn\\(|zinc|metal|iron|magnesium|calcium|manganese|copper", text)) {
      "metal_binding_site"
    } else if (grepl("dna|rna|nucleic acid", text)) {
      "nucleic_acid_binding_site"
    } else {
      "ligand_binding_site"
    },
    catalytic_site = "catalytic_site",
    conserved_site = "conserved_site",
    ptm = if (grepl("phospho", text)) {
      "phosphorylation"
    } else if (grepl("acetyl", text)) {
      "acetylation"
    } else if (grepl("ubiquitin", text)) {
      "ubiquitination"
    } else if (grepl("methyl", text)) {
      "methylation"
    } else if (grepl("glycosyl|n-linked|o-linked", text)) {
      "glycosylation"
    } else if (grepl("lipid|myristoyl|palmitoyl|farnesyl", text)) {
      "lipidation"
    } else if (grepl("disulfide", text)) {
      "disulfide_bond"
    } else if (grepl("cross-link", text)) {
      "cross_link"
    } else {
      "other_modification"
    },
    motif = if (grepl("nuclear localization|localization signal", text)) {
      "localization_signal"
    } else if (grepl("degron|degradation", text)) {
      "degradation_motif"
    } else {
      "linear_motif"
    },
    interaction_region = if (grepl("interaction with (dna|rna)|nucleic acid", text)) {
      "nucleic_acid_interaction_region"
    } else {
      "protein_interaction_region"
    },
    disorder = "intrinsically_disordered_region",
    functional_region = if (type.lower == "coiled coil") {
      "coiled_coil"
    } else if (type.lower == "compositional bias") {
      "compositional_bias"
    } else if (type.lower == "repeat") {
      "repeat_region"
    } else if (type.lower == "transmembrane") {
      "transmembrane_region"
    } else if (type.lower == "topological domain") {
      "topological_domain"
    } else {
      "other_functional_region"
    },
    "other_functional_region"
  )

  list(
    annotation_class = annotation.class,
    feature_category = category,
    feature_subcategory = subcategory,
    classification_rule = paste0("grin3d_annotation_hierarchy_v1:", subcategory)
  )
}

grin3d_bind_rows <- function(rows, columns) {
  if (!length(rows)) return(grin3d_empty_table(columns))
  normalized <- lapply(rows, function(row) {
    values <- setNames(as.list(rep("", length(columns))), columns)
    shared <- intersect(names(row), columns)
    values[shared] <- row[shared]
    as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
  })
  out <- do.call(rbind, normalized)
  rownames(out) <- NULL
  out
}

grin3d_fetch_json <- function(url, cache.file, refresh = FALSE, offline = FALSE) {
  if (file.exists(cache.file) && !refresh) {
    return(jsonlite::read_json(cache.file, simplifyVector = FALSE))
  }
  if (offline) stop("Offline mode requires cached response: ", cache.file)
  dir.create(dirname(cache.file), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(cache.file, ".tmp")
  status <- tryCatch(
    utils::download.file(
      url, temporary, mode = "wb", quiet = TRUE,
      headers = c(Accept = "application/json", `User-Agent` = "GRIN3D-R-annotation/1.0")
    ),
    error = function(error) error
  )
  if (inherits(status, "error") || !file.exists(temporary)) {
    stop("Unable to retrieve ", url, ": ", conditionMessage(status))
  }
  if (!file.rename(temporary, cache.file)) stop("Could not cache: ", cache.file)
  jsonlite::read_json(cache.file, simplifyVector = FALSE)
}

grin3d_file_timestamp <- function(file) {
  modified <- file.info(file)$mtime
  if (is.na(modified)) return(format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  format(modified, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

# Convert exon-based CNA candidates into the residue-coordinate input expected
# by the annotation mapper. Each candidate receives the union of all amino-acid
# positions encoded by its member exons; shared splice-boundary residues are
# represented once. This is a coordinate bridge, not a new statistical test.
grin3d_prepare_cna_annotation_input <- function(
    cna.cluster.file,
    exon.map.file,
    output.file,
    mapping.file = file.path(dirname(output.file), "CNA_cluster_exon_residue_mapping.csv")) {
  clusters <- utils::read.csv(
    cna.cluster.file, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = character()
  )
  exons <- utils::read.csv(
    exon.map.file, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = character()
  )
  required.clusters <- c("protein", "cluster_id", "exon_ids")
  required.exons <- c("exon_id", "aa_start", "aa_end")
  missing.clusters <- setdiff(required.clusters, names(clusters))
  missing.exons <- setdiff(required.exons, names(exons))
  if (length(missing.clusters)) {
    stop("Missing CNA cluster columns: ", paste(missing.clusters, collapse = ", "))
  }
  if (length(missing.exons)) {
    stop("Missing exon-map columns: ", paste(missing.exons, collapse = ", "))
  }
  if (anyDuplicated(exons$exon_id)) stop("Exon IDs must be unique in the exon map")
  exons$aa_start <- suppressWarnings(as.integer(exons$aa_start))
  exons$aa_end <- suppressWarnings(as.integer(exons$aa_end))
  if (anyNA(exons$aa_start) || anyNA(exons$aa_end) ||
      any(exons$aa_start < 1L) || any(exons$aa_end < exons$aa_start)) {
    stop("Exon amino-acid intervals must be positive integers with aa_start <= aa_end")
  }
  value <- function(row, column, default = "") {
    if (column %in% names(row)) row[[column]][[1L]] else default
  }
  annotation.rows <- list(); mapping.rows <- list()
  for (i in seq_len(nrow(clusters))) {
    row <- clusters[i, , drop = FALSE]
    exon.ids <- trimws(strsplit(row$exon_ids[[1L]], ";", fixed = TRUE)[[1L]])
    exon.ids <- unique(exon.ids[nzchar(exon.ids)])
    indices <- match(exon.ids, exons$exon_id)
    if (!length(exon.ids) || anyNA(indices)) {
      missing.ids <- exon.ids[is.na(indices)]
      stop(
        "CNA cluster ", row$cluster_id[[1L]], " has missing exon-map IDs: ",
        paste(missing.ids, collapse = ", ")
      )
    }
    selected <- exons[indices, , drop = FALSE]
    residues <- sort(unique(unlist(Map(
      seq.int, selected$aa_start, selected$aa_end
    ), use.names = FALSE)))
    segments <- paste(selected$aa_start, selected$aa_end, sep = "-")
    annotation.rows[[length(annotation.rows) + 1L]] <- data.frame(
      protein = row$protein[[1L]], cluster_id = row$cluster_id[[1L]],
      candidate_key = value(row, "candidate_key"),
      tree_sources = value(row, "tree_sources"),
      hotspot_class = value(row, "hotspot_class"),
      significant_1d = value(row, "significant_1d"),
      significant_3d = value(row, "significant_3d"),
      significant_any = value(row, "significant_any"),
      n_subjects = value(row, "n_subjects"), n_events = value(row, "n_events"),
      n_residues = length(residues), residues = paste(residues, collapse = ";"),
      p_1d_joint = value(row, "p_1d_joint"),
      p_3d_joint = value(row, "p_3d_joint"),
      p_any_joint = value(row, "p_any_joint"),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    mapping.rows[[length(mapping.rows) + 1L]] <- data.frame(
      protein = row$protein[[1L]], transcript = value(row, "transcript"),
      cluster_id = row$cluster_id[[1L]], analysis_type = value(row, "analysis_type"),
      exon_ids = paste(exon.ids, collapse = ";"),
      exon_orders = value(row, "exon_orders"),
      residue_segments = paste(segments, collapse = ";"),
      n_residues = length(residues), residues = paste(residues, collapse = ";"),
      mapping_interpretation = paste(
        "Union of residues encoded by candidate exons; shared splice-boundary",
        "residues deduplicated"
      ),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }
  annotation.table <- do.call(rbind, annotation.rows)
  mapping.table <- do.call(rbind, mapping.rows)
  dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)
  grin3d_write_csv(annotation.table, output.file)
  grin3d_write_csv(mapping.table, mapping.file)
  invisible(list(clusters = annotation.table, mapping = mapping.table))
}

grin3d_read_clusters <- function(cluster.file, protein = NULL) {
  clusters <- utils::read.csv(
    cluster.file, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = character()
  )
  missing <- setdiff(c("cluster_id", "residues"), names(clusters))
  if (length(missing)) stop("Missing cluster columns: ", paste(missing, collapse = ", "))
  if (!nrow(clusters)) stop("Cluster file contains no rows")
  proteins <- unique(clusters$protein[nzchar(clusters$protein %||% "")])
  if (is.null(protein)) {
    if (length(proteins) != 1L) stop("Supply protein when the file has no unique protein label")
    protein <- proteins[[1L]]
  } else if (length(proteins) && !identical(proteins, protein)) {
    stop("Requested protein does not match the cluster file")
  }
  if (any(!nzchar(clusters$cluster_id)) || anyDuplicated(clusters$cluster_id)) {
    stop("Cluster IDs must be nonempty and unique")
  }
  for (column in setdiff(grin3d_cluster_columns, names(clusters))) clusters[[column]] <- ""
  clusters <- clusters[, grin3d_cluster_columns, drop = FALSE]
  parsed <- lapply(clusters$residues, grin3d_parse_residues)
  clusters$residues <- vapply(parsed, paste, collapse = ";", character(1L))
  list(table = clusters, residues = parsed, protein = protein)
}

grin3d_uniprot_reference_index <- function(data) {
  index <- list()
  for (reference in data$references %||% list()) {
    citation <- reference$citation %||% list()
    crossrefs <- citation$citationCrossReferences %||% list()
    ids <- setNames(
      vapply(crossrefs, function(x) grin3d_as_text(x$id), character(1L)),
      vapply(crossrefs, function(x) grin3d_as_text(x$database), character(1L))
    )
    record <- list(
      pubmed_id = if ("PubMed" %in% names(ids)) ids[["PubMed"]] else "",
      doi = if ("DOI" %in% names(ids)) ids[["DOI"]] else "",
      title = grin3d_as_text(citation$title),
      journal = grin3d_as_text(citation$journal),
      year = substr(grin3d_as_text(citation$publicationDate), 1L, 4L),
      authors = grin3d_join(citation$authors %||% list(), sort.values = FALSE)
    )
    for (crossref in crossrefs) {
      key <- paste(grin3d_as_text(crossref$database), grin3d_as_text(crossref$id), sep = ":")
      index[[key]] <- record
    }
    citation.id <- grin3d_as_text(citation$id)
    if (nzchar(citation.id)) index[[citation.id]] <- record
  }
  index
}

grin3d_uniprot_crossrefs <- function(feature) {
  values <- vapply(feature$featureCrossReferences %||% list(), function(item) {
    paste(grin3d_as_text(item$database), grin3d_as_text(item$id), sep = ":")
  }, character(1L))
  ligand.id <- grin3d_as_text((feature$ligand %||% list())$id)
  grin3d_join(c(values, ligand.id))
}

grin3d_extract_uniprot <- function(data, protein, accession, retrieved.at) {
  audit <- data$entryAudit %||% list()
  version <- grin3d_as_text(audit$entryVersion)
  annotation.date <- grin3d_as_text(audit$lastAnnotationUpdateDate)
  reference.index <- grin3d_uniprot_reference_index(data)
  features <- list()
  publications <- list()

  for (raw in data$features %||% list()) {
    type <- grin3d_as_text(raw$type)
    if (!(type %in% grin3d_uniprot_feature_types)) next
    location <- raw$location %||% list()
    start <- as.integer((location$start %||% list())$value)
    end <- as.integer((location$end %||% list())$value)
    if (is.na(start) || is.na(end)) next
    if (end < start) { temporary <- start; start <- end; end <- temporary }
    description <- grin3d_as_text(raw$description)
    ligand <- raw$ligand %||% list()
    name <- description
    if (!nzchar(name)) name <- grin3d_as_text(ligand$name)
    if (!nzchar(name)) name <- type
    evidence <- raw$evidences %||% list()
    evidence.codes <- grin3d_join(lapply(evidence, `[[`, "evidenceCode"))
    evidence.sources <- grin3d_join(lapply(evidence, `[[`, "source"))
    evidence.ids <- grin3d_join(lapply(evidence, `[[`, "id"))
    crossrefs <- grin3d_uniprot_crossrefs(raw)
    classification <- grin3d_feature_classification(type, name, description)
    feature.id <- paste0(
      "UniProtKB:", accession, ":",
      grin3d_stable_id(accession, type, start, end, name, crossrefs)
    )
    features[[length(features) + 1L]] <- list(
      feature_id = feature.id, protein = protein, uniprot_accession = accession,
      source = "UniProtKB", source_record_id = accession,
      annotation_class = classification$annotation_class,
      feature_category = classification$feature_category,
      feature_subcategory = classification$feature_subcategory,
      classification_rule = classification$classification_rule,
      feature_type = type, feature_name = name, description = description,
      start = start, end = end, feature_segments = paste0(start, "-", end),
      start_modifier = grin3d_as_text((location$start %||% list())$modifier),
      end_modifier = grin3d_as_text((location$end %||% list())$modifier),
      evidence_codes = evidence.codes, evidence_sources = evidence.sources,
      evidence_ids = evidence.ids, cross_references = crossrefs, go_terms = "",
      source_url = sprintf(grin3d_annotation_urls$uniprot.web, accession),
      source_record_version = version, source_annotation_date = annotation.date,
      retrieved_at_utc = retrieved.at, segments = list(c(start, end))
    )
    seen <- character()
    for (item in evidence) {
      evidence.source <- grin3d_as_text(item$source)
      evidence.id <- grin3d_as_text(item$id)
      key <- paste(evidence.source, evidence.id, sep = ":")
      if (!nzchar(evidence.source) || !nzchar(evidence.id) || key %in% seen) next
      seen <- c(seen, key)
      citation <- reference.index[[key]] %||% list()
      is.pubmed <- evidence.source == "PubMed"
      publications[[length(publications) + 1L]] <- list(
        feature_id = feature.id, protein = protein, uniprot_accession = accession,
        source = "UniProtKB", source_record_id = accession,
        evidence_scope = "feature_evidence",
        evidence_code = grin3d_as_text(item$evidenceCode),
        publication_id_type = evidence.source, publication_id = evidence.id,
        pubmed_id = citation$pubmed_id %||% if (is.pubmed) evidence.id else "",
        doi = citation$doi %||% "", title = citation$title %||% "",
        journal = citation$journal %||% "", year = citation$year %||% "",
        authors = citation$authors %||% "",
        publication_url = if (is.pubmed) {
          paste0("https://pubmed.ncbi.nlm.nih.gov/", evidence.id, "/")
        } else ""
      )
    }
  }
  source <- list(
    source = "UniProtKB", source_record_id = accession, protein = protein,
    uniprot_accession = accession,
    api_url = sprintf(grin3d_annotation_urls$uniprot.api, accession),
    web_url = sprintf(grin3d_annotation_urls$uniprot.web, accession),
    record_version = version, annotation_date = annotation.date,
    retrieved_at_utc = retrieved.at,
    raw_response_file = paste0("raw/uniprot_", accession, ".json"),
    notes = "Reviewed UniProtKB feature records; feature-level evidence retained."
  )
  list(features = features, publications = publications, sources = list(source))
}

grin3d_interpro_name <- function(metadata) {
  name <- metadata$name %||% ""
  if (is.list(name)) grin3d_as_text(name$name %||% name$short) else grin3d_as_text(name)
}

grin3d_interpro_description <- function(detail) {
  descriptions <- (detail$metadata %||% list())$description %||% list()
  paste(vapply(descriptions, function(item) {
    grin3d_html_text((item %||% list())$text)
  }, character(1L)), collapse = " ")
}

grin3d_interpro_publications <- function(detail, feature) {
  literature <- (detail$metadata %||% list())$literature %||% list()
  rows <- list()
  for (key in names(literature)) {
    item <- literature[[key]] %||% list()
    pubmed <- grin3d_as_text(item$PMID)
    doi.url <- grin3d_as_text(item$DOI_URL)
    rows[[length(rows) + 1L]] <- list(
      feature_id = feature$feature_id, protein = feature$protein,
      uniprot_accession = feature$uniprot_accession, source = "InterPro",
      source_record_id = feature$source_record_id,
      evidence_scope = "entry_literature",
      evidence_code = "InterPro_entry_literature",
      publication_id_type = if (nzchar(pubmed)) "PubMed" else "InterPro",
      publication_id = if (nzchar(pubmed)) pubmed else key,
      pubmed_id = pubmed, doi = grin3d_normalize_doi(doi.url),
      title = grin3d_as_text(item$title),
      journal = grin3d_as_text(item$ISO_journal %||% item$medline_journal),
      year = grin3d_as_text(item$year),
      authors = grin3d_join(item$authors %||% list(), sort.values = FALSE),
      publication_url = if (nzchar(pubmed)) {
        paste0("https://pubmed.ncbi.nlm.nih.gov/", pubmed, "/")
      } else if (nzchar(grin3d_normalize_doi(doi.url))) {
        paste0("https://doi.org/", grin3d_normalize_doi(doi.url))
      } else ""
    )
  }
  rows
}

grin3d_extract_interpro <- function(matches, details, protein, accession, retrieved.at) {
  features <- list(); publications <- list(); sources <- list()
  for (result in matches$results %||% list()) {
    metadata <- result$metadata %||% list()
    entry <- grin3d_as_text(metadata$accession)
    type <- grin3d_as_text(metadata$type)
    name <- grin3d_interpro_name(metadata)
    detail <- details[[entry]] %||% list()
    description <- grin3d_interpro_description(detail)
    go.terms <- grin3d_join(lapply(metadata$go_terms %||% list(), `[[`, "identifier"))
    member.xrefs <- character()
    for (database in names(metadata$member_databases %||% list())) {
      identifiers <- names(metadata$member_databases[[database]] %||% list())
      member.xrefs <- c(member.xrefs, paste(database, identifiers, sep = ":"))
    }
    location.index <- 0L
    for (protein.match in result$proteins %||% list()) {
      for (location in protein.match$entry_protein_locations %||% list()) {
        fragments <- location$fragments %||% list()
        segments <- lapply(fragments, function(fragment) {
          c(as.integer(fragment$start), as.integer(fragment$end))
        })
        segments <- segments[vapply(segments, function(x) length(x) == 2L && !anyNA(x), logical(1L))]
        if (!length(segments)) next
        location.index <- location.index + 1L
        starts <- vapply(segments, `[[`, integer(1L), 1L)
        ends <- vapply(segments, `[[`, integer(1L), 2L)
        segment.text <- paste(paste(starts, ends, sep = "-"), collapse = ";")
        classification <- grin3d_feature_classification(type, name, description)
        feature <- list(
          feature_id = paste0(
            "InterPro:", entry, ":", location.index, ":",
            grin3d_stable_id(segment.text)
          ),
          protein = protein, uniprot_accession = accession, source = "InterPro",
          source_record_id = entry,
          annotation_class = classification$annotation_class,
          feature_category = classification$feature_category,
          feature_subcategory = classification$feature_subcategory,
          classification_rule = classification$classification_rule,
          feature_type = type, feature_name = name, description = description,
          start = min(starts), end = max(ends), feature_segments = segment.text,
          start_modifier = "EXACT", end_modifier = "EXACT",
          evidence_codes = "database_signature_match",
          evidence_sources = "InterPro", evidence_ids = entry,
          cross_references = grin3d_join(member.xrefs), go_terms = go.terms,
          source_url = sprintf(grin3d_annotation_urls$interpro.web, entry),
          source_record_version = "", source_annotation_date = "",
          retrieved_at_utc = retrieved.at, segments = segments
        )
        features[[length(features) + 1L]] <- feature
        publications <- c(publications, grin3d_interpro_publications(detail, feature))
      }
    }
    sources[[length(sources) + 1L]] <- list(
      source = "InterPro", source_record_id = entry, protein = protein,
      uniprot_accession = accession,
      api_url = sprintf(grin3d_annotation_urls$interpro.entry, entry),
      web_url = sprintf(grin3d_annotation_urls$interpro.web, entry),
      record_version = "", annotation_date = "", retrieved_at_utc = retrieved.at,
      raw_response_file = paste0("raw/interpro_entry_", entry, ".json"),
      notes = paste(
        "InterPro match and entry-level literature; literature is not",
        "protein-match-specific."
      )
    )
  }
  list(features = features, publications = publications, sources = sources)
}

grin3d_pdbe_records <- function(data, accession) {
  record <- data[[accession]] %||% data[[toupper(accession)]] %||% list()
  record$data %||% list()
}

grin3d_pdbe_segments <- function(record) {
  segments <- lapply(record$residues %||% list(), function(residue) {
    if (toupper(grin3d_as_text(residue$indexType)) != "UNIPROT") return(NULL)
    start <- suppressWarnings(as.integer(residue$startIndex))
    end <- suppressWarnings(as.integer(residue$endIndex))
    if (is.na(start) || is.na(end)) return(NULL)
    if (end < start) { temporary <- start; start <- end; end <- temporary }
    c(start, end)
  })
  segments <- segments[!vapply(segments, is.null, logical(1L))]
  if (!length(segments)) return(list())
  keys <- vapply(segments, paste, collapse = "-", character(1L))
  segments[!duplicated(keys)]
}

grin3d_pdbe_pdb_ids <- function(record) {
  ids <- unlist(lapply(record$residues %||% list(), function(residue) {
    c(
      unlist(lapply(residue$interactingPDBEntries %||% list(), `[[`, "pdbId")),
      residue$allPDBEntries %||% list()
    )
  }), recursive = TRUE, use.names = FALSE)
  grin3d_join(toupper(ids))
}

grin3d_extract_pdbe_interfaces <- function(data, protein, accession, retrieved.at) {
  features <- list()
  for (record in grin3d_pdbe_records(data, accession)) {
    segments <- grin3d_pdbe_segments(record)
    if (!length(segments)) next
    partner.name <- grin3d_as_text(record$name)
    partner.accession <- grin3d_as_text(record$accession)
    partner.type <- toupper(grin3d_as_text((record$additionalData %||% list())$type))
    identity <- tolower(paste(partner.name, partner.accession))
    subcategory <- if (partner.type == "UNP" && partner.accession == accession) {
      "protein_self_interface"
    } else if (partner.type == "UNP") {
      "protein_protein_interface"
    } else if (grepl("dna|rna|nucleic", identity)) {
      "nucleic_acid_interface"
    } else {
      "macromolecular_interface"
    }
    starts <- vapply(segments, `[[`, integer(1L), 1L)
    ends <- vapply(segments, `[[`, integer(1L), 2L)
    segment.text <- paste(paste(starts, ends, sep = "-"), collapse = ";")
    pdb.ids <- grin3d_pdbe_pdb_ids(record)
    label <- if (nzchar(partner.name)) partner.name else partner.accession
    features[[length(features) + 1L]] <- list(
      feature_id = paste0(
        "PDBe-KB:interface:", grin3d_stable_id(accession, partner.accession, label, segment.text)
      ),
      protein = protein, uniprot_accession = accession, source = "PDBe-KB/SIFTS",
      source_record_id = partner.accession,
      annotation_class = "interaction", feature_category = "interaction_interface",
      feature_subcategory = subcategory,
      classification_rule = paste0("pdbe_kb_sifts_v1:", subcategory),
      feature_type = "experimental interface", feature_name = paste("Interface with", label),
      description = paste(
        "UniProt-indexed contact residues observed in PDB structures; partner type",
        if (nzchar(partner.type)) partner.type else "unspecified"
      ),
      start = min(starts), end = max(ends), feature_segments = segment.text,
      start_modifier = "EXACT", end_modifier = "EXACT",
      evidence_codes = "experimental_structure_contact",
      evidence_sources = "PDBe-KB;PDB;SIFTS", evidence_ids = pdb.ids,
      cross_references = grin3d_join(c(
        if (nzchar(partner.accession)) paste0("partner:", partner.accession) else "",
        if (nzchar(pdb.ids)) paste0("PDB:", unlist(strsplit(pdb.ids, ";", fixed = TRUE))) else ""
      )),
      go_terms = "", source_url = sprintf(grin3d_annotation_urls$pdbe.web, accession),
      source_record_version = "", source_annotation_date = "",
      retrieved_at_utc = retrieved.at, segments = segments,
      relationship_type = "experimental_interface_contact"
    )
  }
  features
}

grin3d_extract_pdbe_ligands <- function(data, protein, accession, retrieved.at) {
  features <- list(); ligand.rows <- list()
  for (record in grin3d_pdbe_records(data, accession)) {
    additional <- record$additionalData %||% list()
    is.solvent <- isTRUE(additional$isSolvent)
    if (is.solvent) next
    segments <- grin3d_pdbe_segments(record)
    if (!length(segments)) next
    ligand.code <- grin3d_as_text(record$accession)
    ligand.name <- grin3d_as_text(record$name)
    chembl.id <- grin3d_as_text(additional$chemblId)
    drugbank.id <- grin3d_as_text(additional$drugBankId)
    cofactor.id <- grin3d_as_text(additional$coFactorId)
    reaction.id <- grin3d_as_text(additional$reactionId)
    subcategory <- if (nzchar(chembl.id) || nzchar(drugbank.id)) {
      "drug_linked_ligand_site"
    } else if (nzchar(cofactor.id)) {
      "cofactor_binding_site"
    } else {
      "small_molecule_binding_site"
    }
    starts <- vapply(segments, `[[`, integer(1L), 1L)
    ends <- vapply(segments, `[[`, integer(1L), 2L)
    segment.text <- paste(paste(starts, ends, sep = "-"), collapse = ";")
    pdb.ids <- grin3d_pdbe_pdb_ids(record)
    label <- trimws(paste(ligand.name, if (nzchar(ligand.code)) paste0("(", ligand.code, ")") else ""))
    features[[length(features) + 1L]] <- list(
      feature_id = paste0(
        "PDBe-KB:ligand:", grin3d_stable_id(accession, ligand.code, label, segment.text)
      ),
      protein = protein, uniprot_accession = accession, source = "PDBe-KB/SIFTS",
      source_record_id = ligand.code,
      annotation_class = "site", feature_category = "binding_site",
      feature_subcategory = subcategory,
      classification_rule = paste0("pdbe_kb_sifts_v1:", subcategory),
      feature_type = "experimental ligand contact", feature_name = paste("Contact with", label),
      description = "UniProt-indexed protein residues contacting a non-solvent ligand in a PDB structure.",
      start = min(starts), end = max(ends), feature_segments = segment.text,
      start_modifier = "EXACT", end_modifier = "EXACT",
      evidence_codes = "experimental_structure_contact",
      evidence_sources = "PDBe-KB;PDB;SIFTS", evidence_ids = pdb.ids,
      cross_references = grin3d_join(c(
        if (nzchar(ligand.code)) paste0("PDB-CCD:", ligand.code) else "",
        if (nzchar(chembl.id)) paste0("ChEMBL:", chembl.id) else "",
        if (nzchar(drugbank.id)) paste0("DrugBank:", drugbank.id) else ""
      )),
      go_terms = "", source_url = sprintf(grin3d_annotation_urls$pdbe.web, accession),
      source_record_version = "", source_annotation_date = "",
      retrieved_at_utc = retrieved.at, segments = segments,
      relationship_type = "experimental_ligand_contact"
    )
    ligand.rows[[length(ligand.rows) + 1L]] <- list(
      protein = protein, uniprot_accession = accession, ligand_code = ligand.code,
      ligand_name = ligand.name,
      contact_residues = grin3d_join(grin3d_feature_positions(features[[length(features)]]), FALSE),
      pdb_ids = pdb.ids, chembl_id = chembl.id, drugbank_id = drugbank.id,
      cofactor_id = cofactor.id, reaction_id = reaction.id, is_solvent = is.solvent,
      evidence_scope = "residue_level_experimental_contact",
      source_url = sprintf(grin3d_annotation_urls$pdbe.web, accession),
      retrieved_at_utc = retrieved.at
    )
  }
  list(features = features, ligands = ligand.rows)
}

grin3d_extract_chembl <- function(target.data, mechanism.data, molecule.data,
                                  protein, accession, retrieved.at) {
  target.rows <- list(); mechanism.rows <- list()
  targets <- target.data$targets %||% list()
  target.index <- list()
  for (target in targets) {
    target.id <- grin3d_as_text(target$target_chembl_id)
    if (!nzchar(target.id)) next
    components <- target$target_components %||% list()
    accessions <- vapply(components, function(x) grin3d_as_text(x$accession), character(1L))
    relationships <- vapply(components, function(x) grin3d_as_text(x$relationship), character(1L))
    target.index[[target.id]] <- target
    target.rows[[length(target.rows) + 1L]] <- list(
      protein = protein, uniprot_accession = accession, target_chembl_id = target.id,
      target_name = grin3d_as_text(target$pref_name), target_type = grin3d_as_text(target$target_type),
      organism = grin3d_as_text(target$organism), tax_id = grin3d_as_text(target$tax_id),
      component_accessions = grin3d_join(accessions),
      component_relationships = grin3d_join(relationships),
      evidence_scope = "protein_or_complex_target_record",
      source_url = sprintf(grin3d_annotation_urls$chembl.target.web, target.id),
      retrieved_at_utc = retrieved.at
    )
  }
  molecule.index <- list()
  for (molecule in molecule.data$molecules %||% list()) {
    molecule.index[[grin3d_as_text(molecule$molecule_chembl_id)]] <- molecule
  }
  for (mechanism in mechanism.data$mechanisms %||% list()) {
    target.id <- grin3d_as_text(mechanism$target_chembl_id)
    molecule.id <- grin3d_as_text(mechanism$molecule_chembl_id)
    target <- target.index[[target.id]] %||% list()
    molecule <- molecule.index[[molecule.id]] %||% list()
    mechanism.rows[[length(mechanism.rows) + 1L]] <- list(
      protein = protein, uniprot_accession = accession, target_chembl_id = target.id,
      target_name = grin3d_as_text(target$pref_name), target_type = grin3d_as_text(target$target_type),
      molecule_chembl_id = molecule.id, molecule_name = grin3d_as_text(molecule$pref_name),
      molecule_type = grin3d_as_text(molecule$molecule_type), max_phase = grin3d_as_text(molecule$max_phase),
      first_approval = grin3d_as_text(molecule$first_approval),
      action_type = grin3d_as_text(mechanism$action_type),
      mechanism_of_action = grin3d_as_text(mechanism$mechanism_of_action),
      mechanism_comment = grin3d_as_text(mechanism$mechanism_comment),
      direct_interaction = grin3d_as_text(mechanism$direct_interaction),
      molecular_mechanism = grin3d_as_text(mechanism$molecular_mechanism),
      evidence_scope = "target_level_drug_mechanism_not_residue_evidence",
      target_url = sprintf(grin3d_annotation_urls$chembl.target.web, target.id),
      molecule_url = sprintf(grin3d_annotation_urls$chembl.molecule.web, molecule.id),
      retrieved_at_utc = retrieved.at
    )
  }
  list(targets = target.rows, mechanisms = mechanism.rows)
}

grin3d_feature_positions <- function(feature) {
  sort(unique(unlist(lapply(feature$segments, function(segment) {
    seq.int(segment[[1L]], segment[[2L]])
  }))))
}

grin3d_map_annotations <- function(clusters, cluster.residues, features, accession) {
  overlap.rows <- list(); residue.rows <- list(); summary.rows <- list()
  feature.positions <- lapply(features, grin3d_feature_positions)
  for (i in seq_len(nrow(clusters))) {
    cluster <- clusters[i, , drop = FALSE]
    residues <- cluster.residues[[i]]
    matched.features <- list(); annotated <- integer()
    for (j in seq_along(features)) {
      feature <- features[[j]]
      positions <- feature.positions[[j]]
      shared <- intersect(residues, positions)
      if (!length(shared)) next
      relationship <- feature$relationship_type %||% if (length(positions) == 1L) {
        "exact_residue_overlap"
      } else "interval_overlap"
      base <- as.list(cluster[1L, grin3d_cluster_columns, drop = FALSE])
      overlap.rows[[length(overlap.rows) + 1L]] <- c(base, list(
        feature_id = feature$feature_id, annotation_source = feature$source,
        source_record_id = feature$source_record_id,
        annotation_class = feature$annotation_class,
        feature_category = feature$feature_category,
        feature_subcategory = feature$feature_subcategory,
        classification_rule = feature$classification_rule,
        feature_type = feature$feature_type, feature_name = feature$feature_name,
        feature_start = feature$start, feature_end = feature$end,
        feature_segments = feature$feature_segments, relationship = relationship,
        overlap_residues = paste(shared, collapse = ";"),
        n_overlap_residues = length(shared),
        fraction_cluster_residues_overlapping = round(length(shared) / length(residues), 6L),
        fraction_feature_positions_observed_in_cluster = round(length(shared) / length(positions), 6L),
        evidence_codes = feature$evidence_codes, evidence_ids = feature$evidence_ids,
        source_url = feature$source_url
      ))
      for (residue in shared) {
        residue.rows[[length(residue.rows) + 1L]] <- list(
          protein = cluster$protein, uniprot_accession = accession,
          cluster_id = cluster$cluster_id, residue = residue,
          feature_id = feature$feature_id, annotation_source = feature$source,
          source_record_id = feature$source_record_id,
          annotation_class = feature$annotation_class,
          feature_category = feature$feature_category,
          feature_subcategory = feature$feature_subcategory,
          classification_rule = feature$classification_rule,
          feature_type = feature$feature_type, feature_name = feature$feature_name,
          relationship = relationship, evidence_codes = feature$evidence_codes,
          evidence_ids = feature$evidence_ids, source_url = feature$source_url
        )
      }
      matched.features[[length(matched.features) + 1L]] <- feature
      annotated <- union(annotated, shared)
    }
    names.by.category <- function(categories) {
      grin3d_join(vapply(
        matched.features[vapply(matched.features, function(x) x$feature_category %in% categories, logical(1L))],
        `[[`, character(1L), "feature_name"
      ))
    }
    base <- as.list(cluster[1L, grin3d_cluster_columns, drop = FALSE])
    summary.rows[[length(summary.rows) + 1L]] <- c(base, list(
      n_annotation_features = length(matched.features),
      n_annotated_cluster_residues = length(annotated),
      annotated_cluster_residues = paste(sort(annotated), collapse = ";"),
      n_unannotated_cluster_residues = length(setdiff(residues, annotated)),
      unannotated_cluster_residues = paste(sort(setdiff(residues, annotated)), collapse = ";"),
      annotation_sources = grin3d_join(vapply(matched.features, `[[`, character(1L), "source")),
      annotation_classes = grin3d_join(vapply(matched.features, `[[`, character(1L), "annotation_class")),
      feature_categories = grin3d_join(vapply(matched.features, `[[`, character(1L), "feature_category")),
      feature_subcategories = grin3d_join(vapply(matched.features, `[[`, character(1L), "feature_subcategory")),
      domain_features = names.by.category("domain"), motif_features = names.by.category("motif"),
      site_features = names.by.category(c("binding_site", "catalytic_site", "conserved_site")),
      ptm_features = names.by.category("ptm"),
      interaction_features = names.by.category(c("interaction_region", "interaction_interface")),
      disorder_features = names.by.category("disorder")
    ))
  }
  list(overlaps = overlap.rows, residues = residue.rows, summaries = summary.rows)
}

grin3d_write_csv <- function(table, file) {
  utils::write.csv(table, file, row.names = FALSE, na = "")
}

#' Annotate GRIN3D mutation clusters with four independently sourced adapters.
#'
#' @return Invisibly returns all output tables and protein metadata.
annotate_grin3d_clusters <- function(
    cluster.file,
    uniprot.accession,
    output.dir,
    protein = NULL,
    sources = c("uniprot", "interpro", "pdbe", "chembl"),
    refresh = FALSE,
    offline = FALSE) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Install the R package 'jsonlite' before running annotations")
  }
  sources <- unique(tolower(sources))
  unknown <- setdiff(sources, c("uniprot", "interpro", "pdbe", "chembl"))
  if (length(unknown)) stop("Unknown annotation source(s): ", paste(unknown, collapse = ", "))
  parsed <- grin3d_read_clusters(cluster.file, protein)
  protein <- parsed$protein
  accession <- toupper(trimws(uniprot.accession))
  raw.dir <- file.path(output.dir, "raw")
  dir.create(raw.dir, recursive = TRUE, showWarnings = FALSE)
  features <- list(); publications <- list(); source.rows <- list()
  pdbe.ligand.rows <- list(); chembl.target.rows <- list(); chembl.mechanism.rows <- list()
  protein.length <- NA_integer_

  if ("uniprot" %in% sources) {
    uniprot.file <- file.path(raw.dir, paste0("uniprot_", accession, ".json"))
    data <- grin3d_fetch_json(
      sprintf(grin3d_annotation_urls$uniprot.api, accession),
      uniprot.file, refresh, offline
    )
    protein.length <- as.integer((data$sequence %||% list())$length)
    extracted <- grin3d_extract_uniprot(
      data, protein, accession, grin3d_file_timestamp(uniprot.file)
    )
    features <- c(features, extracted$features)
    publications <- c(publications, extracted$publications)
    source.rows <- c(source.rows, extracted$sources)
  }
  invalid <- sort(unique(unlist(parsed$residues)))[sort(unique(unlist(parsed$residues))) > protein.length]
  if (!is.na(protein.length) && length(invalid)) {
    stop("Cluster residues exceed UniProt length ", protein.length, ": ", paste(invalid, collapse = ", "))
  }

  if ("interpro" %in% sources) {
    match.file <- file.path(raw.dir, paste0("interpro_matches_", accession, ".json"))
    matches <- grin3d_fetch_json(
      sprintf(grin3d_annotation_urls$interpro.matches, accession),
      match.file, refresh, offline
    )
    match.retrieved.at <- grin3d_file_timestamp(match.file)
    entries <- sort(unique(vapply(matches$results %||% list(), function(result) {
      grin3d_as_text((result$metadata %||% list())$accession)
    }, character(1L))))
    entries <- entries[nzchar(entries)]
    details <- setNames(vector("list", length(entries)), entries)
    for (entry in entries) {
      details[[entry]] <- grin3d_fetch_json(
        sprintf(grin3d_annotation_urls$interpro.entry, entry),
        file.path(raw.dir, paste0("interpro_entry_", entry, ".json")), refresh, offline
      )
    }
    extracted <- grin3d_extract_interpro(
      matches, details, protein, accession, match.retrieved.at
    )
    features <- c(features, extracted$features)
    publications <- c(publications, extracted$publications)
    match.source <- list(
      source = "InterPro", source_record_id = accession, protein = protein,
      uniprot_accession = accession,
      api_url = sprintf(grin3d_annotation_urls$interpro.matches, accession),
      web_url = paste0("https://www.ebi.ac.uk/interpro/protein/UniProt/", accession, "/"),
      record_version = "", annotation_date = "", retrieved_at_utc = match.retrieved.at,
      raw_response_file = paste0("raw/interpro_matches_", accession, ".json"),
      notes = "Protein-to-InterPro match response."
    )
    source.rows <- c(source.rows, list(match.source), extracted$sources)
  }

  if ("pdbe" %in% sources) {
    interface.file <- file.path(raw.dir, paste0("pdbe_interfaces_", accession, ".json"))
    ligand.file <- file.path(raw.dir, paste0("pdbe_ligands_", accession, ".json"))
    interface.data <- grin3d_fetch_json(
      sprintf(grin3d_annotation_urls$pdbe.interfaces, accession),
      interface.file, refresh, offline
    )
    ligand.data <- grin3d_fetch_json(
      sprintf(grin3d_annotation_urls$pdbe.ligands, accession),
      ligand.file, refresh, offline
    )
    interfaces <- grin3d_extract_pdbe_interfaces(
      interface.data, protein, accession, grin3d_file_timestamp(interface.file)
    )
    ligands <- grin3d_extract_pdbe_ligands(
      ligand.data, protein, accession, grin3d_file_timestamp(ligand.file)
    )
    features <- c(features, interfaces, ligands$features)
    pdbe.ligand.rows <- c(pdbe.ligand.rows, ligands$ligands)
    source.rows <- c(source.rows, list(
      list(
        source = "PDBe-KB/SIFTS", source_record_id = paste0(accession, ":interfaces"),
        protein = protein, uniprot_accession = accession,
        api_url = sprintf(grin3d_annotation_urls$pdbe.interfaces, accession),
        web_url = sprintf(grin3d_annotation_urls$pdbe.web, accession),
        record_version = "", annotation_date = "",
        retrieved_at_utc = grin3d_file_timestamp(interface.file),
        raw_response_file = paste0("raw/pdbe_interfaces_", accession, ".json"),
        notes = paste(
          "Experimental interface residues mapped to UniProt coordinates by SIFTS;",
          "PDB IDs are retained as evidence."
        )
      ),
      list(
        source = "PDBe-KB/SIFTS", source_record_id = paste0(accession, ":ligands"),
        protein = protein, uniprot_accession = accession,
        api_url = sprintf(grin3d_annotation_urls$pdbe.ligands, accession),
        web_url = sprintf(grin3d_annotation_urls$pdbe.web, accession),
        record_version = "", annotation_date = "",
        retrieved_at_utc = grin3d_file_timestamp(ligand.file),
        raw_response_file = paste0("raw/pdbe_ligands_", accession, ".json"),
        notes = paste(
          "Experimental non-solvent ligand-contact residues mapped to UniProt",
          "coordinates by SIFTS; a drug cross-reference is retained but not inferred."
        )
      )
    ))
  }

  if ("chembl" %in% sources) {
    target.file <- file.path(raw.dir, paste0("chembl_targets_", accession, ".json"))
    target.url <- sprintf(grin3d_annotation_urls$chembl.targets, accession)
    target.data <- grin3d_fetch_json(target.url, target.file, refresh, offline)
    target.ids <- sort(unique(vapply(target.data$targets %||% list(), function(target) {
      grin3d_as_text(target$target_chembl_id)
    }, character(1L))))
    target.ids <- target.ids[nzchar(target.ids)]
    mechanism.data <- list(mechanisms = list())
    molecule.data <- list(molecules = list())
    mechanism.url <- ""
    molecule.url <- ""
    if (length(target.ids)) {
      mechanism.url <- sprintf(
        grin3d_annotation_urls$chembl.mechanisms, paste(target.ids, collapse = ",")
      )
      mechanism.file <- file.path(raw.dir, paste0("chembl_mechanisms_", accession, ".json"))
      mechanism.data <- grin3d_fetch_json(mechanism.url, mechanism.file, refresh, offline)
      molecule.ids <- sort(unique(vapply(
        mechanism.data$mechanisms %||% list(),
        function(mechanism) grin3d_as_text(mechanism$molecule_chembl_id),
        character(1L)
      )))
      molecule.ids <- molecule.ids[nzchar(molecule.ids)]
      if (length(molecule.ids)) {
        molecule.url <- sprintf(
          grin3d_annotation_urls$chembl.molecules, paste(molecule.ids, collapse = ",")
        )
        molecule.file <- file.path(raw.dir, paste0("chembl_molecules_", accession, ".json"))
        molecule.data <- grin3d_fetch_json(molecule.url, molecule.file, refresh, offline)
      }
    }
    chembl <- grin3d_extract_chembl(
      target.data, mechanism.data, molecule.data, protein, accession,
      grin3d_file_timestamp(target.file)
    )
    chembl.target.rows <- c(chembl.target.rows, chembl$targets)
    chembl.mechanism.rows <- c(chembl.mechanism.rows, chembl$mechanisms)
    source.rows <- c(source.rows, list(list(
      source = "ChEMBL", source_record_id = accession, protein = protein,
      uniprot_accession = accession, api_url = target.url,
      web_url = if (length(target.ids)) {
        sprintf(grin3d_annotation_urls$chembl.target.web, target.ids[[1L]])
      } else "https://www.ebi.ac.uk/chembl/",
      record_version = "", annotation_date = "",
      retrieved_at_utc = grin3d_file_timestamp(target.file),
      raw_response_file = paste0("raw/chembl_targets_", accession, ".json"),
      notes = paste(
        "ChEMBL target and mechanism records are protein/complex-level context,",
        "not residue-level binding evidence. Related mechanism and molecule responses",
        "are cached separately."
      )
    )))
    if (nzchar(mechanism.url)) {
      source.rows <- c(source.rows, list(list(
        source = "ChEMBL", source_record_id = paste0(accession, ":mechanisms"),
        protein = protein, uniprot_accession = accession, api_url = mechanism.url,
        web_url = "https://www.ebi.ac.uk/chembl/", record_version = "",
        annotation_date = "", retrieved_at_utc = grin3d_file_timestamp(mechanism.file),
        raw_response_file = paste0("raw/chembl_mechanisms_", accession, ".json"),
        notes = "Drug-mechanism records for ChEMBL targets containing this accession."
      )))
    }
    if (nzchar(molecule.url)) {
      source.rows <- c(source.rows, list(list(
        source = "ChEMBL", source_record_id = paste0(accession, ":molecules"),
        protein = protein, uniprot_accession = accession, api_url = molecule.url,
        web_url = "https://www.ebi.ac.uk/chembl/", record_version = "",
        annotation_date = "", retrieved_at_utc = grin3d_file_timestamp(molecule.file),
        raw_response_file = paste0("raw/chembl_molecules_", accession, ".json"),
        notes = "Molecule names and development phases for returned drug mechanisms."
      )))
    }
  }

  order.index <- order(
    vapply(features, `[[`, integer(1L), "start"),
    vapply(features, `[[`, integer(1L), "end"),
    vapply(features, `[[`, character(1L), "source")
  )
  features <- features[order.index]
  mapped <- grin3d_map_annotations(parsed$table, parsed$residues, features, accession)

  feature.table <- grin3d_bind_rows(features, grin3d_feature_columns)
  hierarchy.keys <- unique(feature.table[, c(
    "annotation_class", "feature_category", "feature_subcategory",
    "classification_rule"
  ), drop = FALSE])
  hierarchy.table <- do.call(rbind, lapply(seq_len(nrow(hierarchy.keys)), function(i) {
    key <- hierarchy.keys[i, , drop = FALSE]
    keep <- feature.table$annotation_class == key$annotation_class &
      feature.table$feature_category == key$feature_category &
      feature.table$feature_subcategory == key$feature_subcategory &
      feature.table$classification_rule == key$classification_rule
    data.frame(
      key,
      n_features = sum(keep),
      sources = grin3d_join(feature.table$source[keep]),
      source_feature_types = grin3d_join(feature.table$feature_type[keep]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
  hierarchy.table <- hierarchy.table[order(
    hierarchy.table$annotation_class, hierarchy.table$feature_subcategory
  ), , drop = FALSE]
  publication.table <- unique(grin3d_bind_rows(publications, grin3d_publication_columns))
  source.table <- unique(grin3d_bind_rows(source.rows, grin3d_source_columns))
  pdbe.ligand.table <- unique(grin3d_bind_rows(
    pdbe.ligand.rows, grin3d_pdbe_ligand_columns
  ))
  chembl.target.table <- unique(grin3d_bind_rows(
    chembl.target.rows, grin3d_chembl_target_columns
  ))
  chembl.mechanism.table <- unique(grin3d_bind_rows(
    chembl.mechanism.rows, grin3d_chembl_mechanism_columns
  ))
  overlap.columns <- c(
    grin3d_cluster_columns, "feature_id", "annotation_source", "source_record_id",
    "annotation_class", "feature_category", "feature_subcategory",
    "classification_rule", "feature_type", "feature_name", "feature_start",
    "feature_end", "feature_segments", "relationship", "overlap_residues",
    "n_overlap_residues", "fraction_cluster_residues_overlapping",
    "fraction_feature_positions_observed_in_cluster", "evidence_codes",
    "evidence_ids", "source_url"
  )
  residue.columns <- c(
    "protein", "uniprot_accession", "cluster_id", "residue", "feature_id",
    "annotation_source", "source_record_id", "annotation_class",
    "feature_category", "feature_subcategory", "classification_rule", "feature_type",
    "feature_name", "relationship", "evidence_codes", "evidence_ids", "source_url"
  )
  summary.columns <- c(
    grin3d_cluster_columns, "n_annotation_features", "n_annotated_cluster_residues",
    "annotated_cluster_residues", "n_unannotated_cluster_residues",
    "unannotated_cluster_residues", "annotation_sources", "feature_categories",
    "annotation_classes", "feature_subcategories",
    "domain_features", "motif_features", "site_features", "ptm_features",
    "interaction_features", "disorder_features"
  )
  overlap.table <- grin3d_bind_rows(mapped$overlaps, overlap.columns)
  residue.table <- grin3d_bind_rows(mapped$residues, residue.columns)
  summary.table <- grin3d_bind_rows(mapped$summaries, summary.columns)

  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
  grin3d_write_csv(feature.table, file.path(output.dir, "protein_features.csv"))
  grin3d_write_csv(
    hierarchy.table,
    file.path(output.dir, "annotation_classification_summary.csv")
  )
  grin3d_write_csv(publication.table, file.path(output.dir, "feature_publications.csv"))
  grin3d_write_csv(overlap.table, file.path(output.dir, "cluster_feature_overlaps.csv"))
  grin3d_write_csv(residue.table, file.path(output.dir, "cluster_residue_annotations.csv"))
  grin3d_write_csv(summary.table, file.path(output.dir, "cluster_annotation_summary.csv"))
  grin3d_write_csv(source.table, file.path(output.dir, "annotation_sources.csv"))
  grin3d_write_csv(pdbe.ligand.table, file.path(output.dir, "pdbe_ligand_evidence.csv"))
  grin3d_write_csv(chembl.target.table, file.path(output.dir, "chembl_target_evidence.csv"))
  grin3d_write_csv(
    chembl.mechanism.table,
    file.path(output.dir, "chembl_drug_mechanisms.csv")
  )

  category.counts <- table(feature.table$feature_category)
  message(
    "Annotated ", nrow(parsed$table), " clusters for ", protein, " (", accession,
    ") with ", nrow(feature.table), " normalized features and ",
    nrow(publication.table), " feature-linked publication records."
  )
  message(
    "Feature categories: ",
    paste(names(category.counts), as.integer(category.counts), sep = "=", collapse = ", ")
  )
  result <- list(
    protein = protein, accession = accession, protein_length = protein.length,
    clusters = parsed$table, features = feature.table,
    classification_summary = hierarchy.table,
    publications = publication.table, sources = source.table,
    pdbe_ligands = pdbe.ligand.table, chembl_targets = chembl.target.table,
    chembl_mechanisms = chembl.mechanism.table,
    overlaps = overlap.table, residue_annotations = residue.table,
    summaries = summary.table, output_dir = output.dir
  )
  invisible(result)
}

grin3d_category_colors <- c(
  domain = "#4477AA", conserved_site = "#228833", motif = "#AA3377",
  binding_site = "#EE7733", catalytic_site = "#CC3311", ptm = "#EE3377",
  interaction_region = "#66CCEE", interaction_interface = "#0077BB",
  disorder = "#BBBBBB",
  functional_region = "#999933"
)

grin3d_plot_cluster <- function(cluster, overlaps, protein.length, plot.file) {
  categories <- unique(overlaps$feature_category)
  colors <- grin3d_category_colors[categories]
  colors[is.na(colors)] <- "#777777"
  height <- max(5, 2.5 + 0.31 * nrow(overlaps))
  track.labels <- vapply(seq_len(nrow(overlaps)), function(i) {
    row <- overlaps[i, , drop = FALSE]
    label <- paste0(
      "[", row$feature_category, "] ", row$annotation_source, ": ",
      row$feature_name, " [", row$feature_segments, "]"
    )
    if (nchar(label) > 72L) paste0(substr(label, 1L, 69L), "...") else label
  }, character(1L))
  extension <- tolower(tools::file_ext(plot.file))
  if (extension == "svg") {
    grDevices::svg(plot.file, width = 14, height = height, pointsize = 10)
  } else if (extension == "png") {
    grDevices::png(
      plot.file, width = 14, height = height, units = "in", res = 150,
      pointsize = 10
    )
  } else {
    stop("Cluster plot must use an .svg or .png filename")
  }
  old.par <- graphics::par(mar = c(3.5, 29, 3, 1), xpd = NA)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  n <- nrow(overlaps)
  graphics::plot(
    NA, xlim = c(1, protein.length), ylim = c(0, n + 2.4), axes = FALSE,
    xlab = "UniProt residue position", ylab = "",
    main = paste(cluster$protein, cluster$cluster_id, "annotation map")
  )
  graphics::axis(1)
  graphics::abline(v = pretty(c(1, protein.length)), col = "#EEEEEE", lty = 3)
  graphics::segments(1, n + 1.45, protein.length, n + 1.45, lwd = 3, col = "#444444")
  residues <- grin3d_parse_residues(cluster$residues)
  graphics::segments(residues, n + 1.45, residues, n + 1.9, col = "#7A0019", lwd = 1.3)
  graphics::points(residues, rep(n + 1.95, length(residues)), pch = 21, bg = "#7A0019", col = "white", cex = 0.85)
  graphics::mtext("cluster residues", side = 2, at = n + 1.7, las = 1, line = 1, adj = 1, col = "#7A0019")
  for (i in seq_len(n)) {
    row <- overlaps[n - i + 1L, , drop = FALSE]
    y <- i
    color <- grin3d_category_colors[[row$feature_category]] %||% "#777777"
    segment.parts <- strsplit(row$feature_segments, ";", fixed = TRUE)[[1L]]
    segment.bounds <- lapply(segment.parts, function(part) {
      suppressWarnings(as.integer(strsplit(part, "-", fixed = TRUE)[[1L]]))
    })
    segment.bounds <- segment.bounds[vapply(
      segment.bounds, function(bounds) length(bounds) == 2L && !anyNA(bounds), logical(1L)
    )]
    for (bounds in segment.bounds) {
      if (bounds[[1L]] == bounds[[2L]]) {
        graphics::points(bounds[[1L]], y, pch = 23, bg = color, col = color, cex = 1)
      } else {
        graphics::rect(bounds[[1L]], y - 0.25, bounds[[2L]], y + 0.25, col = color, border = NA)
      }
    }
    shared <- grin3d_parse_residues(row$overlap_residues)
    graphics::points(shared, rep(y, length(shared)), pch = 21, bg = "#7A0019", col = "white", cex = 0.75)
    label <- track.labels[[n - i + 1L]]
    graphics::mtext(label, side = 2, at = y, las = 1, line = 1, adj = 1, cex = 0.72)
  }
  graphics::box(bty = "l")
}

grin3d_publication_links <- function(feature.id, publications, max.links = 5L) {
  rows <- publications[publications$feature_id == feature.id & nzchar(publications$publication_url), , drop = FALSE]
  if (!nrow(rows)) return("")
  rows <- rows[!duplicated(rows$publication_url), , drop = FALSE]
  shown <- utils::head(rows, max.links)
  links <- vapply(seq_len(nrow(shown)), function(i) {
    label <- if (nzchar(shown$pubmed_id[[i]])) paste0("PMID ", shown$pubmed_id[[i]]) else shown$publication_id[[i]]
    paste0('<a href="', grin3d_html_escape(shown$publication_url[[i]]), '">', grin3d_html_escape(label), "</a>")
  }, character(1L))
  remainder <- nrow(rows) - nrow(shown)
  paste0(paste(links, collapse = ", "), if (remainder > 0L) paste0(" (+", remainder, " more)") else "")
}

#' Generate a browsable index, HTML detail page, and SVG map for every cluster.
write_cluster_annotation_reports <- function(annotation.result, report.dir = NULL) {
  if (is.null(report.dir)) report.dir <- file.path(annotation.result$output_dir, "cluster_reports")
  dir.create(report.dir, recursive = TRUE, showWarnings = FALSE)
  summaries <- annotation.result$summaries
  overlaps <- annotation.result$overlaps
  publications <- annotation.result$publications
  css <- paste(
    "body{font-family:system-ui,sans-serif;max-width:1500px;margin:2rem auto;padding:0 1rem;color:#222}",
    "table{border-collapse:collapse;width:100%;font-size:.9rem}th,td{border:1px solid #ddd;padding:.45rem;vertical-align:top}",
    "th{background:#f2f4f7;text-align:left;position:sticky;top:0}tr:nth-child(even){background:#fafafa}",
    ".pill{display:inline-block;background:#e9eef5;border-radius:1rem;padding:.15rem .55rem;margin:.1rem}",
    ".note{background:#fff8df;border-left:4px solid #e0a100;padding:.8rem}img{width:100%;min-width:900px}",
    "a{color:#075985}code{background:#f3f4f6;padding:.1rem .25rem}", sep = ""
  )
  for (i in seq_len(nrow(summaries))) {
    cluster <- summaries[i, , drop = FALSE]
    cluster.overlaps <- overlaps[overlaps$cluster_id == cluster$cluster_id, , drop = FALSE]
    cluster.overlaps <- cluster.overlaps[order(cluster.overlaps$feature_category, cluster.overlaps$feature_start), , drop = FALSE]
    stem <- grin3d_safe_name(cluster$cluster_id)
    svg.name <- paste0(stem, ".svg")
    html.name <- paste0(stem, ".html")
    grin3d_plot_cluster(cluster, cluster.overlaps, annotation.result$protein_length, file.path(report.dir, svg.name))
    rows <- vapply(seq_len(nrow(cluster.overlaps)), function(j) {
      row <- cluster.overlaps[j, , drop = FALSE]
      source.link <- paste0('<a href="', grin3d_html_escape(row$source_url), '">', grin3d_html_escape(row$annotation_source), "</a>")
      literature <- grin3d_publication_links(row$feature_id, publications)
      paste0(
        "<tr><td>", grin3d_html_escape(row$feature_category), "</td><td>",
        grin3d_html_escape(row$feature_name), "</td><td>", grin3d_html_escape(row$feature_segments),
        "</td><td>", grin3d_html_escape(row$overlap_residues), "</td><td>",
        grin3d_html_escape(row$relationship), "</td><td>", source.link,
        "</td><td>", grin3d_html_escape(row$evidence_codes), "</td><td>", literature, "</td></tr>"
      )
    }, character(1L))
    page <- c(
      "<!doctype html><html><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width\">",
      paste0("<title>", grin3d_html_escape(cluster$protein), " ", grin3d_html_escape(cluster$cluster_id), " annotations</title><style>", css, "</style></head><body>"),
      '<p><a href="index.html">&#8592; All clusters</a></p>',
      paste0("<h1>", grin3d_html_escape(cluster$protein), " ", grin3d_html_escape(cluster$cluster_id), "</h1>"),
      paste0('<p><span class="pill">', grin3d_html_escape(cluster$hotspot_class), '</span> ',
             '<span class="pill">', cluster$n_residues, " residues</span> ",
             '<span class="pill">', cluster$n_annotation_features, " overlapping features</span></p>"),
      paste0("<p><b>Cluster residues:</b> <code>", grin3d_html_escape(cluster$residues), "</code></p>"),
      '<p class="note"><b>Interpretation boundary:</b> these are coordinate overlaps added after statistical testing. They do not alter hotspot significance. InterPro entry literature describes the feature and is not direct evidence for this cluster.</p>',
      paste0('<div style="overflow-x:auto"><img src="', svg.name, '" alt="Residue and annotation map"></div>'),
      "<h2>Overlapping annotations and evidence</h2>",
      "<table><thead><tr><th>Category</th><th>Feature</th><th>Feature coordinates</th><th>Cluster residues overlapping</th><th>Relationship</th><th>Source</th><th>Evidence code</th><th>Literature</th></tr></thead><tbody>",
      rows, "</tbody></table></body></html>"
    )
    writeLines(page, file.path(report.dir, html.name), useBytes = TRUE)
  }
  index.rows <- vapply(seq_len(nrow(summaries)), function(i) {
    row <- summaries[i, , drop = FALSE]
    page <- paste0(grin3d_safe_name(row$cluster_id), ".html")
    paste0(
      '<tr><td><a href="', page, '">', grin3d_html_escape(row$cluster_id), "</a></td><td>",
      grin3d_html_escape(row$hotspot_class), "</td><td>", row$n_residues, "</td><td>",
      row$n_annotation_features, "</td><td>", row$n_annotated_cluster_residues,
      "</td><td>", grin3d_html_escape(row$feature_categories), "</td><td>",
      grin3d_html_escape(row$site_features), "</td></tr>"
    )
  }, character(1L))
  index <- c(
    "<!doctype html><html><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width\">",
    paste0("<title>", grin3d_html_escape(annotation.result$protein), " cluster annotations</title><style>", css, "</style></head><body>"),
    paste0("<h1>", grin3d_html_escape(annotation.result$protein), " biological annotation browser</h1>"),
    '<p>Select a cluster to view its residue map, overlapping regional annotations, evidence codes, source records, and literature links.</p>',
    '<p class="note">Annotations are descriptive post-processing and remain separate from statistical significance.</p>',
    "<table><thead><tr><th>Cluster</th><th>Hotspot class</th><th>Residues</th><th>Features</th><th>Annotated residues</th><th>Categories</th><th>Sites</th></tr></thead><tbody>",
    index.rows, "</tbody></table></body></html>"
  )
  writeLines(index, file.path(report.dir, "index.html"), useBytes = TRUE)
  message("Cluster annotation browser written to: ", file.path(report.dir, "index.html"))
  invisible(report.dir)
}

#' Write one self-contained interactive HTML dashboard for all clusters.
#'
#' The dashboard is generated entirely by R and needs no Shiny server. All data,
#' plots, source links, and literature links are embedded in one portable file.
write_cluster_annotation_dashboard <- function(
    annotation.result,
    output.file = file.path(
      annotation.result$output_dir, "cluster_annotation_dashboard.html"
    )) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Install the R package 'jsonlite' before writing the dashboard")
  }
  summaries <- annotation.result$summaries
  overlaps <- annotation.result$overlaps
  publications <- annotation.result$publications

  scalar <- function(row, column, default = "") {
    value <- row[[column]]
    if (is.null(value) || !length(value) || is.na(value[[1L]])) default else value[[1L]]
  }
  cluster.records <- lapply(seq_len(nrow(summaries)), function(i) {
    summary <- summaries[i, , drop = FALSE]
    rows <- overlaps[overlaps$cluster_id == summary$cluster_id, , drop = FALSE]
    rows <- rows[order(
      rows$annotation_class, rows$feature_subcategory, rows$feature_start
    ), , drop = FALSE]
    annotations <- lapply(seq_len(nrow(rows)), function(j) {
      row <- rows[j, , drop = FALSE]
      papers <- publications[
        publications$feature_id == row$feature_id &
          nzchar(publications$publication_url), , drop = FALSE
      ]
      if (nrow(papers)) papers <- papers[!duplicated(papers$publication_url), , drop = FALSE]
      literature <- lapply(seq_len(nrow(papers)), function(k) list(
        label = if (nzchar(papers$pubmed_id[[k]])) {
          paste0("PMID ", papers$pubmed_id[[k]])
        } else papers$publication_id[[k]],
        title = papers$title[[k]],
        scope = papers$evidence_scope[[k]],
        url = papers$publication_url[[k]]
      ))
      list(
        feature_id = scalar(row, "feature_id"),
        annotation_class = scalar(row, "annotation_class"),
        feature_category = scalar(row, "feature_category"),
        feature_subcategory = scalar(row, "feature_subcategory"),
        feature_type = scalar(row, "feature_type"),
        feature_name = scalar(row, "feature_name"),
        start = as.integer(scalar(row, "feature_start", 0L)),
        end = as.integer(scalar(row, "feature_end", 0L)),
        segments = scalar(row, "feature_segments"),
        relationship = scalar(row, "relationship"),
        overlap_residues = scalar(row, "overlap_residues"),
        source = scalar(row, "annotation_source"),
        source_record_id = scalar(row, "source_record_id"),
        source_url = scalar(row, "source_url"),
        evidence_codes = scalar(row, "evidence_codes"),
        evidence_ids = scalar(row, "evidence_ids"),
        classification_rule = scalar(row, "classification_rule"),
        literature = literature
      )
    })
    list(
      summary = list(
        protein = scalar(summary, "protein"),
        cluster_id = scalar(summary, "cluster_id"),
        hotspot_class = scalar(summary, "hotspot_class"),
        tree_sources = scalar(summary, "tree_sources"),
        significant_1d = as.logical(scalar(summary, "significant_1d", FALSE)),
        significant_3d = as.logical(scalar(summary, "significant_3d", FALSE)),
        significant_any = as.logical(scalar(summary, "significant_any", FALSE)),
        n_subjects = as.integer(scalar(summary, "n_subjects", 0L)),
        n_events = as.integer(scalar(summary, "n_events", 0L)),
        n_residues = as.integer(scalar(summary, "n_residues", 0L)),
        residues = scalar(summary, "residues"),
        p_1d_joint = as.numeric(scalar(summary, "p_1d_joint", NA_real_)),
        p_3d_joint = as.numeric(scalar(summary, "p_3d_joint", NA_real_)),
        p_any_joint = as.numeric(scalar(summary, "p_any_joint", NA_real_)),
        n_annotation_features = as.integer(scalar(summary, "n_annotation_features", 0L)),
        n_annotated_cluster_residues = as.integer(scalar(summary, "n_annotated_cluster_residues", 0L)),
        annotation_classes = scalar(summary, "annotation_classes"),
        feature_subcategories = scalar(summary, "feature_subcategories")
      ),
      annotations = annotations
    )
  })
  payload <- list(
    protein = annotation.result$protein,
    accession = annotation.result$accession,
    protein_length = annotation.result$protein_length,
    position_label = annotation.result$position_label %||% "cluster residues",
    clusters = cluster.records
  )
  json <- jsonlite::toJSON(
    payload, auto_unbox = TRUE, dataframe = "rows", null = "null",
    na = "null", digits = 10
  )
  json <- gsub("</", "<\\/", json, fixed = TRUE)

  css <- paste0(
    ":root{--ink:#172033;--muted:#64748b;--line:#dbe3ec;--panel:#fff;--bg:#f5f7fa;--accent:#7a0019}",
    "*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 system-ui,-apple-system,sans-serif}",
    "header{background:#172033;color:#fff;padding:1.2rem 2rem}header h1{margin:0;font-size:1.5rem}header p{margin:.3rem 0 0;color:#cbd5e1}",
    ".layout{display:grid;grid-template-columns:330px minmax(0,1fr);min-height:calc(100vh - 86px)}",
    "aside{border-right:1px solid var(--line);background:#fff;padding:1rem;overflow:auto}main{padding:1rem;overflow:auto}",
    ".controls{display:grid;gap:.65rem}.controls label{font-weight:650;font-size:.82rem}input,select,button{width:100%;padding:.55rem;border:1px solid #bcc8d6;border-radius:6px;background:#fff}",
    ".nav{display:grid;grid-template-columns:1fr 1fr;gap:.5rem}.cluster-list{margin-top:1rem;max-height:60vh;overflow:auto;border:1px solid var(--line);border-radius:7px}",
    ".cluster-row{padding:.6rem;border-bottom:1px solid var(--line);cursor:pointer}.cluster-row:hover,.cluster-row.active{background:#eef5ff}.cluster-row small{display:block;color:var(--muted)}",
    ".cards{display:grid;grid-template-columns:repeat(6,minmax(115px,1fr));gap:.65rem;margin-bottom:1rem}.card,.panel{background:var(--panel);border:1px solid var(--line);border-radius:8px;box-shadow:0 1px 2px #00000008}",
    ".card{padding:.75rem}.card span{display:block;color:var(--muted);font-size:.75rem;text-transform:uppercase}.card strong{font-size:1.1rem}",
    ".panel{padding:1rem;margin-bottom:1rem}.panel h2{margin:0 0 .65rem;font-size:1.05rem}.note{border-left:4px solid #d7a900;background:#fff9df;padding:.7rem;margin-bottom:1rem}",
    ".plot-wrap{overflow-x:auto}svg{min-width:1050px;width:100%;height:auto}.axis{stroke:#52606d}.grid{stroke:#e7edf3;stroke-dasharray:3 3}.residue{fill:#7a0019;stroke:#fff;stroke-width:1.4}",
    "table{border-collapse:collapse;width:100%;font-size:.82rem}th,td{border:1px solid var(--line);padding:.45rem;vertical-align:top}th{background:#eef2f6;text-align:left;position:sticky;top:0}tr:nth-child(even){background:#fafbfd}",
    ".table-wrap{max-height:520px;overflow:auto}.pill{display:inline-block;padding:.16rem .48rem;border-radius:999px;background:#e8eef5;margin:.1rem;font-size:.78rem}.sig{background:#dff6e7;color:#176b38}.nonsig{background:#f1f3f5;color:#52606d}",
    "a{color:#075985}.empty{color:var(--muted)}@media(max-width:1000px){.layout{grid-template-columns:1fr}.cards{grid-template-columns:repeat(2,1fr)}aside{border-right:0;border-bottom:1px solid var(--line)}}"
  )

  javascript <- paste0(
    "const DATA=", json, ";const NS='http://www.w3.org/2000/svg';",
    "const COLORS={domain:'#4477AA',site:'#228833',modification:'#EE3377',interaction:'#66CCEE',motif:'#AA3377',structural_region:'#999999',functional_region:'#999933'};",
    "const $=id=>document.getElementById(id);const esc=s=>String(s??'').replace(/[&<>\"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','\"':'&quot;',\"'\":'&#39;'}[c]));",
    "const fmtP=x=>x==null?'NA':(x<0.001?Number(x).toExponential(2):Number(x).toFixed(4));",
    "let visible=DATA.clusters.map((_,i)=>i),current=0;",
    "function selectedAnnotations(c){const ac=$('annotationFilter').value,sc=$('subclassFilter').value;return c.annotations.filter(a=>(!ac||a.annotation_class===ac)&&(!sc||a.feature_subcategory===sc));}",
    "function filtered(){const q=$('search').value.toLowerCase(),sig=$('sigFilter').value,h=$('classFilter').value,annotationActive=$('annotationFilter').value||$('subclassFilter').value;return DATA.clusters.map((c,i)=>[c,i]).filter(([c])=>{const s=c.summary;return(!q||s.cluster_id.toLowerCase().includes(q)||s.residues.includes(q)||s.hotspot_class.toLowerCase().includes(q))&&(sig==='all'||String(s.significant_any)===sig)&&(!h||s.hotspot_class===h)&&(!annotationActive||selectedAnnotations(c).length)}).map(x=>x[1]);}",
    "function refreshList(keep=true){visible=filtered();const list=$('clusterList');list.innerHTML='';visible.forEach(i=>{const s=DATA.clusters[i].summary,d=document.createElement('div');d.className='cluster-row'+(i===current?' active':'');d.innerHTML='<b>'+esc(s.cluster_id)+'</b><small>'+esc(s.hotspot_class)+'</small><small>'+s.n_residues+' residues · '+s.n_annotation_features+' features</small>';d.onclick=()=>selectCluster(i);list.appendChild(d)});if(!visible.includes(current)&&visible.length)selectCluster(visible[0],false);$('resultCount').textContent=visible.length+' clusters';}",
    "function selectCluster(i,refresh=true){current=i;render();if(refresh)refreshList(false);}",
    "function move(delta){const p=visible.indexOf(current);if(p<0||!visible.length)return;selectCluster(visible[(p+delta+visible.length)%visible.length]);}",
    "function svgEl(tag,a={}){const e=document.createElementNS(NS,tag);Object.entries(a).forEach(([k,v])=>e.setAttribute(k,v));return e;}",
    "function renderPlot(c){const s=c.summary,a=c.annotations,W=1400,L=410,R=35,row=31,top=82,H=Math.max(210,top+a.length*row+55),x=v=>L+(Number(v)-1)/(DATA.protein_length-1)*(W-L-R),svg=$('map');svg.setAttribute('viewBox',`0 0 ${W} ${H}`);svg.innerHTML='';const tick=DATA.protein_length<=500?50:100;for(let v=0;v<=DATA.protein_length;v+=tick){svg.appendChild(svgEl('line',{x1:x(Math.max(1,v)),x2:x(Math.max(1,v)),y1:35,y2:H-30,class:'grid'}));let t=svgEl('text',{x:x(Math.max(1,v)),y:H-10,'text-anchor':'middle',fill:'#52606d'});t.textContent=v;svg.appendChild(t)}svg.appendChild(svgEl('line',{x1:x(1),x2:x(DATA.protein_length),y1:42,y2:42,stroke:'#374151','stroke-width':4}));let residues=s.residues.split(';').map(Number);residues.forEach(r=>{svg.appendChild(svgEl('line',{x1:x(r),x2:x(r),y1:26,y2:42,stroke:'#7a0019'}));svg.appendChild(svgEl('circle',{cx:x(r),cy:24,r:4,class:'residue'}))});let title=svgEl('text',{x:L-15,y:28,'text-anchor':'end',fill:'#7a0019','font-weight':650});title.textContent=DATA.position_label;svg.appendChild(title);a.forEach((d,i)=>{const y=top+i*row,color=COLORS[d.annotation_class]||'#777';const label=`[${d.annotation_class} › ${d.feature_subcategory}] ${d.source}: ${d.feature_name}`;let tx=svgEl('text',{x:L-15,y:y+5,'text-anchor':'end',fill:'#253247'});tx.textContent=label.length>67?label.slice(0,64)+'...':label;svg.appendChild(tx);d.segments.split(';').forEach(seg=>{const b=seg.split('-').map(Number);if(b.length!==2||!b.every(Number.isFinite))return;if(b[0]===b[1])svg.appendChild(svgEl('rect',{x:x(b[0])-3,y:y-6,width:6,height:12,fill:color,opacity:.85}));else svg.appendChild(svgEl('rect',{x:x(b[0]),y:y-9,width:Math.max(2,x(b[1])-x(b[0])),height:18,fill:color,opacity:.9}))});d.overlap_residues.split(';').filter(Boolean).map(Number).forEach(r=>svg.appendChild(svgEl('circle',{cx:x(r),cy:y,r:4,class:'residue'})))});svg.setAttribute('height',H);}",
    "function literature(a){if(!a.literature.length)return '<span class=empty>—</span>';return a.literature.slice(0,5).map(p=>'<a target=_blank title=\"'+esc(p.title)+' ['+esc(p.scope)+']\" href=\"'+esc(p.url)+'\">'+esc(p.label)+'</a>').join(', ')+(a.literature.length>5?' (+'+(a.literature.length-5)+' more)':'');}",
    "function renderTable(c){$('annotationBody').innerHTML=c.annotations.map(a=>'<tr><td><span class=pill>'+esc(a.annotation_class)+'</span></td><td>'+esc(a.feature_subcategory)+'</td><td><b>'+esc(a.feature_name)+'</b><br><small>'+esc(a.feature_type)+' · '+esc(a.segments)+'</small></td><td>'+esc(a.overlap_residues)+'</td><td>'+esc(a.relationship)+'</td><td><a target=_blank href=\"'+esc(a.source_url)+'\">'+esc(a.source)+' '+esc(a.source_record_id)+'</a></td><td>'+esc(a.evidence_codes||'—')+'</td><td>'+literature(a)+'</td></tr>').join('');}",
    "function render(){const c=DATA.clusters[current],s=c.summary,a=selectedAnnotations(c),view={summary:s,annotations:a};$('clusterTitle').textContent=s.cluster_id+' — '+s.hotspot_class;$('clusterMeta').innerHTML='<span class=\"pill '+(s.significant_any?'sig':'nonsig')+'\">'+(s.significant_any?'statistically significant':'not statistically significant')+'</span><span class=pill>trees: '+esc(s.tree_sources)+'</span><span class=pill>'+esc(s.annotation_classes)+'</span>';$('cards').innerHTML=[['residues',s.n_residues],['subjects',s.n_subjects],['events',s.n_events],['features shown',a.length+' / '+s.n_annotation_features],['joint 1D p',fmtP(s.p_1d_joint)],['joint 3D p',fmtP(s.p_3d_joint)]].map(x=>'<div class=card><span>'+x[0]+'</span><strong>'+x[1]+'</strong></div>').join('');$('residueText').textContent=s.residues;renderPlot(view);renderTable(view);}",
    "function updateSubclasses(){const ac=$('annotationFilter').value,currentValue=$('subclassFilter').value,values=[...new Set(DATA.clusters.flatMap(c=>c.annotations).filter(a=>!ac||a.annotation_class===ac).map(a=>a.feature_subcategory))].sort(),select=$('subclassFilter');select.innerHTML='<option value=\"\">All subclasses</option>';values.forEach(v=>{let o=document.createElement('option');o.value=v;o.textContent=v;select.appendChild(o)});if(values.includes(currentValue))select.value=currentValue;}",
    "function init(){[...new Set(DATA.clusters.map(c=>c.summary.hotspot_class))].sort().forEach(v=>{let o=document.createElement('option');o.value=v;o.textContent=v;$('classFilter').appendChild(o)});[...new Set(DATA.clusters.flatMap(c=>c.annotations).map(a=>a.annotation_class))].sort().forEach(v=>{let o=document.createElement('option');o.value=v;o.textContent=v;$('annotationFilter').appendChild(o)});updateSubclasses();['search','sigFilter','classFilter','subclassFilter'].forEach(id=>$(id).addEventListener(id==='search'?'input':'change',()=>{refreshList();render()}));$('annotationFilter').addEventListener('change',()=>{updateSubclasses();refreshList();render()});$('prev').onclick=()=>move(-1);$('next').onclick=()=>move(1);refreshList();render();}init();"
  )

  chembl <- annotation.result$chembl_mechanisms %||% grin3d_empty_table(
    grin3d_chembl_mechanism_columns
  )
  chembl.rows <- if (nrow(chembl)) vapply(seq_len(nrow(chembl)), function(i) {
    row <- chembl[i, , drop = FALSE]
    paste0(
      "<tr><td><a target=_blank href=\"", grin3d_html_escape(row$molecule_url), "\">",
      grin3d_html_escape(if (nzchar(row$molecule_name)) row$molecule_name else row$molecule_chembl_id),
      "</a><br><small>", grin3d_html_escape(row$molecule_chembl_id), "</small></td><td>",
      grin3d_html_escape(row$action_type), "</td><td>",
      grin3d_html_escape(row$mechanism_of_action), "</td><td><a target=_blank href=\"",
      grin3d_html_escape(row$target_url), "\">", grin3d_html_escape(row$target_name),
      "</a><br><small>", grin3d_html_escape(row$target_type), "</small></td><td>",
      grin3d_html_escape(row$max_phase), "</td></tr>"
    )
  }, character(1L)) else '<tr><td colspan="5" class="empty">No ChEMBL mechanism records returned.</td></tr>'
  chembl.panel <- paste0(
    '<section class="panel"><h2>ChEMBL target-level drug context</h2>',
    '<p class="note"><b>Not residue evidence:</b> these mechanisms apply to a protein or ',
    'complex target containing this accession. They are not assigned to the selected cluster.</p>',
    '<div class="table-wrap"><table><thead><tr><th>Molecule</th><th>Action</th>',
    '<th>Mechanism</th><th>Target</th><th>Maximum phase</th></tr></thead><tbody>',
    paste(chembl.rows, collapse = ""), '</tbody></table></div></section>'
  )

  html <- c(
    "<!doctype html><html><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
    paste0("<title>", grin3d_html_escape(annotation.result$protein), " cluster annotation dashboard</title><style>", css, "</style></head><body>"),
    paste0("<header><h1>", grin3d_html_escape(annotation.result$protein), " cluster annotation dashboard</h1><p>UniProt ", grin3d_html_escape(annotation.result$accession), " · statistical results and post-hoc biological annotations</p></header>"),
    "<div class=\"layout\"><aside><div class=\"controls\"><label>Search clusters or residues<input id=\"search\" placeholder=\"cluster_0001 or 248\"></label><label>Statistical status<select id=\"sigFilter\"><option value=\"all\">All clusters</option><option value=\"true\">Significant</option><option value=\"false\">Not significant</option></select></label><label>Hotspot class<select id=\"classFilter\"><option value=\"\">All hotspot classes</option></select></label><label>Annotation class<select id=\"annotationFilter\"><option value=\"\">All annotation classes</option></select></label><label>Annotation subclass<select id=\"subclassFilter\"><option value=\"\">All subclasses</option></select></label><div class=\"nav\"><button id=\"prev\">&#8592; Previous</button><button id=\"next\">Next &#8594;</button></div><b id=\"resultCount\"></b></div><div id=\"clusterList\" class=\"cluster-list\"></div></aside>",
    "<main><div class=\"note\"><b>Interpretation boundary:</b> annotations are exact coordinate overlaps added after GRIN3D testing. They do not alter cluster membership or statistical significance. InterPro entry literature describes the feature, not necessarily this cluster.</div><h2 id=\"clusterTitle\"></h2><div id=\"clusterMeta\"></div><div id=\"cards\" class=\"cards\"></div>",
    paste0("<section class=\"panel\"><h2>Residue and annotation map</h2><p><b>",
      grin3d_html_escape(annotation.result$position_label %||% "Cluster residues"),
      ":</b> <code id=\"residueText\"></code></p><div class=\"plot-wrap\"><svg id=\"map\" role=\"img\" aria-label=\"Cluster residue and annotation map\"></svg></div></section>"),
    "<section class=\"panel\"><h2>Annotation evidence</h2><div class=\"table-wrap\"><table><thead><tr><th>Class</th><th>Subclass</th><th>Source-native feature</th><th>Overlapping residues</th><th>Relationship</th><th>Source record</th><th>Evidence</th><th>Literature</th></tr></thead><tbody id=\"annotationBody\"></tbody></table></div></section>",
    chembl.panel, "</main></div>",
    paste0("<script>", javascript, "</script></body></html>")
  )
  dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, output.file, useBytes = TRUE)
  message("Single-file cluster annotation dashboard written to: ", output.file)
  invisible(output.file)
}
