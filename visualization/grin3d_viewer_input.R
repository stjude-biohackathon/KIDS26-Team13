# GRIN3D viewer input format
# ==========================
# The viewer reads, per protein, three tables named <PROTEIN>_viewer_*.csv in
# visualization/data/, plus one manifest listing the proteins:
#
# proteins.csv           -- one row per protein the viewer can show
#   protein, uniprot, pdb_file (path relative to the repository root),
#   types (";"-separated alteration types present for this protein)
#
# <prefix>_clusters.csv  -- one row per candidate cluster
#   protein, source ("mutation" | "cna"),
#   alteration_type (SNV | INDEL | MUT | HOMDEL | HETDEL | GAIN | AMP | ALL_CNA),
#   MUT is for mutation results that were not split into SNVs and indels.
#   cluster_id, hotspot_class, n_subjects, n_events, n_residues,
#   diameter_1d, diameter_3d, p_1d_size, p_3d_size (size-specific),
#   p_1d_joint, p_3d_joint, p_any_joint (protein-wide joint),
#   significant_any, mean_plddt, exons (";"-separated exon orders, NA for mutations),
#   significant_1d, significant_3d, p_1d_protein, p_3d_protein, tree_sources,
#   residue_min, residue_max (residue range the cluster covers on the protein)
#
# <prefix>_residues.csv  -- one row per residue x alteration_type x cluster
#   protein, residue, alteration_type,
#   cluster_id (NA = background count of all events at that residue/exon),
#   exon_order, n_events, n_subjects, plddt
#
# <prefix>_members.csv   -- one row per cluster x contributing event
#   cluster_id, alteration_type, event_id, subject_id, event_type, detail,
#   coding_fraction (CNA events: share of the event's coding bases inside the
#   cluster's exons; NA for mutation events)
#
# Mutation cluster ids are prefixed with their type ("SNV_", "INDEL_", "MUT_"; CNA
# ids are already typed) so ids stay unique once all results are combined.

split_int <- function(x) as.integer(strsplit(as.character(x), ";", fixed = TRUE)[[1]])

standardize_mutation_results <- function(res, type = "MUT") {
  cl <- res$clusters
  cl$cluster_id <- paste0(type, "_", cl$cluster_id)

  clusters <- data.frame(
    protein = cl$protein, source = "mutation", alteration_type = type,
    cluster_id = cl$cluster_id, hotspot_class = cl$hotspot_class,
    n_subjects = cl$n_subjects, n_events = cl$n_events, n_residues = cl$n_residues,
    diameter_1d = cl$diameter_1d, diameter_3d = cl$diameter_3d,
    p_1d_size = cl$p_1d_size, p_3d_size = cl$p_3d_size,
    p_1d_joint = cl$p_1d_joint, p_3d_joint = cl$p_3d_joint,
    p_any_joint = cl$p_any_joint, significant_any = cl$significant_any,
    exons = NA_character_,
    significant_1d = cl$significant_1d, significant_3d = cl$significant_3d,
    p_1d_protein = cl$p_1d_protein, p_3d_protein = cl$p_3d_protein,
    tree_sources = cl$tree_sources,
    residue_min = cl$residue_min, residue_max = cl$residue_max
  )

  ev <- res$event_residue_mapping
  background <- do.call(rbind, lapply(split(ev, ev$residue), function(d) data.frame(
    residue = d$residue[1], n_events = length(unique(d$event_id)),
    n_subjects = length(unique(d$subject_id)))))
  background$cluster_id <- NA_character_

  in_cluster <- do.call(rbind, lapply(seq_len(nrow(cl)), function(i) {
    r <- split_int(cl$residues[i])
    data.frame(residue = r, cluster_id = cl$cluster_id[i],
               n_events = background$n_events[match(r, background$residue)],
               n_subjects = background$n_subjects[match(r, background$residue)])
  }))

  residues <- rbind(background, in_cluster)
  residues$protein <- cl$protein[1]
  residues$alteration_type <- type
  residues$exon_order <- NA_integer_

  m <- res$cluster_members
  members <- data.frame(
    cluster_id = paste0(type, "_", m$cluster_id), alteration_type = type,
    event_id = m$event_id, subject_id = m$subject_id, event_type = type,
    detail = paste0("residue ", m$affected_residues_in_cluster),
    coding_fraction = NA_real_)

  list(clusters = clusters, residues = residues, members = members)
}

standardize_cna_results <- function(res) {
  cl <- res$clusters
  exon_res <- res$exon_residue_mapping[, c("exon_order", "residue")]

  clusters <- data.frame(
    protein = cl$protein, source = "cna", alteration_type = cl$analysis_type,
    cluster_id = cl$cluster_id, hotspot_class = cl$hotspot_class,
    n_subjects = cl$n_subjects, n_events = cl$n_events,
    n_residues = vapply(cl$exon_orders, function(e)
      length(unique(exon_res$residue[exon_res$exon_order %in% split_int(e)])), 1L),
    diameter_1d = cl$diameter_1d, diameter_3d = cl$diameter_3d,
    p_1d_size = cl$p_1d_size, p_3d_size = cl$p_3d_size,
    p_1d_joint = cl$p_1d_joint, p_3d_joint = cl$p_3d_joint,
    p_any_joint = cl$p_any_joint, significant_any = cl$significant_any,
    exons = cl$exon_orders,
    significant_1d = cl$significant_1d, significant_3d = cl$significant_3d,
    p_1d_protein = cl$p_1d_protein, p_3d_protein = cl$p_3d_protein,
    tree_sources = cl$tree_sources,
    residue_min = vapply(cl$exon_orders, function(e)
      min(exon_res$residue[exon_res$exon_order %in% split_int(e)]), 1L),
    residue_max = vapply(cl$exon_orders, function(e)
      max(exon_res$residue[exon_res$exon_order %in% split_int(e)]), 1L)
  )

  # Background: per-exon event counts for each CNA type, spread over the exon's residues.
  counts <- res$exon_counts_by_type
  types <- sub("^n_events_1d_", "", grep("^n_events_1d_", names(counts), value = TRUE))
  background <- do.call(rbind, lapply(types, function(t) data.frame(
    exon_order = counts$exon_order, alteration_type = t, cluster_id = NA_character_,
    n_events = counts[[paste0("n_events_1d_", t)]],
    n_subjects = counts[[paste0("n_subjects_1d_", t)]])))
  background <- background[background$n_events > 0, ]

  # Cluster rows carry that exon's per-exon counts from events_per_exon ("2:26;5:19").
  per_exon <- function(s, e) {
    kv <- do.call(rbind, strsplit(strsplit(s, ";", fixed = TRUE)[[1]], ":", fixed = TRUE))
    as.integer(kv[match(e, as.integer(kv[, 1])), 2])
  }
  in_cluster <- do.call(rbind, lapply(seq_len(nrow(cl)), function(i) {
    e <- split_int(cl$exon_orders[i])
    data.frame(exon_order = e, alteration_type = cl$analysis_type[i],
               cluster_id = cl$cluster_id[i],
               n_events = per_exon(cl$events_per_exon[i], e),
               n_subjects = per_exon(cl$subjects_per_exon[i], e))
  }))

  residues <- merge(rbind(background, in_cluster), exon_res, by = "exon_order")
  residues$protein <- cl$protein[1]

  # The members table also lists overlapping CNAs of other types; keep only the
  # events that count toward the cluster (their number equals n_events).
  m <- res$cluster_members
  m <- m[m$supports_analysis_type, ]
  members <- data.frame(
    cluster_id = m$cluster_id, alteration_type = m$analysis_type,
    event_id = m$event_id, subject_id = m$subject_id, event_type = m$cna_type,
    detail = paste0("exons ", m$affected_exon_orders_in_cluster),
    coding_fraction = round(m$event_coding_fraction, 3))

  list(clusters = clusters, residues = residues, members = members)
}

# Per-residue pLDDT (C-alpha B-factor column) and C-alpha coordinates from an AlphaFold PDB file.
read_plddt <- function(pdb_file) {
  p <- readLines(pdb_file, warn = FALSE)
  ca <- p[startsWith(p, "ATOM") & substr(p, 13, 16) == " CA "]
  data.frame(residue = as.integer(substr(ca, 23, 26)),
             plddt = as.numeric(substr(ca, 61, 66)),
             x = as.numeric(substr(ca, 31, 38)),
             y = as.numeric(substr(ca, 39, 46)),
             z = as.numeric(substr(ca, 47, 54)))
}

# Writes <data_dir>/<protein>_viewer_*.csv and adds/replaces the protein's row in
# <data_dir>/proteins.csv. mutation_rds is a named vector of mutation results RDS
# files, names = type (SNV, INDEL, or MUT if unsplit); either argument may be NULL.
build_viewer_input <- function(protein, uniprot, pdb_file, mutation_rds = NULL,
                               cna_rds = NULL, data_dir = "visualization/data") {
  parts <- c(
    lapply(names(mutation_rds), function(t) standardize_mutation_results(readRDS(mutation_rds[[t]]), t)),
    if (!is.null(cna_rds)) list(standardize_cna_results(readRDS(cna_rds)))
  )
  stopifnot(length(parts) > 0)
  clusters <- do.call(rbind, lapply(parts, `[[`, "clusters"))
  residue_cols <- c("protein", "residue", "alteration_type", "cluster_id",
                    "exon_order", "n_events", "n_subjects")
  residues <- do.call(rbind, lapply(parts, function(p) p$residues[, residue_cols]))
  members <- do.call(rbind, lapply(parts, `[[`, "members"))

  plddt <- read_plddt(pdb_file)
  residues$plddt <- plddt$plddt[match(residues$residue, plddt$residue)]
  plddt_by_cluster <- tapply(residues$plddt, residues$cluster_id,
                             function(p) mean(unique(p), na.rm = TRUE))
  clusters$mean_plddt <- round(as.numeric(plddt_by_cluster[clusters$cluster_id]), 1)

  stopifnot(
    all(clusters$protein == protein),
    all(residues$residue %in% plddt$residue),
    !anyDuplicated(clusters$cluster_id),
    all(residues$cluster_id[!is.na(residues$cluster_id)] %in% clusters$cluster_id),
    all(members$cluster_id %in% clusters$cluster_id)
  )

  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
  prefix <- file.path(data_dir, paste0(protein, "_viewer"))
  write.csv(clusters, paste0(prefix, "_clusters.csv"), row.names = FALSE)
  write.csv(residues, paste0(prefix, "_residues.csv"), row.names = FALSE)
  write.csv(members, paste0(prefix, "_members.csv"), row.names = FALSE)

  manifest_file <- file.path(data_dir, "proteins.csv")
  present <- unique(c(residues$alteration_type, clusters$alteration_type))
  row <- data.frame(protein = protein, uniprot = uniprot, pdb_file = pdb_file,
                    types = paste(present, collapse = ";"))
  manifest <- if (file.exists(manifest_file)) read.csv(manifest_file) else row[0, ]
  manifest <- rbind(manifest[manifest$protein != protein, names(row)], row)
  write.csv(manifest[order(manifest$protein), ], manifest_file, row.names = FALSE)

  invisible(list(clusters = clusters, residues = residues, members = members))
}

# Biological annotations (UniProt, InterPro, PDBe-KB/SIFTS, ChEMBL) for every
# cluster of a protein, using annotate_grin3d_clusters() from
# development-code/GRIN3D_annotate_protein_clusters.R (PR #1). Runs on the combined
# viewer cluster table (SNV, indel and CNA clusters; CNA clusters carry the residues
# of their exons) and writes that module's standard outputs to
# <data_dir>/annotations/<protein>/. API responses are cached under raw/.
build_viewer_annotations <- function(protein, uniprot, data_dir = "visualization/data",
                                     offline = FALSE) {
  module <- new.env()
  sys.source("development-code/GRIN3D_annotate_protein_clusters.R", envir = module)
  prefix <- file.path(data_dir, paste0(protein, "_viewer"))
  cl <- read.csv(paste0(prefix, "_clusters.csv"))
  res <- read.csv(paste0(prefix, "_residues.csv"))
  res <- res[!is.na(res$cluster_id), ]
  cl$residues <- vapply(cl$cluster_id, function(id)
    paste(sort(unique(res$residue[res$cluster_id == id])), collapse = ";"), "")
  out_dir <- file.path(data_dir, "annotations", protein)
  dir.create(file.path(out_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
  # PDBe-KB answers 404 when a protein has no interface or ligand records
  # (e.g. CDKN2A ligands); the module treats that as a failure, so cache an
  # empty response for those endpoints first.
  if (!offline) for (kind in c("interfaces", "ligands")) {
    cache <- file.path(out_dir, "raw", sprintf("pdbe_%s_%s.json", kind, uniprot))
    url <- sprintf(if (kind == "ligands") module$grin3d_annotation_urls$pdbe.ligands
                   else module$grin3d_annotation_urls$pdbe.interfaces, uniprot)
    if (file.exists(cache)) next
    status <- system2("curl", c("-s", "-o", "/dev/null", "-w", "%{http_code}", shQuote(url)), stdout = TRUE)
    if (identical(status, "404")) writeLines("{}", cache)
  }
  input <- file.path(out_dir, "viewer_annotation_clusters.csv")
  write.csv(cl[, intersect(module$grin3d_cluster_columns, names(cl))], input, row.names = FALSE)
  module$annotate_grin3d_clusters(
    cluster.file = input, uniprot.accession = uniprot, output.dir = out_dir,
    protein = protein, offline = offline)
}

if (sys.nframe() == 0L) {
  # Rscript visualization/grin3d_viewer_input.R   (from the repository root)
  # To add a protein: add a build_viewer_input() call pointing at its results.
  build_viewer_input(
    protein = "PTEN", uniprot = "P60484",
    pdb_file = "examples/PTEN_CNA/input_files/AF-P60484-F1-model_v6.pdb",
    mutation_rds = c(
      SNV = "examples/PTEN_mutations/results_SNV/GRIN3D_mutation_hotspot_results.rds",
      INDEL = "examples/PTEN_mutations/results_INDEL/GRIN3D_mutation_hotspot_results.rds"),
    cna_rds = "examples/PTEN_CNA/results/GRIN3D_exon_CNA_results.rds"
  )
  build_viewer_input(
    protein = "CDKN2A", uniprot = "P42771",
    pdb_file = "examples/PTEN_CNA/input_files/AF-P42771-F1-model_v6.pdb",
    cna_rds = "examples/CDKN2A_CNA/results/GRIN3D_exon_CNA_results.rds"
  )
  build_viewer_input(
    protein = "LEF1", uniprot = "Q9UJU2",
    pdb_file = "examples/LEF1_CNA/input_files/AF-Q9UJU2-F1-model_v6.pdb",
    mutation_rds = c(
      SNV = "examples/LEF1_mutations/results_SNV/GRIN3D_mutation_hotspot_results.rds",
      INDEL = "examples/LEF1_mutations/results_INDEL/GRIN3D_mutation_hotspot_results.rds"),
    cna_rds = "examples/LEF1_CNA/results/GRIN3D_exon_CNA_results.rds"
  )
  # Mutation-only genes from visualization/prepare_mutation_examples.R
  mutation_genes <- c(TP53 = "P04637", FBXW7 = "Q969H0", JAK3 = "P52333", IL7R = "P16871")
  for (g in names(mutation_genes)) {
    dir <- file.path("examples", paste0(g, "_mutations"))
    build_viewer_input(
      protein = g, uniprot = mutation_genes[[g]],
      pdb_file = file.path(dir, "input_files", sprintf("AF-%s-F1-model_v6.pdb", mutation_genes[[g]])),
      mutation_rds = c(
        SNV = file.path(dir, "results_SNV", "GRIN3D_mutation_hotspot_results.rds"),
        INDEL = file.path(dir, "results_INDEL", "GRIN3D_mutation_hotspot_results.rds"))
    )
  }
  proteins <- read.csv("visualization/data/proteins.csv")
  if (file.exists("development-code/GRIN3D_annotate_protein_clusters.R")) {
    for (i in seq_len(nrow(proteins))) build_viewer_annotations(proteins$protein[i], proteins$uniprot[i])
  }
  print(proteins)
}
