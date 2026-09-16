source("development-code/GRIN3D_annotate_protein_clusters.R")

stopifnot(identical(grin3d_parse_residues("248;244;248"), c(244L, 248L)))
stopifnot(grin3d_feature_category("Modified residue") == "ptm")
stopifnot(grin3d_feature_category("Region", "Interaction with MDM2") == "interaction_region")
stopifnot(grin3d_feature_category("Region", "Disordered") == "disorder")
classification <- grin3d_feature_classification("Binding site", "Zn(2+)")
stopifnot(classification$annotation_class == "site")
stopifnot(classification$feature_subcategory == "metal_binding_site")
classification <- grin3d_feature_classification(
  "Modified residue", "Phosphoserine; by ATM"
)
stopifnot(classification$annotation_class == "modification")
stopifnot(classification$feature_subcategory == "phosphorylation")
classification <- grin3d_feature_classification(
  "family", "p53 tumour suppressor family", "DNA-binding proteins"
)
stopifnot(classification$feature_subcategory == "protein_family")
classification <- grin3d_feature_classification(
  "domain", "p53, transactivation domain", "Participates in DNA binding"
)
stopifnot(classification$feature_subcategory == "transactivation_domain")
classification <- grin3d_feature_classification(
  "domain", "p53, tetramerisation domain", "Supports DNA binding"
)
stopifnot(classification$feature_subcategory == "oligomerization_domain")
stopifnot(grin3d_normalize_doi("https://doi.org/10.1038/example") == "10.1038/example")
stopifnot(grin3d_normalize_doi("https://example.org/paper.pdf") == "")

clusters <- data.frame(
  protein = "TEST", cluster_id = "cluster_1", residues = "10;20;30",
  stringsAsFactors = FALSE
)
for (column in setdiff(grin3d_cluster_columns, names(clusters))) clusters[[column]] <- ""
clusters <- clusters[, grin3d_cluster_columns]
features <- list(
  list(
    feature_id = "domain", source = "TestDB", source_record_id = "D1",
    annotation_class = "domain", feature_category = "domain",
    feature_subcategory = "protein_domain",
    classification_rule = "grin3d_annotation_hierarchy_v1:protein_domain",
    feature_type = "domain",
    feature_name = "test domain", start = 15L, end = 25L,
    feature_segments = "15-25", evidence_codes = "ECO:test",
    evidence_ids = "1", source_url = "https://example.test/domain",
    segments = list(c(15L, 25L))
  ),
  list(
    feature_id = "site", source = "TestDB", source_record_id = "S1",
    annotation_class = "site", feature_category = "binding_site",
    feature_subcategory = "ligand_binding_site",
    classification_rule = "grin3d_annotation_hierarchy_v1:ligand_binding_site",
    feature_type = "site",
    feature_name = "test site", start = 30L, end = 30L,
    feature_segments = "30-30", evidence_codes = "", evidence_ids = "",
    source_url = "https://example.test/site", segments = list(c(30L, 30L))
  )
)
mapped <- grin3d_map_annotations(clusters, list(c(10L, 20L, 30L)), features, "P00000")
stopifnot(length(mapped$overlaps) == 2L)
stopifnot(mapped$overlaps[[1L]]$overlap_residues == "20")
stopifnot(mapped$overlaps[[1L]]$relationship == "interval_overlap")
stopifnot(mapped$overlaps[[2L]]$relationship == "exact_residue_overlap")
stopifnot(mapped$summaries[[1L]]$unannotated_cluster_residues == "10")

pdbe.interface.mock <- list(P00000 = list(data = list(list(
  name = "Partner", accession = "P99999",
  residues = list(list(
    startIndex = 20L, endIndex = 20L, indexType = "UNIPROT",
    interactingPDBEntries = list(list(pdbId = "1abc"))
  )),
  additionalData = list(type = "UNP")
))))
pdbe.features <- grin3d_extract_pdbe_interfaces(
  pdbe.interface.mock, "TEST", "P00000", "2026-01-01T00:00:00Z"
)
stopifnot(length(pdbe.features) == 1L)
stopifnot(pdbe.features[[1L]]$feature_subcategory == "protein_protein_interface")
stopifnot(pdbe.features[[1L]]$relationship_type == "experimental_interface_contact")
stopifnot(pdbe.features[[1L]]$evidence_ids == "1ABC")

pdbe.ligand.mock <- list(P00000 = list(data = list(list(
  name = "Compound", accession = "ABC",
  residues = list(list(startIndex = 30L, endIndex = 30L, indexType = "UNIPROT")),
  additionalData = list(drugBankId = "DB00001", isSolvent = FALSE)
))))
pdbe.ligands <- grin3d_extract_pdbe_ligands(
  pdbe.ligand.mock, "TEST", "P00000", "2026-01-01T00:00:00Z"
)
stopifnot(pdbe.ligands$features[[1L]]$feature_subcategory == "drug_linked_ligand_site")
stopifnot(pdbe.ligands$ligands[[1L]]$evidence_scope == "residue_level_experimental_contact")

chembl <- grin3d_extract_chembl(
  list(targets = list(list(
    target_chembl_id = "CHEMBL1", pref_name = "Target", target_type = "SINGLE PROTEIN",
    target_components = list(list(accession = "P00000", relationship = "SINGLE PROTEIN"))
  ))),
  list(mechanisms = list(list(
    target_chembl_id = "CHEMBL1", molecule_chembl_id = "CHEMBL2",
    action_type = "INHIBITOR", mechanism_of_action = "Target inhibitor"
  ))),
  list(molecules = list(list(
    molecule_chembl_id = "CHEMBL2", pref_name = "Drug", max_phase = 3
  ))),
  "TEST", "P00000", "2026-01-01T00:00:00Z"
)
stopifnot(length(chembl$targets) == 1L, length(chembl$mechanisms) == 1L)
stopifnot(chembl$mechanisms[[1L]]$evidence_scope == "target_level_drug_mechanism_not_residue_evidence")

cna.temp <- tempfile(fileext = ".csv")
exon.temp <- tempfile(fileext = ".csv")
annotation.temp <- tempfile(fileext = ".csv")
utils::write.csv(data.frame(
  protein = "TEST", transcript = "TX1", cluster_id = "CNA_1",
  analysis_type = "HOMDEL", candidate_key = "1;2", tree_sources = "3D",
  exon_ids = "E1;E2", exon_orders = "1;2", hotspot_class = "test",
  significant_1d = FALSE, significant_3d = TRUE, significant_any = TRUE,
  n_subjects = 2L, n_events = 2L, p_1d_joint = 0.2,
  p_3d_joint = 0.01, p_any_joint = 0.01
), cna.temp, row.names = FALSE)
utils::write.csv(data.frame(
  exon_id = c("E1", "E2"), aa_start = c(1L, 3L), aa_end = c(3L, 5L)
), exon.temp, row.names = FALSE)
cna.bridge <- grin3d_prepare_cna_annotation_input(
  cna.temp, exon.temp, annotation.temp,
  mapping.file = tempfile(fileext = ".csv")
)
stopifnot(cna.bridge$clusters$residues == "1;2;3;4;5")
stopifnot(cna.bridge$clusters$n_residues == 5L)
stopifnot(cna.bridge$mapping$residue_segments == "1-3;3-5")

message("All GRIN3D annotation tests passed.")
