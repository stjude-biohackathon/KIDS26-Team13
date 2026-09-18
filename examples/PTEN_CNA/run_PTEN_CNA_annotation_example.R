# PTEN CNA example: annotate exon-hotspot candidates after GRIN3D testing.
#
# Run from the repository root:
#   Rscript examples/PTEN_CNA/run_PTEN_CNA_annotation_example.R
#
# CNA candidates are exon sets, whereas the annotation sources use UniProt
# residue coordinates. The bridge below assigns each candidate the union of
# residues encoded by its member exons. It does not change cluster membership,
# p-values, or significance.

analysis.dir <- normalizePath("examples/PTEN_CNA", winslash = "/", mustWork = TRUE)
module.file <- normalizePath(
  "development-code/GRIN3D_annotate_protein_clusters.R",
  winslash = "/",
  mustWork = TRUE
)

source(module.file)

annotation.dir <- file.path(analysis.dir, "annotations")
annotation.input <- file.path(annotation.dir, "CNA_annotation_clusters.csv")

grin3d_prepare_cna_annotation_input(
  cna.cluster.file = file.path(
    analysis.dir, "results", "CNA_exon_hotspot_clusters.csv"
  ),
  exon.map.file = file.path(
    analysis.dir, "input_files", "PTEN_exon_to_protein_map.csv"
  ),
  output.file = annotation.input
)

annotation.results <- annotate_grin3d_clusters(
  cluster.file = annotation.input,
  uniprot.accession = "P60484",
  output.dir = annotation.dir,
  protein = "PTEN",
  sources = c("uniprot", "interpro", "pdbe", "chembl"),
  refresh = FALSE,
  offline = FALSE
)

annotation.results$position_label <- "Residues encoded by candidate exons"
write_cluster_annotation_dashboard(annotation.results)

significant <- annotation.results$summaries$significant_any %in% TRUE
print(annotation.results$summaries[significant, c(
  "cluster_id", "hotspot_class", "n_residues", "n_annotation_features",
  "n_annotated_cluster_residues", "feature_categories"
)])
