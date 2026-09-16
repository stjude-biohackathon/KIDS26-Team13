# TP53 example: add biological annotations after GRIN3D hotspot inference.
#
# Run from the repository root:
#   Rscript examples/TP53_mutations/run_TP53_annotation_example.R

analysis.dir <- normalizePath("examples/TP53_mutations", winslash = "/", mustWork = TRUE)
module.file <- normalizePath(
  "development-code/GRIN3D_annotate_protein_clusters.R",
  winslash = "/",
  mustWork = TRUE
)

source(module.file)

annotation.results <- annotate_grin3d_clusters(
  cluster.file = file.path(analysis.dir, "results", "mutation_hotspot_clusters.csv"),
  uniprot.accession = "P04637",
  output.dir = file.path(analysis.dir, "annotations"),
  protein = "TP53",
  sources = c("uniprot", "interpro", "pdbe", "chembl"),
  refresh = FALSE,
  offline = FALSE
)

write_cluster_annotation_dashboard(annotation.results)

print(annotation.results$summaries[, c(
  "cluster_id", "hotspot_class", "n_residues", "n_annotation_features",
  "n_annotated_cluster_residues", "feature_categories"
)])
