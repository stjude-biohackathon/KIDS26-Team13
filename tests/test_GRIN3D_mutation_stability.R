source("development-code/GRIN3D_mutation_stability.R")

reference <- data.frame(
  cluster_id = c("r1", "r2"), residues = c("1;2;3", "10;11"),
  n_residues = c(3L, 2L), tree_sources = c("1D", "3D"),
  significant_any = c(TRUE, TRUE), p_any_joint = c(0.01, 0.02),
  stringsAsFactors = FALSE
)
comparison <- data.frame(
  cluster_id = c("c1", "c2", "c3"), residues = c("1;2;3", "10;11;12", "50"),
  n_residues = c(3L, 3L, 1L), tree_sources = c("1D", "3D", "1D"),
  significant_any = c(TRUE, FALSE, FALSE), p_any_joint = c(0.02, 0.08, 1),
  stringsAsFactors = FALSE
)

similarity <- grin3d_cluster_similarity_matrix(reference, comparison)
stopifnot(similarity$jaccard["r1", "c1"] == 1)
stopifnot(abs(similarity$jaccard["r2", "c2"] - 2 / 3) < 1e-12)
stopifnot(similarity$containment["r2", "c2"] == 1)

matches <- grin3d_match_clusters(reference, comparison)
stopifnot(nrow(matches) == 2L)
stopifnot(matches$comparison_cluster_id[matches$reference_cluster_id == "r1"] == "c1")
stopifnot(
  matches$membership_interpretation[matches$reference_cluster_id == "r2"] ==
    "stable_nested_core"
)

message("All GRIN3D mutation-stability tests passed.")
