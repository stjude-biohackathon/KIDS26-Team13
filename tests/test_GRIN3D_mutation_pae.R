source("development-code/GRIN3D_mutation_hotspots.R")

pae <- matrix(c(
  0, 2, 20,
  4, 0, 8,
  25, 6, 0
), nrow = 3L, byrow = TRUE)

clusters <- data.frame(
  cluster_id = c("supported", "uncertain", "singleton"),
  residues = c("1;2", "1;2;3", "3"),
  stringsAsFactors = FALSE
)

annotated <- annotate_clusters_with_pae(
  clusters, pae, cutoff = 10, minimum.reliable.fraction = 0.80
)
stopifnot(annotated$pae_n_pairs[[1L]] == 1L)
stopifnot(annotated$pae_median[[1L]] == 4)
stopifnot(isTRUE(annotated$pae_supported[[1L]]))
stopifnot(abs(annotated$pae_reliable_pair_fraction[[2L]] - 2 / 3) < 1e-12)
stopifnot(identical(annotated$pae_supported[[2L]], FALSE))
stopifnot(is.na(annotated$pae_supported[[3L]]))
stopifnot(annotated$pae_interpretation[[3L]] == "singleton_no_pairwise_pae")

if (requireNamespace("jsonlite", quietly = TRUE)) {
  path <- tempfile(fileext = ".json")
  jsonlite::write_json(
    list(list(predicted_aligned_error = unname(split(pae, row(pae))))),
    path, auto_unbox = TRUE
  )
  imported <- read_alphafold_pae(path, expected.residues = 1:3)
  stopifnot(identical(dim(imported), c(3L, 3L)))
  stopifnot(unname(imported[3L, 1L]) == 25)
  unlink(path)
}

message("All GRIN3D mutation PAE tests passed.")
