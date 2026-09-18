# Build SNV and indel hotspot examples from the public T-ALL mutation set.
# SNVs and indels are analysed separately (separate GRIN3D runs), each placed on
# the protein's AlphaFold model. Run from the repository root:
#   Rscript visualization/prepare_mutation_examples.R            # all genes
#   Rscript visualization/prepare_mutation_examples.R PTEN SNV   # one gene / one type
#   Rscript visualization/prepare_mutation_examples.R PTEN SNV 10000   # with more simulations

genes <- list(
  PTEN = list(transcript = "NM_000314",
              pdb = "examples/PTEN_CNA/input_files/AF-P60484-F1-model_v6.pdb"),
  TP53 = list(transcript = "NM_000546",
              pdb = "examples/TP53_mutations/input_files/AF-P04637-F1-model_v6.pdb"),
  FBXW7 = list(transcript = "NM_033632",
               pdb = "examples/FBXW7_mutations/input_files/AF-Q969H0-F1-model_v6.pdb"),
  JAK3 = list(transcript = "NM_000215",
              pdb = "examples/JAK3_mutations/input_files/AF-P52333-F1-model_v6.pdb"),
  IL7R = list(transcript = "NM_002185",
              pdb = "examples/IL7R_mutations/input_files/AF-P16871-F1-model_v6.pdb"),
  LEF1 = list(transcript = "NM_016269",
              pdb = "examples/LEF1_CNA/input_files/AF-Q9UJU2-F1-model_v6.pdb")
)
mutation.xlsx <- "datasets/T_ALL_public_data/mutations_dataset.xlsx"
n.simulations <- 1000L  # for speed; raise to 10000 for final results (joint p floor is ~0.1 at 1000)

source("development-code/GRIN3D_mutation_hotspots.R")

args <- commandArgs(trailingOnly = TRUE)
run.genes <- if (length(args) >= 1) args[1] else names(genes)
run.types <- if (length(args) >= 2) args[2] else c("SNV", "INDEL")
if (length(args) >= 3) n.simulations <- as.integer(args[3])

# C-alpha coordinates and pLDDT straight from the AlphaFold PDB file.
read_ca <- function(pdb.file) {
  p <- readLines(pdb.file, warn = FALSE)
  ca <- p[startsWith(p, "ATOM") & substr(p, 13, 16) == " CA "]
  data.frame(residue = as.integer(substr(ca, 23, 26)),
             x = as.numeric(substr(ca, 31, 38)),
             y = as.numeric(substr(ca, 39, 46)),
             z = as.numeric(substr(ca, 47, 54)),
             plddt = as.numeric(substr(ca, 61, 66)))
}

# Protein start/end from aa_change: p.N31H, p.P30fs, p.132_133del, p.G251_D252delinsAVX.
aa_position <- function(aa) {
  ok <- grepl("^p\\.[A-Z*]?\\d+", aa)
  start <- as.integer(ifelse(ok, sub("^p\\.[A-Z*]?(\\d+).*$", "\\1", aa), NA))
  ranged <- grepl("^p\\.[A-Z*]?\\d+_[A-Z*]?\\d+", aa)
  end <- ifelse(ranged, as.integer(sub("^p\\.[A-Z*]?\\d+_[A-Z*]?(\\d+).*$", "\\1", aa)), start)
  data.frame(start = start, end = end)
}

stopifnot(
  identical(aa_position(c("p.N31H", "p.P30fs", "p.132_133del", "p.G251_D252delinsAVX", "."))$start,
            c(31L, 30L, 132L, 251L, NA)),
  identical(aa_position(c("p.N31H", "p.132_133del", "p.G251_D252delinsAVX"))$end, c(31L, 133L, 252L))
)

muts <- as.data.frame(readxl::read_excel(mutation.xlsx))

for (gene in run.genes) {
  g <- genes[[gene]]
  d <- muts[muts$gene == gene & muts$NCBI.accession.number == g$transcript, ]
  d <- cbind(d, aa_position(d$aa_change))
  d$type <- ifelse(d$var_type == "SNV", "SNV", "INDEL")
  d <- d[!is.na(d$start), ]
  d <- d[!duplicated(d[, c("sample", "position", "ref", "alt")]), ]
  example.dir <- file.path("examples", paste0(gene, "_mutations"))
  dir.create(file.path(example.dir, "input_files"), recursive = TRUE, showWarnings = FALSE)
  coordinate.file <- file.path(example.dir, "input_files", paste0(gene, "_alphafold_ca_plddt.csv"))
  write.csv(read_ca(g$pdb), coordinate.file, row.names = FALSE)

  for (type in intersect(run.types, unique(d$type))) {
    x <- d[d$type == type, ]
    lesions <- data.frame(id = seq_len(nrow(x)), ID = x$sample, start = x$start, end = x$end)
    lesion.file <- file.path(example.dir, "input_files", sprintf("%s_%s_lesions.csv", gene, type))
    write.csv(lesions, lesion.file, row.names = FALSE)
    message(sprintf("%s %s: %d events from %d subjects", gene, type,
                    nrow(lesions), length(unique(lesions$ID))))

    run_grin3d_mutation_hotspots(
      lesion.file = lesion.file,
      coordinate.file = coordinate.file,
      results.dir = file.path(example.dir, paste0("results_", type)),
      protein = gene,
      lesion.columns = list(subject = "ID", event = "id", start = "start", end = "end"),
      coordinate.columns = list(residue = "residue", x = "x", y = "y", z = "z",
                                confidence = "plddt", opportunity = NULL),
      min.subjects = 2L,
      min.events = 2L,
      min.confidence = 70,
      require.complete.event.mapping = TRUE,
      avoid.within.subject.overlap = TRUE,
      alpha = 0.05,
      n.sim = n.simulations,
      calibration.fraction = 0.50,
      random.seed = 20260828L,
      progress.every = 1000L
    )
  }
}
