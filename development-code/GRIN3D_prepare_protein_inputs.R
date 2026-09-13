
# GRIN3D: create exon-to-protein map and retreive AlphaFold C-alpha input files.
#
# Set the gene, transcript, and output location in the SETTINGS section below.
# The script compares complete Ensembl and AlphaFold protein sequences before
# writing either CSV. A matching gene name or protein length alone is not enough.
# The exon output has one row per coding exon; split residues appear twice.
#
# Install the required Bioconductor packages once if needed:
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install(c(
#   "AnnotationHub", "AnnotationFilter", "ensembldb", "GenomicRanges",
#   "IRanges", "S4Vectors"
# ))
# install.packages(c("jsonlite", "bio3d"))

# ==============================================================================
# SETTINGS: change these for another protein
# ==============================================================================

gene.name <- "PTEN"
transcript.id <- "ENST00000371953"
expected.protein.id <- "ENSP00000361021"  # Set to NULL if unknown.
uniprot.accession <- "P60484"             # AlphaFold accession for this isoform.
expected.protein.length <- 403L            # Set to NULL to skip this check.
ensembl.release <- 110L
genome.build <- "GRCh38"
output.dir <- getwd()                       # Folder for both output CSV files.
pdb.dir <- file.path(output.dir, "AlphaFold_PDB")

# The script uses the sequence supplied by the selected EnsDb annotation.
# Check that the genome build matches the coordinates in your CNA input.

required.packages <- c(
  "AnnotationHub", "AnnotationFilter", "ensembldb", "GenomicRanges",
  "IRanges", "S4Vectors", "jsonlite", "bio3d"
)

missing.packages <- required.packages[
  !vapply(required.packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing.packages) > 0L) {
  stop(
    "Install the following Bioconductor package(s) before running this script: ",
    paste(missing.packages, collapse = ", ")
  )
}

strip_ensembl_version <- function(x) {
  sub("\\.[0-9]+$", "", as.character(x))
}

# Retrieve one versioned human EnsDb resource from AnnotationHub.
get_human_ensdb <- function(ensembl.release = 110L, genome = "GRCh38") {
  hub <- AnnotationHub::AnnotationHub()
  hits <- AnnotationHub::query(
    hub,
    c("EnsDb", "Homo sapiens", as.character(ensembl.release))
  )

  if (length(hits) == 0L) {
    stop("No human EnsDb resource was found for Ensembl release ", ensembl.release)
  }

  metadata <- as.data.frame(S4Vectors::mcols(hits))
  expected.title <- paste0(
    "Ensembl ", ensembl.release, " EnsDb for Homo sapiens"
  )
  keep <- tolower(metadata$title) == tolower(expected.title) &
    metadata$genome == genome

  if (sum(keep) != 1L) {
    available <- paste(
      paste0(metadata$title, " [", metadata$genome, "]"),
      collapse = "; "
    )
    stop(
      "Expected exactly one matching EnsDb resource but found ", sum(keep),
      ". Available candidates: ", available
    )
  }

  hits[[which(keep)]]
}

# Keep the first n nucleotides of a CDS in transcript (5' to 3') order.
# This removes the stop codon when cdsBy() includes it. Exon metadata are kept.
keep_translated_cds <- function(cds.parts, translated.nt) {
  if (!inherits(cds.parts, "GRanges")) {
    stop("cds.parts must be a GRanges object")
  }
  if (length(cds.parts) == 0L || translated.nt < 1L) {
    stop("No translated CDS nucleotides were supplied")
  }

  exon.rank <- S4Vectors::mcols(cds.parts)$exon_rank
  if (is.null(exon.rank) || anyNA(exon.rank)) {
    stop("The CDS annotation does not contain complete exon_rank values")
  }

  cds.parts <- cds.parts[order(exon.rank)]
  remaining <- as.integer(translated.nt)
  kept <- vector("list", length(cds.parts))
  n.kept <- 0L

  for (i in seq_along(cds.parts)) {
    if (remaining == 0L) {
      break
    }

    current <- cds.parts[i]
    take <- min(IRanges::width(current), remaining)

    # Transcript order proceeds from low to high genomic coordinates on the
    # plus strand and from high to low coordinates on the minus strand.
    current.strand <- as.character(GenomicRanges::strand(current))
    if (current.strand == "+") {
      IRanges::end(current) <- IRanges::start(current) + take - 1L
    } else if (current.strand == "-") {
      IRanges::start(current) <- IRanges::end(current) - take + 1L
    } else {
      stop("CDS strand must be '+' or '-'")
    }

    n.kept <- n.kept + 1L
    kept[[n.kept]] <- current
    remaining <- remaining - take
  }

  if (remaining > 0L) {
    stop(
      "The annotated CDS is shorter than the translated protein by ",
      remaining, " nucleotide(s)"
    )
  }

  do.call(c, kept[seq_len(n.kept)])
}

# Map every coding exon of one transcript to its inclusive amino-acid range.
map_transcript_exons_to_protein <- function(
    edb,
    transcript.id,
    expected.protein.id = NULL,
    require.cds.ok = TRUE) {

  transcript.id <- strip_ensembl_version(transcript.id)
  if (!is.null(expected.protein.id)) {
    expected.protein.id <- strip_ensembl_version(expected.protein.id)
  }

  if (!ensembldb::hasProteinData(edb)) {
    stop("The supplied EnsDb does not contain protein annotation data")
  }
  tx.filter <- AnnotationFilter::TxIdFilter(transcript.id)

  # Obtain the transcript-specific protein sequence and identifier.
  protein.annotation <- ensembldb::proteins(
    edb,
    filter = tx.filter,
    columns = c(
      "gene_id", "gene_name", "tx_id", "protein_id", "protein_sequence"
    )
  )
  protein.data <- as.data.frame(protein.annotation)

  if (nrow(protein.data) == 0L) {
    stop("No encoded protein was found for transcript ", transcript.id)
  }

  protein.data$protein_id <- strip_ensembl_version(protein.data$protein_id)
  if (!is.null(expected.protein.id)) {
    protein.data <- protein.data[
      protein.data$protein_id == expected.protein.id,
      ,
      drop = FALSE
    ]
  }

  protein.ids <- unique(protein.data$protein_id)
  if (length(protein.ids) != 1L) {
    stop(
      "Transcript ", transcript.id,
      " did not resolve to exactly one Ensembl protein ID"
    )
  }
  protein.id <- protein.ids[[1L]]

  protein.sequence <- unique(as.character(
    protein.data$protein_sequence[protein.data$protein_id == protein.id]
  ))
  protein.sequence <- unique(sub("\\*$", "", protein.sequence))
  protein.sequence <- protein.sequence[nzchar(protein.sequence)]
  if (length(protein.sequence) != 1L) {
    stop("A unique protein sequence was not available for ", protein.id)
  }
  protein.length <- nchar(protein.sequence)
  translated.nt <- 3L * protein.length

  # Retrieve CDS pieces grouped by transcript. Unlike complete exons, these
  # ranges exclude UTR sequence and retain exon_id and exon_rank metadata.
  cds.by.tx <- ensembldb::cdsBy(
    edb,
    by = "tx",
    filter = tx.filter
  )
  if (length(cds.by.tx) != 1L) {
    stop(
      "Expected one CDS record for ", transcript.id,
      " but found ", length(cds.by.tx)
    )
  }

  cds.parts.all <- cds.by.tx[[1L]]
  cds.total.nt <- sum(IRanges::width(cds.parts.all))
  if (cds.total.nt < translated.nt) {
    stop(
      "CDS length (", cds.total.nt,
      ") is shorter than 3 x protein length (", translated.nt, ")"
    )
  }

  # Retain exactly the nucleotides translated into amino acids. This safely
  # removes a terminal stop codon if the EnsDb CDS representation includes it.
  cds.parts <- keep_translated_cds(cds.parts.all, translated.nt)

  # genomeToProtein() uses CDS phase, transcript strand and protein annotation
  # to map each genomic CDS segment to protein-relative positions.
  # Restrict the mapping resources to the requested transcript. Without these
  # prefiltered objects, genomeToProtein() also evaluates overlapping isoforms
  # and can emit warnings about unrelated incomplete or noncoding transcripts.
  selected.exons <- ensembldb::exonsBy(
    edb,
    by = "tx",
    filter = tx.filter
  )
  selected.transcripts <- ensembldb::transcripts(
    edb,
    filter = tx.filter
  )
  mapped <- ensembldb::genomeToProtein(
    cds.parts,
    edb,
    proteins = protein.annotation,
    exons = selected.exons,
    transcripts = selected.transcripts
  )

  rows <- lapply(seq_along(cds.parts), function(i) {
    one <- mapped[[i]]
    if (length(one) == 0L) {
      stop("CDS part ", i, " could not be mapped to a protein")
    }

    metadata <- as.data.frame(S4Vectors::mcols(one))
    mapped.protein <- strip_ensembl_version(names(one))
    mapped.tx <- strip_ensembl_version(metadata$tx_id)
    keep <- mapped.protein == protein.id & mapped.tx == transcript.id &
      IRanges::start(one) > 0L

    if (!any(keep)) {
      stop(
        "CDS part ", i, " did not map to transcript ", transcript.id,
        " and protein ", protein.id
      )
    }

    one <- one[keep]
    metadata <- as.data.frame(S4Vectors::mcols(one))
    cds.ok.values <- metadata$cds_ok
    cds.ok <- if (all(is.na(cds.ok.values))) {
      NA
    } else {
      all(cds.ok.values %in% TRUE)
    }

    data.frame(
      gene_id = unique(protein.data$gene_id)[1L],
      gene_name = unique(protein.data$gene_name)[1L],
      transcript_id = transcript.id,
      protein_id = protein.id,
      exon_id = as.character(S4Vectors::mcols(cds.parts)$exon_id[i]),
      exon_rank = as.integer(S4Vectors::mcols(cds.parts)$exon_rank[i]),
      chrom = as.character(GenomicRanges::seqnames(cds.parts)[i]),
      strand = as.character(GenomicRanges::strand(cds.parts)[i]),
      cds_genomic_start = IRanges::start(cds.parts)[i],
      cds_genomic_end = IRanges::end(cds.parts)[i],
      cds_length_nt = IRanges::width(cds.parts)[i],
      aa_start = min(IRanges::start(one)),
      aa_end = max(IRanges::end(one)),
      cds_ok = cds.ok,
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, rows)
  result <- result[order(result$exon_rank), , drop = FALSE]
  rownames(result) <- NULL
  result$n_aa_touched <- result$aa_end - result$aa_start + 1L

  # Identify residues encoded across an exon junction. Such residues correctly
  # appear in both exon residue ranges and must not be double counted.
  result$shared_with_previous_exon <- NA_character_
  result$shared_with_next_exon <- NA_character_
  if (nrow(result) > 1L) {
    for (i in 2:nrow(result)) {
      overlap.start <- max(result$aa_start[i - 1L], result$aa_start[i])
      overlap.end <- min(result$aa_end[i - 1L], result$aa_end[i])
      if (overlap.start <= overlap.end) {
        shared <- paste(seq.int(overlap.start, overlap.end), collapse = ";")
        result$shared_with_next_exon[i - 1L] <- shared
        result$shared_with_previous_exon[i] <- shared
      }
    }
  }

  mapped.residues <- unique(unlist(Map(
    seq.int,
    result$aa_start,
    result$aa_end
  )))
  expected.residues <- seq_len(protein.length)
  complete.mapping <- identical(sort(mapped.residues), expected.residues)

  if (!complete.mapping) {
    stop("The exon mappings do not cover residues 1 through ", protein.length)
  }
  if (require.cds.ok && any(is.na(result$cds_ok) | !result$cds_ok)) {
    stop("At least one mapped CDS segment has cds_ok != TRUE")
  }

  attr(result, "protein_length") <- protein.length
  attr(result, "protein_sequence") <- protein.sequence
  attr(result, "total_cds_nt_including_stop_if_present") <- cds.total.nt
  attr(result, "translated_cds_nt") <- translated.nt
  result
}

# Download the PDB linked to the requested UniProt accession.
download_alphafold_pdb <- function(uniprot.accession, pdb.dir) {
  api.url <- paste0(
    "https://alphafold.ebi.ac.uk/api/prediction/",
    utils::URLencode(uniprot.accession, reserved = TRUE)
  )
  metadata <- tryCatch(
    jsonlite::fromJSON(api.url),
    error = function(e) stop("AlphaFold lookup failed: ", conditionMessage(e))
  )
  if (!is.data.frame(metadata) || nrow(metadata) < 1L ||
      !"pdbUrl" %in% names(metadata) ||
      is.na(metadata$pdbUrl[1L]) || !nzchar(metadata$pdbUrl[1L])) {
    stop("No AlphaFold PDB was found for ", uniprot.accession)
  }
  pdb.url <- metadata$pdbUrl[1L]
  dir.create(pdb.dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(pdb.dir)) stop("Could not create PDB directory: ", pdb.dir)
  pdb.file <- file.path(pdb.dir, basename(pdb.url))

  # Retain a previously downloaded model for reproducible reruns.
  if (!file.exists(pdb.file) || file.info(pdb.file)$size == 0) {
    temporary.file <- tempfile("alphafold_", fileext = ".pdb")
    on.exit(unlink(temporary.file), add = TRUE)
    status <- tryCatch(
      utils::download.file(pdb.url, temporary.file, mode = "wb", quiet = TRUE),
      error = function(e) stop("PDB download failed: ", conditionMessage(e))
    )
    if (!identical(status, 0L) || !file.exists(temporary.file) ||
        file.info(temporary.file)$size == 0) {
      stop("Incomplete PDB download from: ", pdb.url)
    }
    replace.empty <- file.exists(pdb.file) && file.info(pdb.file)$size == 0
    if (!file.copy(temporary.file, pdb.file, overwrite = replace.empty)) {
      stop("Cannot save downloaded PDB: ", pdb.file)
    }
  }
  pdb.file
}

# Keep each residue's C-alpha coordinates and pLDDT (PDB B-factor field).
extract_ca_coordinates <- function(pdb.file) {
  atoms <- bio3d::read.pdb(pdb.file)$atom
  ca <- atoms[trimws(toupper(atoms$elety)) == "CA", , drop = FALSE]
  if (!nrow(ca)) stop("No C-alpha atoms found in: ", pdb.file)
  result <- data.frame(
    residue = as.integer(ca$resno),
    amino_acid = as.character(ca$resid),
    chain = as.character(ca$chain),
    x = as.numeric(ca$x), y = as.numeric(ca$y), z = as.numeric(ca$z),
    plddt = as.numeric(ca$b), stringsAsFactors = FALSE
  )
  if (anyNA(result$residue) ||
      anyNA(result[c("x", "y", "z", "plddt")])) {
    stop("PDB contains missing C-alpha positions or pLDDT values")
  }
  if (length(unique(result$chain)) != 1L ||
      anyDuplicated(result$residue)) {
    stop("Expected one chain with exactly one C-alpha atom per residue")
  }
  result <- result[order(result$residue), , drop = FALSE]
  rownames(result) <- NULL
  result
}

# Fail when the structure is numbered differently or encodes another isoform.
check_matching_sequences <- function(exon.map, coordinates) {
  ensembl.sequence <- toupper(attr(exon.map, "protein_sequence"))
  if (length(ensembl.sequence) != 1L || is.na(ensembl.sequence) ||
      !nzchar(ensembl.sequence)) {
    stop("An Ensembl protein sequence is required for comparison")
  }
  expected.positions <- seq_len(nchar(ensembl.sequence))
  if (!identical(coordinates$residue, expected.positions)) {
    stop(
      "The AlphaFold PDB does not contain precisely residues 1-",
      length(expected.positions), " of the Ensembl protein. Check the isoform."
    )
  }

  # PDB amino-acid names use three letters; the Ensembl sequence uses one.
  residue.codes <- c(
    ALA = "A", ARG = "R", ASN = "N", ASP = "D", CYS = "C",
    GLN = "Q", GLU = "E", GLY = "G", HIS = "H", ILE = "I",
    LEU = "L", LYS = "K", MET = "M", PHE = "F", PRO = "P",
    SER = "S", THR = "T", TRP = "W", TYR = "Y", VAL = "V"
  )
  pdb.residues <- toupper(trimws(coordinates$amino_acid))
  pdb.letters <- unname(residue.codes[pdb.residues])
  if (anyNA(pdb.letters)) {
    stop(
      "Unknown PDB amino acid at residue(s): ",
      paste(head(coordinates$residue[is.na(pdb.letters)], 10L), collapse = ", ")
    )
  }
  ensembl.letters <- strsplit(ensembl.sequence, "", fixed = TRUE)[[1L]]
  mismatch <- which(pdb.letters != ensembl.letters)
  if (length(mismatch)) {
    stop(
      "Ensembl and AlphaFold protein sequences differ at ",
      length(mismatch), " residue(s); first mismatch at position ",
      mismatch[[1L]], " (Ensembl ", ensembl.letters[mismatch[[1L]]],
      ", AlphaFold ", pdb.letters[mismatch[[1L]]], "). ",
      "Select a matching UniProt isoform before running GRIN3D."
    )
  }
  invisible(TRUE)
}

# ==============================================================================
# RUN AND WRITE BOTH MATCHED CSV FILES
# ==============================================================================

if (!is.character(gene.name) || length(gene.name) != 1L ||
    is.na(gene.name) || !nzchar(gene.name)) {
  stop("Set gene.name to one nonempty gene symbol")
}
if (!is.character(transcript.id) || length(transcript.id) != 1L ||
    is.na(transcript.id) || !nzchar(transcript.id)) {
  stop("Set transcript.id to one Ensembl transcript ID")
}
if (!is.null(expected.protein.length) &&
    (length(expected.protein.length) != 1L || is.na(expected.protein.length) ||
     expected.protein.length < 1L)) {
  stop("expected.protein.length must be a positive number or NULL")
}
if (!dir.exists(output.dir)) {
  stop("Output directory does not exist: ", output.dir)
}
if (!is.character(uniprot.accession) || length(uniprot.accession) != 1L ||
    is.na(uniprot.accession) || !nzchar(uniprot.accession)) {
  stop("Set uniprot.accession to the matching UniProt protein accession")
}

edb <- get_human_ensdb(
  ensembl.release = ensembl.release,
  genome = genome.build
)

exon.protein.map <- map_transcript_exons_to_protein(
  edb = edb,
  transcript.id = transcript.id,
  expected.protein.id = expected.protein.id
)

mapped.genes <- unique(exon.protein.map$gene_name)
if (!identical(mapped.genes, gene.name)) {
  stop(
    "The transcript maps to gene ", paste(mapped.genes, collapse = ", "),
    ", not the requested gene ", gene.name
  )
}

# Check the expected protein length only when one was supplied above.
protein.length <- attr(exon.protein.map, "protein_length")
if (!is.null(expected.protein.length) &&
    protein.length != expected.protein.length) {
  stop(
    "Expected ", expected.protein.length, " amino acids, but ",
    transcript.id, " encodes ", protein.length
  )
}

pdb.file <- download_alphafold_pdb(uniprot.accession, pdb.dir)
ca.coordinates <- extract_ca_coordinates(pdb.file)
check_matching_sequences(exon.protein.map, ca.coordinates)

exon.file <- file.path(output.dir, paste0(gene.name, "_exon_to_protein_map.csv"))
coordinate.file <- file.path(
  output.dir, paste0(gene.name, "_alphafold_coordinates.csv")
)
utils::write.csv(exon.protein.map, exon.file, row.names = FALSE, na = "")
utils::write.csv(ca.coordinates, coordinate.file, row.names = FALSE, na = "")

message(
  "Exact protein sequence match: ", gene.name, " ", transcript.id,
  " / ", unique(exon.protein.map$protein_id), " / ", uniprot.accession,
  " (", protein.length, " residues)."
)
message("Exon map: ", exon.file)
message("C-alpha coordinates: ", coordinate.file)
message("Downloaded model: ", pdb.file)

# The CNA prototype uses exon_id, exon_rank, genomic CDS coordinates,
# aa_start, and aa_end. Adjacent exons can share a split amino acid;
# shared_with_previous_exon and shared_with_next_exon mark these residues.
