# --------------------------------------------------------------------
# CDKN2A p16 example for GRIN3D_exon_CNA_hotspots.R
# --------------------------------------------------------------------
# Transcript: ENST00000304494 (MANE Select)
# Protein:    ENSP00000307101
# UniProt:    P42771
# Length:     156 amino acids
#
# Expected input files
# --------------------
# CDKN2A_CNA_lesions.csv
#   ID, chrom, loc.start, loc.end, lsn.type
#
# CDKN2A_exon_to_protein_map.csv
#   gene_id, gene_name, transcript_id, protein_id, exon_id, exon_rank,
#   chrom, strand, cds_genomic_start, cds_genomic_end, aa_start, aa_end
#
# CDKN2A_alphafold_coordinates.csv
#   residue, x, y, z
#   plddt is required only when use.plddt.filter is TRUE.
#
# The exon map and AlphaFold coordinate file can be generated using
# GRIN3D_prepare_protein_inputs.R.
#
# Save the input files in an "input_files" subfolder inside analysis.dir.
# Results will be written to the "results" subfolder.

# ==============================================================================
# 1. DIRECTORIES AND INPUT FILES
# ==============================================================================

analysis.dir <- normalizePath(
  "C:/Users/aelsayed/Desktop/GRIN3D-biohackathon/GRIN3D/CDKN2A_CNA",
  winslash = "/",
  mustWork = TRUE
)

input.dir <- file.path(analysis.dir, "input_files")
results.dir <- file.path(analysis.dir, "results")

module.file <- file.path(
  analysis.dir,
  "GRIN3D_exon_CNA_hotspots.R"
)

cna.file <- file.path(
  input.dir,
  "CDKN2A_CNA_lesions.csv"
)

exon.file <- file.path(
  input.dir,
  "CDKN2A_exon_to_protein_map.csv"
)

coordinate.file <- file.path(
  input.dir,
  "CDKN2A_alphafold_coordinates.csv"
)

required.files <- c(
  module.file,
  cna.file,
  exon.file,
  coordinate.file
)

missing.files <- required.files[!file.exists(required.files)]

if (length(missing.files)) {
  stop(
    "Missing required file(s): ",
    paste(missing.files, collapse = "; ")
  )
}

dir.create(
  results.dir,
  recursive = TRUE,
  showWarnings = FALSE
)

source(module.file)

print(
  grin3d_cna_input_spec(),
  row.names = FALSE
)


# ==============================================================================
# 2. ANALYSIS SETTINGS: change for each run
# ==============================================================================

protein.name <- "CDKN2A"
transcript.id <- "ENST00000304494"

# GRCh38 length of chromosome 9, where CDKN2A is located.
chromosome.length <- 138394717L

# CNA types to analyze separately.
# Use NULL to analyze every CNA type present in the input file.
cna.types <- c("HOMDEL", "HETDEL", "GAIN", "AMP")

# FALSE uses all available coordinates.
# TRUE restricts the 3D analysis to residues meeting minimum.plddt.
use.plddt.filter <- FALSE
minimum.plddt <- 70

# Use 1,000 simulations for development and at least 10,000 for stable results.
n.coverage.simulations <- 1000L
n.structural.simulations <- 1000L
simulation.seed <- 20260907L
progress.interval <- 200L

# ==============================================================================
# 3. OPTIONAL INPUT-COLUMN CHECK
# ==============================================================================

cna.check <- utils::read.csv(
  cna.file,
  nrows = 1L,
  check.names = FALSE
)

exon.check <- utils::read.csv(
  exon.file,
  nrows = 1L,
  check.names = FALSE
)

coordinate.check <- utils::read.csv(
  coordinate.file,
  nrows = 1L,
  check.names = FALSE
)

require_grin3d_columns(
  cna.check,
  c("ID", "chrom", "loc.start", "loc.end", "lsn.type"),
  basename(cna.file)
)

require_grin3d_columns(
  exon.check,
  c(
    "gene_id", "gene_name", "transcript_id", "protein_id",
    "exon_id", "exon_rank", "chrom", "strand",
    "cds_genomic_start", "cds_genomic_end",
    "aa_start", "aa_end"
  ),
  basename(exon.file)
)

required.coordinate.columns <- c("residue", "x", "y", "z")

if (use.plddt.filter) {
  required.coordinate.columns <- c(
    required.coordinate.columns,
    "plddt"
  )
}

require_grin3d_columns(
  coordinate.check,
  required.coordinate.columns,
  basename(coordinate.file)
)


# ==============================================================================
# 4. RUN THE ANALYSIS
# ==============================================================================

cdkn2a.cna.results <- run_grin3d_exon_cna_hotspots(
  cna.file = cna.file,
  exon.file = exon.file,
  coordinate.file = coordinate.file,
  results.dir = results.dir,
  
  protein = protein.name,
  transcript = transcript.id,
  
  # Use NULL to analyze every CNA type observed in the input file.
  analysis.types = cna.types,
  
  # Add an ALL_CNA analysis combining the selected CNA types.
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
    confidence = if (use.plddt.filter) "plddt" else NULL
  ),
  
  # Report conditional results for CNAs affecting at least 50%, 75%, 90%,
  # or 100% of the selected transcript's coding sequence.
  coverage.thresholds = c(0.50, 0.75, 0.90, 1.00),
  
  # Prevent simulated genomic intervals from extending beyond chromosome 9.
  chromosome.length = chromosome.length,
  
  # Count any positive coding-exon overlap.
  min.exon.coding.overlap = 0,
  
  # Retain partial CNAs in localized tests. Whole-coding-sequence events remain
  # in gene-level summaries but are excluded from localized hotspot testing.
  structural.max.coding.fraction = 0.999999,
  
  # Apply this threshold only when use.plddt.filter is TRUE.
  min.confidence = minimum.plddt,
  
  # Require at least 50% of an exon's residues and at least one residue to have
  # eligible coordinates before including that exon in the 3D analysis.
  min.exon.mapped.fraction = 0.50,
  min.exon.mapped.residues = 1L,
  
  # Define exon proximity as the mean of the five smallest eligible C-alpha
  # distances between residues encoded by two exons.
  n.closest.cross.exon.distances = 5L,
  
  # When possible, exclude residue pairs separated by fewer than 10 residues
  # in sequence to reduce trivial protein-backbone proximity.
  min.sequence.separation = 10L,
  
  # ---------------------------------------
  # Settings for relocating partial CNAs
  # ---------------------------------------
  
  # Prefer exon windows with coding-residue spans similar to the observed CNA.
  residue.span.temperature = 0.50,
  
  # Assign 75% of placement weight by span similarity and 25% uniformly.
  uniform.window.mixture = 0.25,
  
  # Flag an event with fewer than 1.25 effectively available exon windows.
  min.effective.windows = 1.25,
  
  # Require at least 5% probability of placement outside the observed window.
  min.alternative.probability = 0.05,
  
  # Test clusters supported by at least two subjects and two events.
  min.subjects = 2L,
  min.events = 2L,
  
  # Significance threshold for corrected p-values.
  alpha = 0.05,
  
  # Use the simulation counts specified above.
  n.sim.coverage = n.coverage.simulations,
  n.sim.structural = n.structural.simulations,
  
  # Split structural simulations equally between size-specific calibration
  # and protein-wide/joint search correction.
  calibration.fraction = 0.50,
  
  random.seed = simulation.seed,
  progress.every = progress.interval
)


# ==============================================================================
# 5. INSPECT THE MAIN OUTPUTS
# ==============================================================================

print(cdkn2a.cna.results$structural_testability)
print(cdkn2a.cna.results$gene_coverage_results)
print(cdkn2a.cna.results$exon_recurrence_results)
print(cdkn2a.cna.results$clusters)

message(
  "Results were saved to: ",
  normalizePath(results.dir, winslash = "/")
)