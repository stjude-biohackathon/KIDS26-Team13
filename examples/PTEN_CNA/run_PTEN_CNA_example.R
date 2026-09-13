# --------------------------------------------------------------------
# PTEN example for GRIN3D_exon_CNA_hotspot_prototype_scratch-code.R
# --------------------------------------------------------------------
# Transcript: ENST00000371953 (MANE Select)
# Protein:    ENSP00000361021
# UniProt:    P60484
#
# Expected input files
# --------------------
# PTEN_CNA_lesions.csv
#   ID, chrom, loc.start, loc.end, lsn.type
#
# PTEN_exon_to_protein_map.csv
#   gene_id, gene_name, transcript_id, protein_id, exon_id, exon_rank,
#   chrom, strand, cds_genomic_start, cds_genomic_end, aa_start, aa_end
#
# PTEN_alphafold_coordinates.csv
#   residue, x, y, z
#   plddt is required only when use.plddt.filter is TRUE.
#
# The exon map and AlphaFold coordinate file can be generated using
# GRIN3D_prepare_protein_inputs.R.
#
# Save the input files in an "input_files" subfolder inside analysis.dir.
# Results will be written to the "results" subfolder.
#
# ==============================================================================
# 1. DIRECTORIES AND INPUT FILES
# ==============================================================================

analysis.dir <- normalizePath(
  "C:/Users/aelsayed/Desktop/GRIN3D-biohackathon/GRIN3D/PTEN_CNA",
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
  "PTEN_CNA_lesions.csv"
)

exon.file <- file.path(
  input.dir,
  "PTEN_exon_to_protein_map.csv"
)

coordinate.file <- file.path(
  input.dir,
  "PTEN_alphafold_coordinates.csv"
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

# FALSE uses all available coordinates.
# TRUE applies the pLDDT threshold specified by min.confidence.
use.plddt.filter <- FALSE

print(
  grin3d_cna_input_spec(),
  row.names = FALSE
)

# ==============================================================================
# 2. OPTIONAL INPUT-COLUMN CHECK
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

require_grin3d_columns(
  coordinate.check,
  c(
    "residue", "x", "y", "z",
    if (use.plddt.filter) "plddt"
  ),
  basename(coordinate.file)
)


# ==============================================================================
# 3. ANALYSIS SETTINGS: change for each run
# ==============================================================================

protein.name <- "PTEN"
transcript.id <- "ENST00000371953"

# GRCh38 length of chromosome 10, where PTEN is located.
chromosome.length <- 133797422L

# CNA types to analyze separately.
# Use NULL to analyze every type present in the input file.
cna.types <- c("HOMDEL", "HETDEL", "GAIN", "AMP")

# Structural-confidence filter.
use.plddt.filter <- FALSE
minimum.plddt <- 70

# Simulation settings.
# Use 1,000 for development and at least 10,000 for stable comparisons.
n.coverage.simulations <- 1000L
n.structural.simulations <- 1000L
simulation.seed <- 20260907L
progress.interval <- 200L

# -----------------------------------
# 4. run the analysis
# -----------------------------------

pten.cna.results <- run_grin3d_exon_cna_hotspots(
  cna.file = cna.file,
  exon.file = exon.file,
  coordinate.file = coordinate.file,
  results.dir = results.dir,
  
  protein = protein.name,
  transcript = transcript.id,
  
  # Use NULL to analyze every CNA type observed in the input file.
  analysis.types = cna.types,
  
  # Add an ALL_CNA analysis combining the CNA types selected above.
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
  
  # Report conditional coverage results for CNAs affecting at least 50%, 75%,
  # 90% or 100% of the selected transcript's coding sequence.
  coverage.thresholds = c(0.50, 0.75, 0.90, 1.00),
  
  # Chromosome length for the selected genome build. This prevents simulated
  # genomic intervals from extending beyond the chromosome boundary.
  chromosome.length = chromosome.length,
  
  # Retain an exon as affected when the CNA overlaps any positive fraction of its
  # coding sequence. Increase this value only for breakpoint-sensitivity analyses.
  min.exon.coding.overlap = 0,
  
  # Allow partial CNAs affecting less than 100% of the coding sequence in the
  # localized 1D/3D analysis. Whole-coding-sequence events remain in gene-level
  # coverage results but are excluded from structural-hotspot testing.
  structural.max.coding.fraction = 0.999999,
  
  # When pLDDT filtering is enabled, retain coordinates with pLDDT greater than
  # or equal to minimum.plddt. This setting is ignored when
  # coordinate.columns$confidence is NULL.
  min.confidence = minimum.plddt,
  
  # Require at least 50% of an exon's encoded residues to have eligible structural
  # coordinates before using that exon in 3D analysis when pLDDT is specified.
  min.exon.mapped.fraction = 0.50,
  
  # Also require at least one eligible structural residue per exon.
  min.exon.mapped.residues = 1L,
  
  # Define exon proximity as the mean of the five smallest eligible C-alpha
  # distances between residues encoded by two exons.
  n.closest.cross.exon.distances = 5L,
  
  # When possible, exclude residue pairs fewer than 10 amino acids apart in the
  # linear sequence to reduce trivial proximity along the protein backbone.
  min.sequence.separation = 10L,
  
  # ---------------------------------------
  # Settings for relocating partial CNAs among eligible exon windows
  # ---------------------------------------
  
  # Controls preference for windows with a coding-residue span similar to the
  # observed CNA. Smaller values give stronger preference to similar spans.
  residue.span.temperature = 0.50,
  
  # Mixes span-based sampling with uniform sampling. A value of 0.25 gives
  # 75% weight to span similarity and 25% to equal placement probabilities.
  uniform.window.mixture = 0.25,
  
  # Flag an event when its sampling probabilities represent fewer than 1.25
  # effectively available exon windows.
  min.effective.windows = 1.25,
  
  # Require at least 5% total probability of placing the CNA outside its
  # observed exon window; otherwise flag the event as poorly permutable.
  min.alternative.probability = 0.05,
  
  # Test only candidate clusters affecting at least two subjects and two events.
  min.subjects = 2L,
  min.events = 2L,
  
  # Use 0.05 as the significance threshold for the corrected p-values.
  alpha = 0.05,
  
  # Use the simulation counts specified in the analysis-settings section.
  # Use 1,000 for development and at least 10,000 for stable comparisons.
  n.sim.coverage = n.coverage.simulations,
  n.sim.structural = n.structural.simulations,
  
  # Use 50% of structural simulations to calibrate size-specific p-values and
  # the remaining 50% to calculate protein-wide and joint 1D-3D corrections.
  calibration.fraction = 0.50,
  
  # Use the seed and progress interval specified in the settings section.
  random.seed = simulation.seed,
  progress.every = progress.interval
)

# ==============================================================================
# 5. INSPECT THE MAIN OUTPUTS
# ==============================================================================

print(pten.cna.results$structural_testability)
print(pten.cna.results$gene_coverage_results)
print(pten.cna.results$exon_recurrence_results)
print(pten.cna.results$clusters)

message("Results were saved to: ", normalizePath(results.dir, winslash = "/"))
