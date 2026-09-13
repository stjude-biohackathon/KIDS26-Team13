
# TP53 example for GRIN3D_mutation_hotspot_prototype_scratch-code.R
#
# Expected input files
# --------------------
# tp53_start_end.csv
#   ID, id, start, end
#
# tp53_coordinates.csv
#   residue, x, y, z
#   plddt and opportunity are optional.
#   this file can be generated for selected proteins using GRIN3D_exon_to_protein_download_alphafold.R script

# ==============================================================================
# 1. SETTINGS TO CHANGE FOR EACH PROTEIN
# ==============================================================================

analysis.dir <- "C:/Users/aelsayed/Desktop/GRIN3D-biohackathon/GRIN3D/TP53_mutations"
protein.name <- "TP53"

module.name <- "GRIN3D_mutation_hotspots.R"
lesion.name <- "tp53_start_end.csv"
coordinate.name <- "tp53_coordinates.csv"

# Column names in the input files.
lesion.columns <- list(
  subject = "ID",
  event = "id",       # Set to NULL if there is no unique event-ID column.
  start = "start",
  end = "end"
)

use.plddt.filter <- FALSE
use.opportunity.weights <- FALSE

coordinate.columns <- list(
  residue = "residue",
  x = "x",
  y = "y",
  z = "z",
  confidence = if (use.plddt.filter) "plddt" else NULL,
  opportunity = if (use.opportunity.weights) "opportunity" else NULL
)

n.simulations <- 10000L

# ==============================================================================
# 2. INPUT FILES AND ANALYSIS DIRECTORY
# ==============================================================================

analysis.dir <- normalizePath(
  analysis.dir,
  winslash = "/",
  mustWork = TRUE
)

input.dir <- file.path(analysis.dir, "input_files")
results.dir <- file.path(analysis.dir, "results")

module.file <- file.path(analysis.dir, module.name)
lesion.file <- file.path(input.dir, lesion.name)
coordinate.file <- file.path(input.dir, coordinate.name)

required.files <- c(module.file, lesion.file, coordinate.file)
missing.files <- required.files[!file.exists(required.files)]

if (length(missing.files)) {
  stop("Missing required file(s): ", paste(missing.files, collapse = "; "))
}

source(module.file)

# ==============================================================================
# 3. RUN THE MUTATION-HOTSPOT ANALYSIS
# ==============================================================================

mutation.results <- run_grin3d_mutation_hotspots(
  lesion.file = lesion.file,
  coordinate.file = coordinate.file,
  results.dir = results.dir,
  protein = protein.name,
  
  lesion.columns = lesion.columns,
  coordinate.columns = coordinate.columns,
  
  min.subjects = 2L,
  min.events = 2L,
  
  # Applied only when the confidence column is set to "plddt".
  min.confidence = 70,
  require.complete.event.mapping = TRUE,
  
  # Avoid overlap between a subject's simulated events when possible.
  avoid.within.subject.overlap = TRUE,
  
  alpha = 0.05,
  n.sim = n.simulations,
  calibration.fraction = 0.50,
  random.seed = 20260828L,
  progress.every = 200L
)

# ==============================================================================
# 4. INSPECT THE MAIN OUTPUTS
# ==============================================================================

print(mutation.results$clusters)
print(mutation.results$size_summary)
print(mutation.results$event_mapping_summary)

message("Results were saved to: ", normalizePath(results.dir, winslash = "/"))
