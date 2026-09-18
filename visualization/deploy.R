# Deploy the viewer to shinyapps.io. Run from the repository root:
#   Rscript visualization/deploy.R
# Set the account once (get the values from shinyapps.io > Account > Tokens):
#   rsconnect::setAccountInfo(name = "<account>", token = "<token>", secret = "<secret>")
#
# Only the files the running app reads are uploaded: the viewer tables, the
# AlphaFold models, and the three annotation files per protein. The rest of
# data/annotations (per-residue tables, cached API responses) is build output
# the app never opens, and is ~40 MB.

app_dir <- "visualization"
account <- Sys.getenv("SHINYAPPS_ACCOUNT", unset = NA)
proteins <- read.csv(file.path(app_dir, "data", "proteins.csv"))

# Copy each protein's PDB into the app folder; app.R prefers this copy when present.
structures <- file.path(app_dir, "data", "structures")
dir.create(structures, showWarnings = FALSE)
for (f in proteins$pdb_file) file.copy(f, file.path(structures, basename(f)), overwrite = TRUE)

annotation_files <- function(protein) {
  f <- file.path("data", "annotations", protein,
                 c("protein_features.csv", "cluster_feature_overlaps.csv", "chembl_drug_mechanisms.csv"))
  f[file.exists(file.path(app_dir, f))]
}

app_files <- c(
  "app.R", "annotation_track.R", "grin3d_viewer_input.R", "cluster_export.R",
  file.path("data", c("proteins.csv", sprintf("%s_viewer_%s.csv",
    rep(proteins$protein, each = 3), c("clusters", "residues", "members")))),
  file.path("data", "structures", basename(proteins$pdb_file)),
  unlist(lapply(proteins$protein, annotation_files)))

missing <- app_files[!file.exists(file.path(app_dir, app_files))]
if (length(missing)) stop("missing: ", paste(missing, collapse = ", "))
message(sprintf("%d files, %.1f MB", length(app_files),
                sum(file.size(file.path(app_dir, app_files))) / 1e6))

rsconnect::deployApp(
  appDir = app_dir, appFiles = app_files,
  appName = "grin3d", appTitle = "GRIN3D structure viewer",
  account = if (is.na(account)) NULL else account,
  forceUpdate = TRUE, launch.browser = FALSE)
