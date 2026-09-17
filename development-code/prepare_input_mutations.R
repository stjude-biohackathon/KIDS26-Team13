
# ==============================================================================
# CREATE GRIN3D MUTATION INPUT FOR ONE PROTEIN
# ==============================================================================

# Change this for another gene.
gene.name <- "PTEN"
output.file <- paste0(gene.name, "_start_end.csv")

# Use the dataset already loaded in R.
mutations <- mutations_dataset

if (!is.data.frame(mutations)) {
  stop("mutations_dataset must be a data frame.")
}

required.columns <- c("gene", "var_type", "aa_change", "sample")
missing.columns <- setdiff(required.columns, names(mutations))

if (length(missing.columns)) {
  stop(
    "Missing required column(s): ",
    paste(missing.columns, collapse = ", ")
  )
}

# Select records annotated to this gene.
is.gene <- !is.na(mutations$gene) &
  toupper(trimws(as.character(mutations$gene))) == toupper(gene.name)

gene.records <- mutations[is.gene, , drop = FALSE]

# Select SNVs.
is.snv <- !is.na(gene.records$var_type) &
  toupper(trimws(as.character(gene.records$var_type))) == "SNV"

snvs <- gene.records[is.snv, , drop = FALSE]

# Accept single-residue substitutions (p.N31H) and stop-gain SNVs (p.R233X).
# Rows with "." or another unrecognized protein annotation are excluded.
aa.change <- toupper(trimws(as.character(snvs$aa_change)))

has.position <- !is.na(aa.change) & grepl(
  "^P\\.[ACDEFGHIKLMNPQRSTVWY][0-9]+[ACDEFGHIKLMNPQRSTVWYX]$",
  aa.change
)

selected <- snvs[has.position, , drop = FALSE]
excluded <- snvs[!has.position, , drop = FALSE]
selected.aa <- aa.change[has.position]

# Extract the amino-acid position.
residue <- as.integer(sub(
  "^P\\.[ACDEFGHIKLMNPQRSTVWY]([0-9]+)[ACDEFGHIKLMNPQRSTVWYX]$",
  "\\1",
  selected.aa
))

# Each row remains one event. A subject may have more than one event.
grin3d.input <- data.frame(
  id = seq_len(nrow(selected)),
  ID = as.character(selected$sample),
  start = residue,
  end = residue,
  stringsAsFactors = FALSE
)

if (
  anyNA(grin3d.input$ID) ||
  any(!nzchar(trimws(grin3d.input$ID)))
) {
  stop("Some selected mutations have a missing sample ID.")
}

if (anyNA(grin3d.input$start)) {
  stop("Could not extract the amino-acid position for some selected SNVs.")
}

write.csv(grin3d.input, output.file, row.names = FALSE)

message("All ", gene.name, " records: ", nrow(gene.records))
message("SNVs: ", nrow(snvs))
message("SNVs with usable amino-acid positions: ", nrow(grin3d.input))
message("SNVs excluded because the position could not be read: ", nrow(excluded))
message("Saved: ", normalizePath(output.file, winslash = "/"))

print(head(grin3d.input, 10))

# Show what types of records were present before selecting SNVs.
print(table(gene.records$var_type, useNA = "ifany"))

# Review the excluded SNVs if needed.
print(
  excluded[, c("aa_change", "sample"), drop = FALSE],
  n = nrow(excluded)
)