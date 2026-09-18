# Viewer smoke test: every protein in visualization/data/proteins.csv renders its
# layer legend, structure, feature map, cluster table and cluster panel in each
# backbone mode. Run from the repository root:
#   Rscript tests/test_grin3d_viewer.R
suppressMessages(library(shiny))
setwd("visualization")
genes <- read.csv("data/proteins.csv")$protein
testServer(shinyAppDir("."), {
  for (g in genes) {
    session$setInputs(protein = g, view = "all", backbone = "hotspots", cluster_table_rows_selected = NULL,
                      layers_mut = c("SNV", "INDEL", "MUT"),
                      layers_cna = c("HOMDEL", "HETDEL", "GAIN", "AMP"), layers_pooled = "ALL_CNA")
    session$elapse(300)
    stopifnot(nzchar(output$layer_controls$html), grepl("<svg", output$feature_map$html))
    invisible(output$cluster_table); invisible(output$inspector)
    for (b in c("domains", "plddt", "hotspots")) { session$setInputs(backbone = b); invisible(output$structure) }
    session$setInputs(view = "lr3d"); invisible(output$cluster_table)
    session$setInputs(view = "all", cluster_table_rows_selected = 1)
    stopifnot(grepl(table_rows()$cluster_id[1], output$inspector$html), grepl("fm-sel", output$feature_map$html))
    invisible(output$structure); invisible(output$member_table)
    if (!is.null(dat()$ann)) stopifnot(nrow(annotation_rows(dat()$ann)) > 0,
                                       grepl("features", output$annotation_sub$html))
    session$setInputs(view = "conv", cluster_table_rows_selected = NULL)
    session$setInputs(cluster_table_rows_selected = if (nrow(table_rows())) 1)
    invisible(output$cluster_table); invisible(output$structure); invisible(output$context_key)
  }
  # PTEN convergent sites: the long-range indel cluster 195/246/247 shares its site
  # with SNV clusters; the 3D-specific HOMDEL exon cluster 1/2/5 contains the R130 SNV cluster.
  session$setInputs(protein = "PTEN", view = "all", cluster_table_rows_selected = NULL)
  session$elapse(300)
  pick <- function(id) session$setInputs(cluster_table_rows_selected = match(id, table_rows()$cluster_id))
  pick("INDEL_cluster_0033")
  stopifnot("SNV" %in% focus()$co$alteration_type, all(c(246, 247) %in% focus()$co$residue),
            grepl("At this 3D site", output$inspector$html), grepl("C2", output$context_key$html))
  # Annotations moved to their own table under the feature map.
  ann_rows <- annotation_rows(dat()$ann, "INDEL_cluster_0033")
  stopifnot(nrow(ann_rows) > 0, "C2 tensin-type" %in% ann_rows$Type | any(grepl("C2", ann_rows$Feature)),
            grepl("INDEL_cluster_0033", output$annotation_sub$html))
  invisible(output$annotation_table); invisible(output$annotation_note)
  invisible(output$structure)
  pick("HOMDEL_cluster_0001")
  stopifnot("SNV_cluster_0001" %in% focus()$co$cluster_id, all(focus()$co$residue <= 164))
  invisible(output$structure)
  # Export: a 1x1 PNG stands in for the browser snapshot; the zip holds a
  # well-formed SVG and a one-row stats CSV.
  png <- "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNkYAAAAAYAAjCB0C8AAAAASUVORK5CYII="
  session$setInputs(structure_png = list(uri = png, w = 1000, h = 560))
  zipped <- output$export_zip
  files <- utils::unzip(zipped, exdir = tempfile())
  svg <- xml2::read_xml(grep("svg$", files, value = TRUE))
  stats <- read.csv(grep("csv$", files, value = TRUE))
  stopifnot(nrow(stats) == 1, stats$cluster_id == "HOMDEL_cluster_0001",
            grepl("SNV_cluster_0001", stats$colocated_clusters),
            grepl("HOMDEL_cluster_0001", xml2::xml_text(svg)))
  session$setInputs(view = "conv", cluster_table_rows_selected = NULL)
  stopifnot(all(c("INDEL_cluster_0033", "HOMDEL_cluster_0001") %in% table_rows()$cluster_id))
  invisible(output$structure)
})
source("annotation_track.R")
d <- select_domains(read.csv("data/annotations/FBXW7/protein_features.csv"))
stopifnot(identical(d$feature_name, c("F-box", paste("WD", 1:7))))
cat("viewer tests passed\n")
