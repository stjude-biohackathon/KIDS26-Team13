# Protein annotations for the viewer: reads the outputs of
# development-code/GRIN3D_annotate_protein_clusters.R (written per protein to
# data/annotations/<PROTEIN>/ by grin3d_viewer_input.R) and draws the linear
# feature map shown under the structure.

# Muted, so domains read as context rather than as another alteration layer.
domain_palette <- c("#B7A6DD", "#8CC4CC", "#D9C29E", "#C4A9C3", "#A7BFD9", "#C9C99A")

feature_positions <- function(segments) {
  unlist(lapply(strsplit(segments, ";", fixed = TRUE)[[1]], function(s) {
    b <- as.integer(strsplit(s, "-", fixed = TRUE)[[1]])
    seq(b[1], b[length(b)])
  }))
}

# One row per domain drawn on the map and in the 3D "Domains" mode. UniProt's
# curated domains/repeats come first; InterPro domains are added only where they
# are not already (mostly) covered, so each region is drawn once.
select_domains <- function(f) {
  uni <- f[f$source == "UniProtKB" &
             f$feature_type %in% c("Domain", "Repeat", "Zinc finger", "DNA binding", "Transmembrane"), ]
  ipr <- f[f$source == "InterPro" & f$feature_type == "domain", ]
  ipr <- ipr[order(ipr$start - ipr$end), ]  # longest first
  covered <- unlist(Map(seq, uni$start, uni$end))
  keep <- logical(nrow(ipr))
  for (i in seq_len(nrow(ipr))) {
    span <- seq(ipr$start[i], ipr$end[i])
    if (mean(span %in% covered) < 0.5) { keep[i] <- TRUE; covered <- c(covered, span) }
  }
  d <- rbind(uni, ipr[keep, ])
  d <- d[order(d$start), c("feature_id", "source", "feature_type", "feature_name", "start", "end", "source_url")]
  d$feature_name[d$feature_type == "Transmembrane"] <- "Transmembrane"
  d$group <- sub("\\s*\\d+$", "", d$feature_name)  # "WD 1".."WD 7" share a color
  groups <- unique(d$group)
  d$color <- domain_palette[(match(d$group, groups) - 1) %% length(domain_palette) + 1]
  d
}

load_annotations <- function(protein, data_dir = "data") {
  dir <- file.path(data_dir, "annotations", protein)
  if (!file.exists(file.path(dir, "protein_features.csv"))) return(NULL)
  f <- read.csv(file.path(dir, "protein_features.csv"))
  ov <- read.csv(file.path(dir, "cluster_feature_overlaps.csv"))
  ov <- ov[, c("cluster_id", "feature_id", "annotation_class", "feature_subcategory",
               "feature_type", "feature_name", "annotation_source", "relationship",
               "overlap_residues", "n_overlap_residues", "source_url")]
  chembl_file <- file.path(dir, "chembl_drug_mechanisms.csv")
  list(
    features = f,
    domains = select_domains(f),
    overlaps = ov,
    chembl = if (file.exists(chembl_file) && length(readLines(chembl_file, n = 2)) > 1)
      read.csv(chembl_file) else NULL
  )
}

# Functional sites drawn in the 3D focus view: curated active/binding sites and
# PDBe ligand contacts, plus protein-interaction interfaces. Violet passes the
# palette validator against all six alteration colors (cyan clashes with HETDEL).
site_color <- "#8A5CD0"
functional_sites <- function(ann) {
  none <- data.frame(residue = integer(), kind = character(), name = character())
  if (is.null(ann)) return(none)
  f <- ann$features[ann$features$annotation_class %in% c("site", "interaction"), ]
  if (!nrow(f)) return(none)
  seg <- ifelse(is.na(f$feature_segments) | f$feature_segments == "",
                paste(f$start, f$end, sep = "-"), f$feature_segments)
  do.call(rbind, lapply(seq_len(nrow(f)), function(i) data.frame(
    residue = feature_positions(seg[i]),
    kind = if (f$annotation_class[i] == "site") "Ligand or active site" else "Interaction interface",
    name = f$feature_name[i])))
}

# Domain names overlapping each cluster, for the table.
cluster_domains <- function(ann, cluster_ids) {
  if (is.null(ann)) return(rep("", length(cluster_ids)))
  d <- ann$domains
  ov <- ann$overlaps[ann$overlaps$feature_id %in% d$feature_id, ]
  vapply(cluster_ids, function(id) {
    hit <- d[d$feature_id %in% ov$feature_id[ov$cluster_id == id], ]
    if (!nrow(hit)) return("")
    n <- table(hit$group)[unique(hit$group)]
    paste(ifelse(n > 1, sprintf("%s \u00d7%d", names(n), n), hit$feature_name[match(names(n), hit$group)]), collapse = ", ")
  }, "")
}

esc <- function(x) htmltools::htmlEscape(as.character(x), attribute = TRUE)

# Linear feature map: sequence axis, annotation rows, exons (CNA proteins) and
# mutation lollipops; the selected cluster is shaded across every row.
feature_map_svg <- function(len, ann, residues, visible_mut, colors, selected = NULL) {
  W <- 1000; x0 <- 118; x1 <- 988
  sx <- function(r) x0 + (r - 0.5) / len * (x1 - x0)
  rows <- list(); y <- 22
  add_row <- function(label, h, body) { rows[[length(rows) + 1]] <<- list(label = label, y = y, h = h, body = body); y <<- y + h + 6 }
  bar <- function(s, e, yy, h, fill, title, label = "", stroke = "none") {
    w <- max(sx(e + 0.5) - sx(s - 0.5), 1.5)
    fit <- floor((w - 6) / 6.2)  # characters that fit; shorten with an ellipsis rather than drop
    if (nzchar(label) && nchar(label) > fit) label <- if (fit >= 4) paste0(substr(label, 1, fit - 1), "\u2026") else ""
    txt <- if (nzchar(label))
      sprintf('<text x="%.1f" y="%.1f" class="fm-in">%s</text>', sx(s - 0.5) + w / 2, yy + h / 2 + 4, esc(label)) else ""
    sprintf('<g><title>%s</title><rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="3" fill="%s" stroke="%s"/>%s</g>',
            esc(title), sx(s - 0.5), yy, w, h, fill, stroke, txt)
  }
  ticks <- function(pos, yy, h, fill, title) {
    if (!length(pos)) return("")
    paste(sprintf('<rect x="%.1f" y="%.1f" width="2" height="%.1f" fill="%s"><title>%s %d</title></rect>',
                  sx(pos) - 1, yy, h, fill, esc(title), pos), collapse = "")
  }

  # Mutation lollipops (visible mutation types; the busier type colors the head).
  mut <- residues[is.na(residues$cluster_id) & residues$alteration_type %in% visible_mut, ]
  if (nrow(mut)) {
    tot <- tapply(mut$n_events, mut$residue, sum)
    mut <- mut[order(mut$residue, -mut$n_events), ]
    mut <- mut[!duplicated(mut$residue), ]
    h <- 58; base <- y + h
    hs <- 6 + 44 * sqrt(tot[as.character(mut$residue)] / max(tot))
    add_row("Mutations", h, paste0(
      paste(sprintf('<line x1="%.1f" x2="%.1f" y1="%.1f" y2="%.1f" class="fm-stem"/>', sx(mut$residue), sx(mut$residue), base, base - hs), collapse = ""),
      paste(sprintf('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" stroke="#fff" stroke-width="1"><title>Residue %d: %d events</title></circle>',
                    sx(mut$residue), base - hs, 2.5 + 2 * sqrt(tot[as.character(mut$residue)] / max(tot)),
                    colors[mut$alteration_type], mut$residue, tot[as.character(mut$residue)]), collapse = "")))
  }

  # Exons, for proteins with CNA results.
  ex <- unique(residues[!is.na(residues$exon_order), c("exon_order", "residue")])
  if (nrow(ex)) {
    rng <- do.call(rbind, lapply(split(ex$residue, ex$exon_order), range))
    add_row("Exons", 14, paste(vapply(seq_len(nrow(rng)), function(i)
      bar(rng[i, 1], rng[i, 2], y, 14, if (i %% 2) "#E6E8EC" else "#D3D7DE",
          sprintf("Exon %s: aa %d-%d", rownames(rng)[i], rng[i, 1], rng[i, 2]), rownames(rng)[i]), ""), collapse = ""))
  }

  if (!is.null(ann)) {
    f <- ann$features; d <- ann$domains
    if (nrow(d)) add_row("Domains", 18, paste(vapply(seq_len(nrow(d)), function(i)
      bar(d$start[i], d$end[i], y, 18, d$color[i],
          sprintf("%s (%s %s, aa %d-%d)", d$feature_name[i], d$source[i], d$feature_type[i], d$start[i], d$end[i]),
          d$feature_name[i]), ""), collapse = ""))
    reg <- f[f$annotation_class == "motif" | f$feature_subcategory == "intrinsically_disordered_region" |
               f$feature_type == "Topological domain", ]
    if (nrow(reg)) add_row("Regions", 12, paste(vapply(seq_len(nrow(reg)), function(i)
      bar(reg$start[i], reg$end[i], y, 12,
          if (reg$annotation_class[i] == "motif") "#6F7A8C" else if (reg$feature_type[i] == "Topological domain") "#EEF0F3" else "#C9CCD2",
          sprintf("%s: %s (aa %d-%d)", reg$feature_type[i], reg$feature_name[i], reg$start[i], reg$end[i]),
          stroke = if (reg$feature_type[i] == "Topological domain") "#C9CCD2" else "none"), ""), collapse = ""))
    site <- f[f$annotation_class == "site" & f$source != "PDBe-KB/SIFTS", ]
    if (nrow(site)) add_row("Sites", 12, paste(vapply(seq_len(nrow(site)), function(i)
      bar(site$start[i], site$end[i], y, 12, "#8A5A00",
          sprintf("%s: %s (aa %d-%d)", site$feature_type[i], site$feature_name[i], site$start[i], site$end[i])), ""), collapse = ""))
    ptm <- f[f$annotation_class == "modification", ]
    if (nrow(ptm)) add_row("PTMs", 12, paste(vapply(seq_len(nrow(ptm)), function(i)
      sprintf('<circle cx="%.1f" cy="%.1f" r="3" fill="#5B5A56"><title>%s at %d</title></circle>',
              sx(ptm$start[i]), y + 6, esc(ptm$feature_name[i]), ptm$start[i]), ""), collapse = ""))
    pdbe <- f[f$source == "PDBe-KB/SIFTS", ]
    if (nrow(pdbe)) {
      pos_of <- function(rows) sort(unique(unlist(lapply(rows$feature_segments, feature_positions))))
      iface <- pos_of(pdbe[pdbe$annotation_class == "interaction", ])
      lig <- pos_of(pdbe[pdbe$annotation_class == "site", ])
      add_row("Contacts", 14, paste0(ticks(iface, y, 6, "#7D8CA3", "Interface contact (PDBe-KB) at"),
                                     ticks(lig, y + 8, 6, "#B07A2A", "Ligand contact (PDBe-KB) at")))
    }
  }

  H <- y + 18
  # Axis at the bottom.
  step <- c(10, 25, 50, 100, 200)[findInterval(len / 8, c(0, 20, 45, 90, 180))]
  at <- c(1, seq(step, len, by = step))
  axis <- paste0(sprintf('<line x1="%d" x2="%d" y1="%d" y2="%d" class="fm-axis"/>', x0, x1, y, y),
                 paste(sprintf('<text x="%.1f" y="%d" class="fm-tick">%d</text>', sx(at), y + 13, at), collapse = ""))
  shade <- ""
  if (!is.null(selected)) {
    r <- sort(unique(residues$residue[residues$cluster_id %in% selected$cluster_id]))
    shade <- paste0(
      sprintf('<rect x="%.1f" y="8" width="%.1f" height="%.1f" class="fm-sel"/>',
              sx(min(r) - 0.5), max(sx(max(r) + 0.5) - sx(min(r) - 0.5), 2), y - 8),
      paste(sprintf('<line x1="%.1f" x2="%.1f" y1="8" y2="%.1f" class="fm-selres"/>', sx(r), sx(r), y), collapse = ""))
  }
  labels <- paste(vapply(rows, function(r) sprintf('<text x="8" y="%.1f" class="fm-label">%s</text>', r$y + r$h / 2 + 4, r$label), ""), collapse = "")
  bodies <- paste(vapply(rows, `[[`, "", "body"), collapse = "")
  sprintf('<svg viewBox="0 0 %d %.0f" class="feature-map" role="img" aria-label="Protein feature map">%s%s%s%s</svg>',
          W, H, shade, labels, bodies, axis)
}

table_class_order <- c("site", "interaction", "functional_region", "motif",
                       "modification", "structural_region", "domain")
ann_groups <- c(domain = "Domains", site = "Sites and ligand contacts", interaction = "Interactions",
                modification = "Modifications", motif = "Motifs",
                functional_region = "Functional regions", structural_region = "Structural regions")
source_short <- c(UniProtKB = "UniProt", InterPro = "InterPro", `PDBe-KB/SIFTS` = "PDBe-KB")

compact_residues <- function(x) {
  r <- sort(unique(as.integer(unlist(strsplit(x, ";", fixed = TRUE)))))
  if (length(r) > 4) sprintf("%d res. (aa %d–%d)", length(r), min(r), max(r)) else paste("aa", paste(r, collapse = ", "))
}

# Cluster panel section: features from the annotation module that overlap the
# cluster's residues, grouped by annotation class, each linked to its source.
annotation_section <- function(ann, cluster_id, n_residues, is_cna) {
  if (is.null(ann)) return(NULL)
  ov <- ann$overlaps[ann$overlaps$cluster_id == cluster_id, ]
  hits <- length(unique(unlist(strsplit(ov$overlap_residues, ";", fixed = TRUE))))
  tagList(
    div(class = "eyebrow", style = "margin-top:14px", "Biological annotations"),
    if (!nrow(ov)) p(class = "hint", "No annotated features overlap this cluster.") else tagList(
      p(class = "hint", sprintf("%d of %d %s fall in %d annotated features.", hits, n_residues,
                                if (is_cna) "exon-encoded residues" else "residues", length(unique(ov$feature_id)))),
      lapply(intersect(names(ann_groups), unique(ov$annotation_class)), function(g) {
        x <- ov[ov$annotation_class == g, ]
        x <- x[order(-x$n_overlap_residues, x$feature_name), ]
        item <- function(i) {
          dom <- match(x$feature_id[i], ann$domains$feature_id)
          div(class = "ann-item",
              if (!is.na(dom)) span(class = "dom-sw", style = sprintf("background:%s", ann$domains$color[dom])),
              tags$a(href = x$source_url[i], target = "_blank", rel = "noopener", x$feature_name[i]),
              span(class = "sub", source_short[x$annotation_source[i]]),
              span(class = "res", compact_residues(x$overlap_residues[i])))
        }
        shown <- seq_len(min(nrow(x), 5))
        div(class = "ann-group",
            div(class = "hdr", span(ann_groups[[g]]), span(nrow(x))),
            lapply(shown, item),
            if (nrow(x) > 5) tags$details(tags$summary(class = "hint", sprintf("%d more", nrow(x) - 5)),
                                          lapply(setdiff(seq_len(nrow(x)), shown), item)))
      })),
    p(class = "scope-note", "UniProt, InterPro and PDBe-KB features overlapping the cluster.",
      if (is_cna) "For CNA clusters this means residues encoded by the cluster's exons.",
      "Annotations add context only; they do not change significance.")
  )
}

# Rows for the annotation table under the feature map: the features overlapping
# the selected cluster, or every annotated feature when nothing is selected.
annotation_rows <- function(ann, cluster_id = NULL) {
  if (is.null(ann)) return(NULL)
  link <- function(name, url) sprintf('<a href="%s" target="_blank" rel="noopener">%s</a>', esc(url), esc(name))
  dot <- function(ids) ifelse(is.na(match(ids, ann$domains$feature_id)), "",
    sprintf('<i class="dom-sw" style="background:%s"></i> ', ann$domains$color[match(ids, ann$domains$feature_id)]))
  if (is.null(cluster_id)) {
    f <- ann$features[order(ann$features$start, ann$features$end), ]
    return(data.frame(
      Class = ann_groups[f$annotation_class],
      Feature = paste0(dot(f$feature_id), link(f$feature_name, f$source_url)),
      Type = f$feature_type, Source = source_short[f$source],
      Position = ifelse(f$start == f$end, paste("aa", f$start), sprintf("aa %d–%d", f$start, f$end)),
      Residues = vapply(f$feature_segments, function(x) length(feature_positions(x)), 1L),
      check.names = FALSE, row.names = NULL))
  }
  ov <- ann$overlaps[ann$overlaps$cluster_id == cluster_id, ]
  if (!nrow(ov)) return(data.frame())
  # Most specific first: a cluster's sites and interactions say more than the
  # domain rows, which repeat the same region across InterPro entries.
  ov <- ov[order(match(ov$annotation_class, table_class_order), -ov$n_overlap_residues, ov$feature_name), ]
  data.frame(
    Class = ann_groups[ov$annotation_class],
    Feature = paste0(dot(ov$feature_id), link(ov$feature_name, ov$source_url)),
    Type = ov$feature_type, Source = source_short[ov$annotation_source],
    `Overlapping residues` = vapply(ov$overlap_residues, compact_residues, ""),
    Residues = ov$n_overlap_residues,
    check.names = FALSE, row.names = NULL)
}

# Protein-level drug context for the empty cluster panel.
chembl_section <- function(ann) {
  m <- ann$chembl
  if (is.null(m) || !nrow(m)) return(NULL)
  m <- m[order(-suppressWarnings(as.numeric(m$max_phase)), m$molecule_name), ]
  m <- m[!duplicated(m$molecule_chembl_id), ]
  tagList(
    div(class = "eyebrow", style = "margin-top:18px", sprintf("Drugs targeting this protein (ChEMBL, %d)", nrow(m))),
    lapply(seq_len(min(nrow(m), 8)), function(i) div(class = "ann-item",
      tags$a(href = m$molecule_url[i], target = "_blank", rel = "noopener", tools::toTitleCase(tolower(m$molecule_name[i]))),
      span(class = "sub", tolower(m$action_type[i])),
      span(class = "res", if (is.na(m$max_phase[i])) "" else sprintf("phase %s", m$max_phase[i])))),
    p(class = "scope-note", "Target-level evidence for the whole protein, not for any cluster.")
  )
}
