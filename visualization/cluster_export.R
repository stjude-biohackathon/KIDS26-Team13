# Cluster export: one SVG figure (3D snapshot + legend + statistics + feature
# map) and a one-row CSV of the cluster's statistics, zipped together.
# The 3D view is WebGL, so it is embedded as a PNG taken from the viewer's
# canvas; everything else in the figure is vector.

# Feature-map styles, inlined so the exported SVG renders on its own.
fm_export_css <- "
  .fm-label { font-size:11px; fill:#6B6A66; }
  .fm-in { font-size:10.5px; fill:#1F1F1D; text-anchor:middle; }
  .fm-tick { font-size:10px; fill:#8A8984; text-anchor:middle; }
  .fm-axis { stroke:#C9CCD2; stroke-width:1; }
  .fm-stem { stroke:#B9BCC3; stroke-width:1; }
  .fm-sel { fill:#2F3B52; fill-opacity:.07; }
  .fm-selres { stroke:#2F3B52; stroke-opacity:.35; stroke-width:1; }
"

fmt_p <- function(x) formatC(x, digits = 2, format = "g")

# One row: every cluster column plus what the viewer derives for it.
cluster_stats_row <- function(d, f, location) {
  s <- f$s
  co <- unique(f$co[, c("alteration_type", "cluster_id")])
  sites <- if (nrow(f$sites)) vapply(split(f$sites, f$sites$name), function(x)
    sprintf("%s (%s)", x$name[1], compact_residues(paste(unique(x$residue), collapse = ";"))), "") else character()
  span <- s$residue_max - s$residue_min
  data.frame(
    protein = d$info$protein, uniprot = d$info$uniprot,
    s[, setdiff(names(s), c("protein", "span_aa", "compaction", "colocated"))],
    location = location,
    linear_span_aa = span,
    compaction = if (s$source == "mutation" && s$diameter_3d > 0) round(span * 3.8 / s$diameter_3d, 2) else NA,
    domains = cluster_domains(d$ann, s$cluster_id),
    colocated_clusters = paste(sprintf("%s:%s", co$alteration_type, co$cluster_id), collapse = "; "),
    functional_sites_nearby = paste(sites, collapse = "; "),
    near_cutoff_angstrom = near_cutoff,
    exported = format(Sys.time(), "%Y-%m-%d %H:%M"),
    check.names = FALSE, row.names = NULL)
}

cluster_figure_svg <- function(d, f, png, fm_svg, location) {
  s <- f$s
  e <- esc
  W <- 1400
  out <- character()
  add <- function(...) out <<- c(out, sprintf(...))
  text <- function(x, y, label, size = 13, fill = ink, weight = 400, anchor = "start", family = "") add(
    '<text x="%.1f" y="%.1f" font-size="%g" fill="%s" font-weight="%d" text-anchor="%s"%s>%s</text>',
    x, y, size, fill, weight, anchor, if (nzchar(family)) sprintf(' font-family="%s"', family) else "", e(label))
  swatch <- function(x, y, t) {
    if (t %in% mut_types) add('<circle cx="%.1f" cy="%.1f" r="6" fill="%s"/>', x + 7, y - 4, layer_colors[[t]])
    else add('<rect x="%.1f" y="%.1f" width="16" height="7" rx="3.5" fill="%s"/>', x, y - 8, layer_colors[[t]])
  }
  mono <- "ui-monospace, Menlo, monospace"
  muted <- "#6B6A66"

  # Title
  text(24, 40, sprintf("%s · %s", d$info$protein, s$cluster_id), size = 24, weight = 700)
  text(24, 64, paste(c(sprintf("AlphaFold %s", d$info$uniprot), location, s$hotspot_class), collapse = "  ·  "), size = 14, fill = muted)

  # 3D snapshot, scaled to an 860 px column
  img_w <- 860
  img_h <- img_w * png$h / png$w
  y0 <- 88
  add('<rect x="24" y="%.1f" width="%.1f" height="%.1f" rx="10" fill="#FFFFFF" stroke="#E2E4E8"/>', y0, img_w, img_h)
  add('<image x="24" y="%.1f" width="%.1f" height="%.1f" href="%s" preserveAspectRatio="xMidYMid meet"/>', y0, img_w, img_h, png$uri)

  # Legend over the snapshot's top-left corner
  types <- unique(c(s$alteration_type, f$co$alteration_type))
  rows <- c(length(types), NROW(f$dom), length(unique(f$sites$kind)))
  lh <- 22
  add('<rect x="36" y="%.1f" width="250" height="%.1f" rx="8" fill="#FFFFFF" fill-opacity="0.94" stroke="#E2E4E8"/>',
      y0 + 12, 34 + sum(rows) * lh)
  text(48, y0 + 32, "AT THIS SITE", size = 10, fill = muted, weight = 600)
  y <- y0 + 54
  for (t in types) {
    swatch(48, y, t)
    text(72, y, layer_labels[[t]])
    if (t == s$alteration_type) text(276, y, "selected", size = 11, fill = muted, anchor = "end")
    y <- y + lh
  }
  for (i in seq_len(NROW(f$dom))) {
    add('<rect x="48" y="%.1f" width="14" height="14" rx="4" fill="%s" fill-opacity="0.8"/>', y - 11, f$dom$color[i])
    text(72, y, f$dom$feature_name[i]); y <- y + lh
  }
  for (k in unique(f$sites$kind)) {
    add('<rect x="48" y="%.1f" width="14" height="14" rx="4" fill="%s" fill-opacity="0.8"/>', y - 11, site_color)
    text(72, y, k); y <- y + lh
  }

  # Statistics column
  x <- 910
  y <- y0 + 8
  swatch(x, y + 14, s$alteration_type)
  text(x + 24, y + 14, layer_labels[[s$alteration_type]], size = 15, weight = 600)
  text(W - 24, y + 14, if (s$significant_any) "Significant" else "Not significant", size = 13, anchor = "end",
       fill = if (s$significant_any) "#11642F" else muted, weight = 600)
  y <- y + 40
  kpi <- list(c(s$n_subjects, "Subjects"), c(s$n_events, "Events"), c(round(s$diameter_3d, 1), "3D diameter (Å)"))
  for (i in seq_along(kpi)) {
    kx <- x + (i - 1) * 158
    add('<rect x="%.1f" y="%.1f" width="148" height="62" rx="8" fill="#F7F8FA"/>', kx, y)
    text(kx + 12, y + 30, kpi[[i]][1], size = 22, weight = 650)
    text(kx + 12, y + 50, kpi[[i]][2], size = 11, fill = muted)
  }
  y <- y + 92
  span <- s$residue_max - s$residue_min
  kv <- list(
    c("Location", location),
    c("Residues", sprintf("%d (aa %d–%d)", s$n_residues, s$residue_min, s$residue_max)),
    c("Linear span", sprintf("%d aa", span)),
    c("Compaction", if (s$source == "cna") "n/a for exon clusters" else if (s$diameter_3d > 0)
      sprintf("%.1f× (%d aa, %.1f Å across)", span * 3.8 / s$diameter_3d, span, s$diameter_3d) else "single position"),
    c("Significant in", paste(c(if (s$significant_1d) "sequence (1D)", if (s$significant_3d) "structure (3D)"), collapse = " and ")),
    c("Mean pLDDT", as.character(s$mean_plddt)),
    c("Domain", cluster_domains(d$ann, s$cluster_id)))
  for (r in kv) {
    if (!nzchar(r[2])) r[2] <- "—"
    text(x, y, r[1], fill = muted)
    text(x + 130, y, r[2])
    y <- y + 22
  }
  y <- y + 14
  text(x, y, "EMPIRICAL P-VALUES", size = 10, fill = muted, weight = 600)
  y <- y + 22
  text(W - 110, y, "1D", size = 12, fill = muted, anchor = "end"); text(W - 24, y, "3D", size = 12, fill = muted, anchor = "end")
  for (r in list(c("Size-specific", fmt_p(s$p_1d_size), fmt_p(s$p_3d_size)),
                 c("Joint", fmt_p(s$p_1d_joint), fmt_p(s$p_3d_joint)))) {
    y <- y + 22
    text(x, y, r[1]); text(W - 110, y, r[2], anchor = "end"); text(W - 24, y, r[3], anchor = "end")
  }
  y <- y + 22
  text(x, y, "Joint, either", weight = 700); text(W - 24, y, fmt_p(s$p_any_joint), anchor = "end", weight = 700)

  # Other clusters and sites at this 3D site
  co <- unique(f$co[, c("alteration_type", "cluster_id")])
  site_names <- unique(f$sites$name)
  if (nrow(co) || length(site_names)) {
    y <- y + 36
    text(x, y, "AT THIS 3D SITE", size = 10, fill = muted, weight = 600)
    for (i in seq_len(nrow(co))) {
      y <- y + 22
      swatch(x, y, co$alteration_type[i])
      text(x + 24, y, layer_labels[[co$alteration_type[i]]])
      text(W - 24, y, co$cluster_id[i], size = 12, anchor = "end", family = mono)
    }
    for (n in head(site_names, 6)) {
      y <- y + 22
      add('<rect x="%.1f" y="%.1f" width="12" height="12" rx="3" fill="%s"/>', x + 2, y - 10, site_color)
      text(x + 24, y, n)
    }
  }

  # Feature map, full width, below both columns
  fy <- max(y0 + img_h, y) + 28
  vb <- as.numeric(regmatches(fm_svg, regexpr('viewBox="0 0 [0-9.]+ [0-9.]+"', fm_svg)) |>
                     sub(pattern = 'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', replacement = "\\1 \\2") |>
                     strsplit(" ") |> unlist())
  fm_w <- W - 48
  fm_h <- fm_w * vb[2] / vb[1]
  fm <- sub("<svg ", sprintf('<svg x="24" y="%.1f" width="%.1f" height="%.1f" ', fy, fm_w, fm_h), fm_svg, fixed = TRUE)
  H <- fy + fm_h + 36
  text(24, H - 14, sprintf("GRIN3D · exported %s · members of the site are each significant on their own; the combined site is not a separate test",
                           format(Sys.Date())), size = 11, fill = muted)

  sprintf(paste0('<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="%d" height="%.0f" viewBox="0 0 %d %.0f" ',
                 'font-family="system-ui, -apple-system, Helvetica, Arial, sans-serif">',
                 '<style>%s</style><rect width="100%%" height="100%%" fill="#FFFFFF"/>%s%s</svg>'),
          W, H, W, H, fm_export_css, paste(out, collapse = ""), fm)
}
