# GRIN3D structure viewer: mutation and CNA-exon hotspots on the AlphaFold model.
# Run from the repository root:  shiny::runApp("visualization")
# Inputs are the standard tables written by grin3d_viewer_input.R, one set per
# protein listed in data/proteins.csv.

library(shiny)
library(r3dmol)
library(DT)
library(bslib)

source("grin3d_viewer_input.R")  # read_plddt()
source("annotation_track.R")       # protein annotations and the feature map
source("cluster_export.R", local = TRUE)  # SVG figure + stats CSV; local so it sees the app globals below

# Proteins come from data/proteins.csv, written by grin3d_viewer_input.R.
proteins <- read.csv("data/proteins.csv")
proteins$types <- strsplit(proteins$types, ";", fixed = TRUE)

near_cutoff <- 8  # Å between C-alpha atoms: "at the same 3D site"

load_protein <- function(protein) {
  m <- proteins[proteins$protein == protein, ]
  prefix <- file.path("data", paste0(protein, "_viewer"))
  # Deployed bundles carry the structures in data/structures (see deploy.R);
  # locally they are read from the examples folders the manifest points at.
  bundled <- file.path("data", "structures", basename(m$pdb_file))
  pdb_file <- if (file.exists(bundled)) bundled else file.path("..", m$pdb_file)
  coords <- read_plddt(pdb_file)
  near <- as.matrix(dist(coords[, c("x", "y", "z")])) <= near_cutoff
  dimnames(near) <- list(coords$residue, coords$residue)
  ann <- load_annotations(protein)
  list(
    info = m,
    clusters = read.csv(paste0(prefix, "_clusters.csv")),
    residues = read.csv(paste0(prefix, "_residues.csv")),
    members = read.csv(paste0(prefix, "_members.csv")),
    coords = coords,
    near = near,
    ann = ann,
    sites = functional_sites(ann),
    pdb_text = paste(readLines(pdb_file, warn = FALSE), collapse = "\n")
  )
}

# Residues within near_cutoff of any of the given residues.
near_residues <- function(d, resi) {
  resi <- intersect(as.character(resi), rownames(d$near))
  as.integer(colnames(d$near)[colSums(d$near[resi, , drop = FALSE]) > 0])
}

# Significant clusters of other alteration types at the same 3D site as cluster
# s, as (cluster_id, alteration_type, residue) rows to draw alongside it:
# - mutation seed: other mutation types' cluster residues within near_cutoff of
#   the seed, and whole CNA clusters whose exons contain a seed residue;
# - CNA seed: clusters of any other type lying wholly inside the seed's exons.
# Descriptive only: each member is significant on its own, the site is not tested.
colocated <- function(d, s, types) {
  cl <- d$clusters
  rs <- d$residues
  seed <- unique(rs$residue[rs$cluster_id %in% s$cluster_id])
  cand <- cl$cluster_id[cl$significant_any & cl$alteration_type %in% setdiff(types, c(s$alteration_type, "ALL_CNA"))]
  r <- unique(rs[rs$cluster_id %in% cand, c("cluster_id", "alteration_type", "residue")])
  is_mut <- r$alteration_type %in% mut_types
  keep <- if (s$source == "mutation") {
    ifelse(is_mut, r$residue %in% near_residues(d, seed),
           r$cluster_id %in% r$cluster_id[!is_mut & r$residue %in% seed])
  } else {
    inside <- tapply(r$residue, r$cluster_id, function(x) all(x %in% seed))
    inside[r$cluster_id]
  }
  r[keep, ]
}

# Alteration colors. Copy number is a diverging scale (blue = loss, rose = gain,
# darker = more copies changed). SNVs (amber) and indels (green) sit outside that
# scale and are drawn as a different mark (spheres, not backbone). The six colors
# pass the dataviz palette validator all-pairs on white (CVD dE >= 8.6, normal
# dE >= 16.5); SNV vs indel, which share a mark, is CVD dE 11.
layer_colors <- c(
  SNV = "#EDA100", INDEL = "#1A9A5C",
  MUT = "#EDA100",  # unsplit mutation results; never shown beside SNV/INDEL
  HOMDEL = "#184F95", HETDEL = "#6DA7EC", GAIN = "#EF8AA0", AMP = "#B3264A",
  ALL_CNA = "#5B5A56"  # pooled analysis: neutral, never shown beside the others
)
layer_labels <- c(
  SNV = "SNV", INDEL = "Indel", MUT = "Mutation (SNV + indel)",
  HOMDEL = "Homozygous deletion", HETDEL = "Heterozygous deletion",
  GAIN = "Gain", AMP = "Amplification", ALL_CNA = "Pooled CNAs"
)
mut_types <- c("SNV", "INDEL", "MUT")
cna_types <- c("HOMDEL", "HETDEL", "GAIN", "AMP")  # copy-number order, loss to gain

protein_choices <- setNames(proteins$protein, sprintf(
  "%s (%s) - %s", proteins$protein, proteins$uniprot,
  vapply(proteins$types, function(t) paste(c(
    if (any(t %in% mut_types)) "mutations",
    if (any(t %in% c(cna_types, "ALL_CNA"))) "CNAs"), collapse = " + "), "")))
ink <- "#1F1F1D"
# AlphaFold pLDDT bands
plddt_bins <- list(c(90, 101, "#0053D6"), c(70, 90, "#65CBF3"),
                   c(50, 70, "#FFDB13"), c(0, 50, "#FF7D45"))
backbone_grey <- "#D9D9D9"
focus_grey <- "#EFEFEF"  # backbone outside the selected cluster

# One legend row per layer, doubling as its toggle: swatch shaped like the mark
# (dot = mutation sphere, bar = CNA backbone), name, and cluster counts.
layer_choice <- function(t, clusters) {
  d <- clusters[clusters$alteration_type == t, ]
  tags$span(class = "layer",
    tags$i(class = if (t %in% mut_types) "sw dot" else "sw bar",
           style = sprintf("background:%s", layer_colors[[t]])),
    tags$span(class = "name", layer_labels[[t]]),
    tags$span(class = "count", if (nrow(d)) sprintf("%d sig / %d", sum(d$significant_any), nrow(d)) else "none"))
}
type_chip <- function(t) sprintf(
  '<span class="chip"><i class="sw %s" style="background:%s"></i>%s</span>',
  ifelse(t %in% mut_types, "dot", "bar"), layer_colors[t], layer_labels[t])

# A legend section of layer toggles, listing only the types this protein has.
layer_section <- function(id, label, types, clusters, ...) {
  if (!length(types)) return(NULL)
  div(class = "legend-section",
    div(class = "eyebrow", label), ...,
    checkboxGroupInput(id, NULL, choiceNames = lapply(types, layer_choice, clusters),
                       choiceValues = types, selected = types))
}

default_protein <- if ("PTEN" %in% proteins$protein) "PTEN" else proteins$protein[1]

# Neutral UI chrome; the data palette above is the only color on the page.
ui_theme <- bs_theme(
  version = 5, bg = "#F4F5F7", fg = ink, primary = "#2F3B52",
  base_font = font_collection("system-ui", "-apple-system", "Segoe UI", "sans-serif"),
  code_font = font_collection("ui-monospace", "SFMono-Regular", "Menlo", "monospace"),
  "card-bg" = "#FFFFFF", "border-color" = "#E2E4E8", "border-radius" = "10px",
  "font-size-base" = "0.9rem"
)

viewer_css <- tags$style(HTML(sprintf("
  body { padding-bottom: 32px; }
  .appbar { display:flex; align-items:center; gap:20px; flex-wrap:wrap; padding:14px 4px 6px; }
  .wordmark { font-size:22px; font-weight:700; letter-spacing:.02em; margin:0; }
  .appbar .form-group { margin:0; min-width:300px; }
  .appbar .selectize-input { border-radius:8px; }
  .protein-meta { color:#6B6A66; font-size:13px; display:flex; gap:14px; flex-wrap:wrap; }
  .protein-meta b { color:%1$s; font-weight:600; }
  .card-header { background:#FFFFFF; display:flex; align-items:center; gap:12px; flex-wrap:wrap; font-weight:600; }
  .card-header .sub { color:#6B6A66; font-weight:400; font-size:13px; }
  .card-header .tools { margin-left:auto; display:flex; gap:18px; align-items:center; font-weight:400; }
  .card-header .form-check { margin:0; }
  .card-header .shiny-input-container { margin:0; width:auto; }
  .eyebrow { font-size:11px; letter-spacing:.07em; text-transform:uppercase; color:#6B6A66; font-weight:600; margin:0 0 4px; }
  .viewer { position:relative; }
  .legend { position:absolute; top:12px; left:12px; z-index:5; width:250px; max-height:calc(100%% - 24px); overflow:auto;
            background:rgba(255,255,255,.94); border:1px solid #E2E4E8; border-radius:10px;
            box-shadow:0 4px 16px rgba(20,24,31,.08); padding:10px 12px; }
  .legend summary { cursor:pointer; font-weight:600; list-style:none; display:flex; align-items:center; justify-content:space-between; }
  .legend summary::after { content:'Hide'; font-weight:400; font-size:12px; color:#6B6A66; }
  .legend:not([open]) summary::after { content:'Show'; }
  .legend-section { margin-top:10px; }
  .legend .form-group { margin-bottom:0; }
  .legend .checkbox { margin:2px 0; }
  .legend .checkbox label { display:flex; align-items:center; gap:8px; width:100%%; padding:0; }
  .legend .checkbox input { position:static; margin:0; flex:none; }
  .legend .checkbox label > span { flex:1; min-width:0; }
  .hint { color:#6B6A66; font-size:12px; line-height:1.35; margin:2px 0 0; }
  .sw { display:inline-block; flex:none; background:#999; }
  .sw.dot { width:12px; height:12px; border-radius:50%%; box-shadow:0 0 0 2px #fff, 0 0 0 3px rgba(0,0,0,.12); }
  .sw.bar { width:18px; height:8px; border-radius:4px; }
  .layer { display:inline-flex; align-items:center; gap:8px; width:calc(100%% - 4px); }
  .layer .count { margin-left:auto; white-space:nowrap; color:#6B6A66; font-size:12px; font-variant-numeric:tabular-nums; }
  .cn-scale { display:grid; gap:2px; margin:4px 0 2px; }
  .cn-scale i { height:6px; }
  .cn-scale i:first-child { border-radius:3px 0 0 3px; } .cn-scale i:last-child { border-radius:0 3px 3px 0; }
  .cn-ends { display:flex; justify-content:space-between; font-size:11px; color:#6B6A66; }
  .chip { display:inline-flex; align-items:center; gap:6px; white-space:nowrap; color:%1$s; }
  .mono { font-family:ui-monospace, SFMono-Regular, Menlo, monospace; font-size:.85em; }
  .status { white-space:nowrap; display:inline-flex; align-items:center; gap:5px; padding:2px 9px; border-radius:999px; font-size:12px; font-weight:600; border:1px solid; }
  .status.yes { color:#11642F; border-color:#9ED3AE; background:#EEF8F1; }
  .status.no { color:#6B6A66; border-color:#E2E4E8; background:#F7F7F8; }
  .kpis { display:grid; grid-template-columns:repeat(3, 1fr); gap:8px; margin:12px 0; }
  .kpi { background:#F7F8FA; border-radius:8px; padding:8px 10px; }
  .kpi .v { font-size:20px; font-weight:650; font-variant-numeric:tabular-nums; }
  .kpi .l { font-size:11px; color:#6B6A66; }
  .kv { display:grid; grid-template-columns:auto 1fr; gap:4px 14px; font-size:13px; margin:0; }
  .kv dt { color:#6B6A66; font-weight:400; } .kv dd { margin:0; font-variant-numeric:tabular-nums; }
  .pgrid { display:grid; grid-template-columns:1fr auto auto; gap:4px 28px; font-size:13px; font-variant-numeric:tabular-nums; margin:2px 0 14px; }
  .pgrid > :not(:nth-child(3n+1)) { text-align:right; }
  .pgrid .h { color:#6B6A66; font-size:12px; }
  .dt-buttons { margin-bottom:6px; }
  .dt-buttons .btn { background:#fff; color:%1$s; border:1px solid #D5D8DE; padding:3px 10px; font-size:13px; }
  .dt-buttons .btn:hover { background:#F1F2F4; }
  .empty { color:#6B6A66; }
  .empty h5 { color:%1$s; font-weight:600; }
  .sig-list { list-style:none; padding:0; margin:10px 0 0; display:grid; gap:6px; }
  .sig-list li { display:flex; align-items:center; justify-content:space-between; }
  table.dataTable td { vertical-align:middle; }
  .seg .shiny-options-group { display:inline-flex; border:1px solid #D5D8DE; border-radius:8px; overflow:hidden; }
  .seg .radio-inline, .seg .radio-inline + .radio-inline { margin:0 !important; padding:0 !important; }
  .seg .shiny-input-container { margin:0; width:auto; }
  .seg input { position:absolute; opacity:0; pointer-events:none; }
  .seg .radio-inline span { display:block; padding:4px 12px; font-size:13px; cursor:pointer; border-left:1px solid #D5D8DE; }
  .seg .radio-inline:first-child span { border-left:0; }
  .seg input:checked + span { background:#2F3B52; color:#fff; }
  .seg input:focus-visible + span { outline:2px solid #2F3B52; outline-offset:-2px; }
  .fm-wrap { padding:6px 12px 10px; border-top:1px solid #E2E4E8; }
  .feature-map { width:100%%; height:auto; display:block; font-family:system-ui, -apple-system, sans-serif; }
  .fm-label { font-size:11px; fill:#6B6A66; }
  .fm-in { font-size:10.5px; fill:%1$s; text-anchor:middle; pointer-events:none; }
  .fm-tick { font-size:10px; fill:#8A8984; text-anchor:middle; }
  .fm-axis { stroke:#C9CCD2; stroke-width:1; }
  .fm-stem { stroke:#B9BCC3; stroke-width:1; }
  .fm-sel { fill:#2F3B52; fill-opacity:.07; }
  .fm-selres { stroke:#2F3B52; stroke-opacity:.35; stroke-width:1; }
  .ann-group { margin:6px 0 10px; }
  .ann-group .hdr { display:flex; justify-content:space-between; font-size:12px; color:#6B6A66; margin-bottom:3px; }
  .ann-item { display:flex; align-items:baseline; gap:8px; font-size:13px; line-height:1.35; padding:1px 0; }
  .ann-item .res { margin-left:auto; color:#6B6A66; font-size:12px; white-space:nowrap; font-variant-numeric:tabular-nums; }
  .ann-item a { color:%1$s; text-decoration:none; border-bottom:1px solid #D5D8DE; }
  .ann-item a:hover { border-bottom-color:%1$s; }
  .dom-sw { width:10px; height:10px; border-radius:2px; flex:none; display:inline-block; transform:translateY(1px); }
  .scope-note { font-size:11.5px; color:#6B6A66; margin:4px 0 0; }
  .frac { display:inline-block; width:104px; height:8px; border-radius:4px; background:#E9EAEE; overflow:hidden; vertical-align:middle; }
  .frac i { display:block; height:100%%; border-radius:4px; }
  table.dataTable td.frac-cell { width:120px; }
  .stack { display:flex; flex-direction:column; gap:var(--bslib-spacer, 1rem); }
  .view-note { background:#F7F8FA; border-radius:8px; padding:8px 12px; font-size:13px; color:#4A4A47; margin-bottom:8px; }
  table.dataTable td.cls { max-width:220px; white-space:normal; line-height:1.3; }
  .sw.env { width:14px; height:14px; border-radius:5px; opacity:.8; }
  .key-row { display:flex; align-items:center; gap:8px; margin:3px 0; font-size:13px; }
  .key-row .count { margin-left:auto; color:#6B6A66; font-size:12px; }
  .focus-key { padding-bottom:8px; border-bottom:1px solid #E2E4E8; }
", ink)))

ui <- page_fluid(
  theme = ui_theme,
  title = "GRIN3D",
  div(style = "display:none", downloadLink("export_zip", "download")),  # clicked by grin3dExport()
  tags$head(viewer_css, tags$script(HTML("
    // Export: snapshot the 3D canvas (re-rendered first so its buffer is fresh),
    // hand it to the server, then fetch the zip once the server has it.
    function grin3dExport() {
      var el = document.getElementById('structure');
      var c = el && el.querySelector('canvas');
      if (!c) return;
      el.widget.render();
      Shiny.setInputValue('structure_png', {uri: c.toDataURL('image/png'), w: c.width, h: c.height}, {priority: 'event'});
    }
    Shiny.addCustomMessageHandler('grin3d-download', function(id) { document.getElementById(id).click(); });
  "))),
  div(class = "appbar",
    h1(class = "wordmark", "GRIN3D"),
    selectizeInput("protein", NULL, choices = protein_choices, selected = default_protein,
                   options = list(placeholder = "Search a gene or UniProt ID")),
    uiOutput("protein_meta", class = "protein-meta")
  ),
  layout_columns(
    col_widths = breakpoints(sm = 12, lg = c(8, 4)),
    div(class = "stack",
    card(
      full_screen = TRUE,
      card_header(
        "Structure", uiOutput("structure_sub", inline = TRUE, class = "sub"),
        div(class = "tools", span(class = "sub", "Backbone"),
            div(class = "seg", radioButtons("backbone", NULL, inline = TRUE, selected = "hotspots",
              choiceNames = list("Hotspots", "Domains", "pLDDT"),
              choiceValues = c("hotspots", "domains", "plddt"))))
      ),
      card_body(
        padding = 0,
        div(class = "viewer",
          r3dmolOutput("structure", height = "560px"),
          tags$details(class = "legend", open = NA,
            tags$summary("Layers"),
            uiOutput("context_key"),
            uiOutput("layer_controls"))
        ),
        div(class = "fm-wrap", uiOutput("feature_map"))
      )
    ),
    card(
      card_header("Biological annotations", uiOutput("annotation_sub", inline = TRUE, class = "sub")),
      card_body(fillable = FALSE, DTOutput("annotation_table"), uiOutput("annotation_note"))
    )),
    card(
      card_header("Cluster", div(class = "tools", uiOutput("clear_button", inline = TRUE))),
      card_body(uiOutput("inspector"))
    )
  ),
  card(
    fill = FALSE,
    card_header(
      "Candidate clusters", uiOutput("table_sub", inline = TRUE, class = "sub"),
      div(class = "tools seg", radioButtons("view", NULL, inline = TRUE, selected = "conv",
        choiceNames = list("All", "Significant",
                           tags$span(title = "Significant in 3D but not along the sequence: residues far apart in the chain that fold close together", "Long-range 3D"),
                           tags$span(title = "Significant clusters that share a 3D site with significant clusters of other alteration types", "Convergent")),
        choiceValues = c("all", "sig", "lr3d", "conv")))
    ),
    card_body(fillable = FALSE, uiOutput("view_note"), DTOutput("cluster_table"))
  )
)

server <- function(input, output, session) {

  dat <- reactive(load_protein(req(input$protein)))
  # The layer checkboxes are rebuilt per protein and report their values a moment
  # after the protein changes; debouncing lets the structure and table render once
  # with the new values instead of first rendering with no layers.
  layers <- debounce(reactive(list(
    mut = input$layers_mut, cna = input$layers_cna, pooled = input$layers_pooled,
    all_cna = isTRUE(input$all_cna))), 250)
  visible <- reactive(with(layers(), c(mut, cna, pooled)))

  # Only the alteration types this protein actually has get a legend row;
  # sections with none are left out entirely.
  output$layer_controls <- renderUI({
    present <- dat()$info$types[[1]]
    cl <- dat()$clusters
    muts <- intersect(mut_types, present)
    cnas <- intersect(cna_types, present)
    tagList(
      layer_section("layers_mut", "Mutations", muts, cl),
      if (length(muts)) p(class = "hint", "Spheres on mutated residues, sized by events."),
      layer_section("layers_cna", "Copy number", cnas, cl,
        div(class = "cn-scale", style = sprintf("grid-template-columns:repeat(%d, 1fr)", length(cnas)),
            lapply(cnas, function(t) tags$i(style = sprintf("background:%s", layer_colors[[t]])))),
        div(class = "cn-ends", span("Loss"), span("Gain"))),
      if (length(cnas)) tagList(
        p(class = "hint", "Backbone of exons in significant hotspots."),
        div(style = "margin-top:6px", input_switch("all_cna", "All altered exons", FALSE))),
      layer_section("layers_pooled", "Pooled analysis", intersect("ALL_CNA", present), cl),
      if ("ALL_CNA" %in% present) p(class = "hint", "All CNA types tested together; drawn when selected.")
    )
  })

  table_rows <- reactive({
    d <- dat()$clusters
    d <- d[d$alteration_type %in% visible(), ]
    d$span_aa <- d$residue_max - d$residue_min
    # Linear chain length (3.8 Å per residue) over the 3D diameter: how many times
    # more compact the cluster is in the fold than stretched along the sequence.
    # Mutation clusters only: a CNA cluster's 3D diameter is an exon-to-exon
    # proximity (closest residues), not the width of the residue set.
    d$compaction <- ifelse(d$source == "mutation" & d$diameter_3d > 0,
                           d$span_aa * 3.8 / d$diameter_3d, NA)
    # Other significant alteration types at the same 3D site (significant rows only).
    d$colocated <- vapply(seq_len(nrow(d)), function(i) if (!d$significant_any[i]) "" else
      paste(unique(colocated(dat(), d[i, ], visible())$alteration_type), collapse = ";"), "")
    view <- if (is.null(input$view)) "all" else input$view
    if (view == "sig") d <- d[d$significant_any, ]
    if (view == "conv") {
      # 3D-significant clusters sharing their site with other alteration types:
      # long-range 3D first, then by number of types and compaction.
      d <- d[nzchar(d$colocated) & d$significant_3d, ]
      lr <- d$significant_3d & !d$significant_1d
      return(d[order(!lr, -lengths(strsplit(d$colocated, ";")), -ifelse(is.na(d$compaction), 0, d$compaction)), ])
    }
    if (view == "lr3d") {
      d <- d[d$significant_3d & !d$significant_1d, ]
      return(d[order(-d$compaction, d$p_3d_joint), ])
    }
    d[order(d$p_any_joint, -d$n_subjects), ]
  })

  # The selection is held as a cluster id, not a row index, so hiding a layer or
  # switching the table view keeps the structure focused on the same cluster
  # instead of zooming back out. It is cleared only by "Show all layers" or by
  # changing protein.
  selected_id <- reactiveVal(NULL)
  observeEvent(input$cluster_table_rows_selected, {
    i <- input$cluster_table_rows_selected
    d <- table_rows()
    if (length(i) && i <= nrow(d)) selected_id(d$cluster_id[i])
  })
  observeEvent(input$protein, selected_id(NULL))

  selected <- reactive({
    id <- selected_id()
    d <- dat()$clusters
    r <- d[d$cluster_id %in% id, ]
    if (!nrow(r)) return(NULL)  # the id can briefly belong to the previous protein
    r$span_aa <- r$residue_max - r$residue_min
    r$compaction <- ifelse(r$source == "mutation" & r$diameter_3d > 0, r$span_aa * 3.8 / r$diameter_3d, NA)
    r
  })

  # Keep the table's highlight on the selected cluster when it is still listed.
  observe({
    d <- table_rows()
    i <- match(selected_id(), d$cluster_id)
    if (!is.null(selected_id()) && !identical(i, input$cluster_table_rows_selected))
      selectRows(dataTableProxy("cluster_table"), if (is.na(i)) NULL else i)
  })

  # "Exons 1, 2, 5" for CNA clusters, "aa 232-246" for mutation clusters.
  cluster_location <- function(d) {
    r <- dat()$residues
    rng <- tapply(r$residue, r$cluster_id, function(x) sprintf("aa %d–%d", min(x), max(x)))
    ifelse(d$source == "cna", paste("Exons", gsub(";", ", ", d$exons)),
           sub("aa (\\d+)–\\1$", "aa \\1", rng[d$cluster_id]))
  }

  output$protein_meta <- renderUI({
    m <- dat()$info
    tagList(span("UniProt ", tags$b(m$uniprot)),
            span(tags$b(nrow(dat()$coords)), " residues"),
            span(tags$b(sum(dat()$clusters$significant_any)), " significant of ",
                 tags$b(nrow(dat()$clusters)), " clusters"))
  })
  # Names the focused cluster (and its exons) here rather than as a 3D callout that collides with residue labels.
  output$structure_sub <- renderUI({
    f <- focus()
    base <- sprintf("%s · AlphaFold %s", dat()$info$protein, dat()$info$uniprot)
    if (is.null(f)) return(base)
    tagList(base, " · ", tags$b(style = sprintf("color:%s", ink), f$s$cluster_id),
            if (length(f$exons)) sprintf(" · exon%s %s", if (length(f$exons) > 1) "s" else "", paste(f$exons, collapse = ", ")))
  })
  output$feature_map <- renderUI({
    HTML(feature_map_svg(nrow(dat()$coords), dat()$ann, dat()$residues,
                         intersect(layers()$mut, mut_types), layer_colors, selected()))
  })

  output$view_note <- renderUI(switch(if (is.null(input$view)) "all" else input$view,
    lr3d = div(class = "view-note",
      "Clusters significant in 3D but not along the sequence, sorted by compaction: ",
      "linear span × 3.8 Å per residue ÷ 3D diameter (mutation clusters; CNA clusters follow). 10× means the residues sit ten times closer in the fold than along the stretched chain."),
    conv = div(class = "view-note",
      sprintf("Clusters significant in 3D that share their site with significant clusters of other alteration types (mutations within %g Å; CNA clusters by shared exons). ", near_cutoff),
      "Long-range 3D clusters first. Select one to draw every type at that site. Each member is significant on its own; the combined site is not a separate test."),
    NULL))

  output$table_sub <- renderUI({
    d <- table_rows()
    sprintf("%d shown · %d significant · select a row to focus", nrow(d), sum(d$significant_any))
  })

  output$cluster_table <- renderDT({
    d <- table_rows()
    out <- data.frame(
      Cluster = sprintf('<span class="mono">%s</span>', d$cluster_id),
      Type = type_chip(d$alteration_type),
      Sig = ifelse(d$significant_any, '<span class="status yes">✓ Yes</span>', '<span class="status no">No</span>'),
      Class = d$hotspot_class,
      Subjects = d$n_subjects, Events = d$n_events,
      Location = cluster_location(d),
      Domain = cluster_domains(dat()$ann, d$cluster_id),
      `Also at site` = vapply(strsplit(d$colocated, ";"), function(t) paste(type_chip(t), collapse = " "), ""),
      `Linear span (aa)` = d$span_aa,
      `3D diam. (Å)` = d$diameter_3d,
      Compaction = d$compaction,
      `p size 1D` = d$p_1d_size, `p size 3D` = d$p_3d_size, `p joint` = d$p_any_joint,
      pLDDT = d$mean_plddt,
      # Remaining GRIN3D cluster-table columns, hidden until picked under "Columns".
      `Sig 1D` = d$significant_1d, `Sig 3D` = d$significant_3d,
      `p joint 1D` = d$p_1d_joint, `p joint 3D` = d$p_3d_joint,
      `p protein 1D` = d$p_1d_protein, `p protein 3D` = d$p_3d_protein,
      `1D span (tree)` = d$diameter_1d, Residues = d$n_residues,
      `Start aa` = d$residue_min, `End aa` = d$residue_max, Trees = d$tree_sources,
      check.names = FALSE)
    hidden <- which(names(out) %in% c("Sig 1D", "Sig 3D", "p joint 1D", "p joint 3D", "p protein 1D",
                                      "p protein 3D", "1D span (tree)", "Residues", "Start aa", "End aa", "Trees")) - 1
    datatable(out, rownames = FALSE, selection = "single", escape = -c(1, 2, 3, which(names(out) == "Also at site")),
              style = "bootstrap5", class = "table table-sm table-hover", extensions = "Buttons",
              options = list(pageLength = 10, scrollX = TRUE, dom = "Bftip",
                             buttons = list(list(extend = "colvis", text = "Columns"),
                                            list(extend = "csv", text = "Download CSV")),
                             language = list(search = "", searchPlaceholder = "Filter clusters"),
                             columnDefs = list(list(className = "cls", targets = 3),
                                               list(visible = FALSE, targets = hidden),
                                               # "36.3×", or a muted dash where compaction does not apply
                                               # (CNA exon clusters, single-residue clusters); sorts numerically.
                                               list(targets = which(names(out) == "Compaction") - 1,
                                                    render = JS("function(d, type) { if (type !== 'display') return d === null ? -1 : d;",
                                                                "return d === null ? '<span title=\"Not applicable: CNA exon cluster or single residue\" style=\"color:#9A9994\">—</span>' : d.toFixed(1) + '×'; }"))))) |>
      formatRound(c("3D diam. (Å)", "1D span (tree)"), 1) |>
      formatSignif(c("p size 1D", "p size 3D", "p joint", "p joint 1D", "p joint 3D",
                     "p protein 1D", "p protein 3D"), 2)
  })

  output$structure <- renderR3dmol({
    clusters <- dat()$clusters
    residues <- dat()$residues
    coords <- dat()$coords
    v <- r3dmol() |>
      m_add_model(data = dat()$pdb_text, format = "pdb")

    if (!is.null(selected())) return(draw_focus(v))

    v <- v |> m_set_style(style = m_style_cartoon(color = backbone_grey, opacity = 0.9))

    mode <- if (is.null(input$backbone)) "hotspots" else input$backbone
    if (mode == "plddt") {
      for (b in plddt_bins) {
        r <- coords$residue[coords$plddt >= as.numeric(b[1]) & coords$plddt < as.numeric(b[2])]
        if (length(r)) v <- v |> m_set_style(sel = m_sel(resi = r),
                                             style = m_style_cartoon(color = b[3]))
      }
    }

    # CNA layers: exons in significant hotspot clusters (or, optionally, every
    # altered exon); each residue takes the visible CNA type with the most events.
    # ponytail: one color per residue; switch to per-layer surfaces if overlap display matters
    if (mode == "domains" && !is.null(dat()$ann)) {
      # Each domain as colored backbone inside a translucent envelope, labeled at its middle.
      d <- dat()$ann$domains
      for (i in seq_len(nrow(d))) {
        r <- seq(d$start[i], d$end[i])
        v <- v |>
          m_set_style(sel = m_sel(resi = r), style = m_style_cartoon(color = d$color[i])) |>
          m_add_surface(type = "MS", atomsel = m_sel(resi = r),
                        style = m_style_surface(opacity = 0.35, colorScheme = NULL, color = d$color[i])) |>
          m_add_label(text = d$feature_name[i], sel = m_sel(resi = r[ceiling(length(r) / 2)], atom = "CA"),
                      style = m_style_label(fontSize = 11, fontColor = ink, backgroundColor = "#FFFFFF",
                                            backgroundOpacity = 0.9, borderColor = d$color[i], borderThickness = 2))
      }
    }
    sig_ids <- clusters$cluster_id[clusters$significant_any]
    if (mode == "hotspots") bg <- residues[residues$alteration_type %in% intersect(layers()$cna, cna_types) &
                     (if (layers()$all_cna) is.na(residues$cluster_id)
                      else residues$cluster_id %in% sig_ids), ]
    if (mode == "hotspots" && nrow(bg)) {
      bg <- bg[order(bg$residue, -bg$n_events), ]
      bg <- bg[!duplicated(bg$residue), ]
      for (t in unique(bg$alteration_type)) {
        v <- v |> m_set_style(sel = m_sel(resi = bg$residue[bg$alteration_type == t]),
                              style = m_style_cartoon(color = layer_colors[[t]]))
      }
    }

    # Mutation layers: C-alpha spheres sized by events of the visible types; a
    # residue hit by both SNVs and indels takes the type with more events.
    mut <- residues[is.na(residues$cluster_id) & residues$alteration_type %in% layers()$mut, ]
    if (nrow(mut)) {
      total <- tapply(mut$n_events, mut$residue, sum)
      mut <- mut[order(mut$residue, -mut$n_events), ]
      mut <- mut[!duplicated(mut$residue), ]
      for (i in seq_len(nrow(mut))) {
        v <- v |> m_add_style(
          sel = m_sel(resi = mut$residue[i], atom = "CA"),
          style = m_style_sphere(colorScheme = NULL, color = layer_colors[[mut$alteration_type[i]]],
                                 opacity = 0.9,
                                 radius = min(3, 0.7 + 0.3 * sqrt(total[[as.character(mut$residue[i])]]))))  # capped so hotspot residues (e.g. FBXW7 R465, 135 events) don't swallow their neighbors
      }
    }

    v |> m_zoom_to()
  })

  # Everything drawn around the selected cluster, shared by the structure, the
  # legend key and the cluster panel.
  focus <- reactive({
    s <- selected()
    if (is.null(s)) return(NULL)
    d <- dat()
    cr <- d$residues[d$residues$cluster_id %in% s$cluster_id, ]
    seed <- unique(cr$residue)
    co <- colocated(d, s, visible())
    dom <- if (is.null(d$ann)) NULL else d$ann$domains
    if (!is.null(dom)) dom <- dom[mapply(function(a, b) any(seed >= a & seed <= b), dom$start, dom$end), , drop = FALSE]
    zone <- if (s$source == "mutation") near_residues(d, seed) else seed
    sites <- d$sites[d$sites$residue %in% zone, ]
    ex <- unique(d$residues[!is.na(d$residues$exon_order), c("residue", "exon_order")])
    # A residue on an exon boundary maps to both exons, so CNA clusters use their own exon list.
    exons <- if (s$source == "cna") as.integer(strsplit(s$exons, ";")[[1]])
             else sort(unique(ex$exon_order[ex$residue %in% seed]))
    list(s = s, cr = cr, seed = seed, co = co, dom = dom, sites = sites, exons = exons)
  })

  # Focus view for one cluster, drawn like a figure: the rest of the protein
  # fades to a pale backbone; the enclosing domain is a translucent envelope;
  # mutated residues are sphere clumps and CNA exons thick backbone, each in its
  # alteration color, including significant clusters of other types at the same
  # site; functional sites nearby are a violet surface.
  draw_focus <- function(v) {
    f <- focus()
    s <- f$s
    d <- dat()
    col <- layer_colors[[s$alteration_type]]
    # Labels keep ink text on a light plate; the colored border carries identity.
    label <- function(v, text, resi, border = col, size = 12, align = "topLeft") m_add_label(v,
      text = text, sel = m_sel(resi = resi, atom = "CA"),
      style = m_style_label(fontSize = size, fontColor = ink, backgroundColor = "#FFFFFF",
                            backgroundOpacity = 0.92, borderColor = border, borderThickness = 2,
                            alignment = align))
    centroid <- colMeans(d$coords[d$coords$residue %in% f$seed, c("x", "y", "z")])
    dist_to_seed <- sqrt(colSums((t(d$coords[, c("x", "y", "z")]) - centroid)^2))
    # Zoom on the functional site when there is one (the pocket is the point),
    # else on the cluster and the other types drawn with it.
    # A CNA cluster covers whole exons, so the site alone frames it; a mutation
    # cluster is small enough to keep in frame beside the site.
    zoom <- if (!nrow(f$sites)) union(f$seed, f$co$residue[f$co$alteration_type %in% mut_types])
            else if (s$source == "cna") unique(f$sites$residue)
            else union(f$seed, f$sites$residue)

    v <- v |> m_set_style(style = m_style_cartoon(color = focus_grey))

    # Domains: backbone in the domain color, labeled at a residue ~16 Å out from
    # the cluster so the label stays in view. No envelope here: the focus view
    # zooms inside the domain, where a surface around the camera fogs the scene.
    # The envelope belongs to the whole-protein "Domains" backbone mode.
    for (i in seq_len(NROW(f$dom))) {
      r <- seq(f$dom$start[i], f$dom$end[i])
      at <- r[which.min(abs(dist_to_seed[match(r, d$coords$residue)] - 16))]
      v <- v |>
        m_set_style(sel = m_sel(resi = r), style = m_style_cartoon(color = f$dom$color[i], opacity = 0.95)) |>
        label(sprintf("%s (aa %d–%d)", f$dom$feature_name[i], f$dom$start[i], f$dom$end[i]), at,
              border = f$dom$color[i], size = 11)
    }

    # Copy number: whole exons as thick backbone, alternating shade so boundaries show.
    x <- f$co[!f$co$alteration_type %in% mut_types, ]
    x$exon_order <- d$residues$exon_order[match(paste(x$cluster_id, x$residue), paste(d$residues$cluster_id, d$residues$residue))]
    cna <- rbind(if (s$source == "cna") f$cr[, c("residue", "exon_order", "alteration_type")],
                 x[, c("residue", "exon_order", "alteration_type")])
    cna <- cna[!duplicated(cna$residue), ]  # a CNA seed keeps its own color on shared exons
    for (t in unique(cna$alteration_type)) {
      tc <- layer_colors[[t]]
      light <- sum(col2rgb(tc) * c(.299, .587, .114)) > 150
      tint <- colorRampPalette(c(tc, if (light) ink else "white"))(4)[2]
      exons <- sort(unique(cna$exon_order[cna$alteration_type == t]))
      for (k in seq_along(exons)) {
        r <- cna$residue[cna$alteration_type == t & cna$exon_order == exons[k]]
        v <- v |>
          m_set_style(sel = m_sel(resi = r), style = m_style_cartoon(color = if (k %% 2) tc else tint, thickness = 1)) |>
          label(paste("Exon", exons[k]), r[ceiling(length(r) / 2)], border = tc)
      }
    }

    # Functional sites: a violet surface over the site residues.
    for (k in unique(f$sites$kind)) {
      r <- unique(f$sites$residue[f$sites$kind == k])
      v <- v |>
        m_add_surface(type = "MS", atomsel = m_sel(resi = r),
                      style = m_style_surface(opacity = 0.55, colorScheme = NULL, color = site_color)) |>
        label(k, r[which.max(dist_to_seed[match(r, d$coords$residue)])], border = site_color, size = 11)
    }

    # Mutations: every atom of each mutated residue as a sphere clump; other
    # types' residues get smaller spheres so the selected cluster stands out. A
    # residue hit by both keeps the selected color on its side chain and takes
    # the other color on its backbone.
    bg <- d$residues[is.na(d$residues$cluster_id), ]
    events <- function(t, r) { e <- bg$n_events[bg$alteration_type == t & bg$residue == r]; if (length(e)) e[1] else 0 }
    muts <- rbind(
      if (s$source == "mutation") data.frame(residue = f$seed, type = s$alteration_type),
      unique(data.frame(residue = f$co$residue, type = f$co$alteration_type)[f$co$alteration_type %in% mut_types, ]))
    muts$events <- vapply(seq_len(nrow(muts)), function(i) as.numeric(events(muts$type[i], muts$residue[i])), 0)
    if (s$source == "mutation") muts$events[seq_along(f$seed)] <- f$cr$n_events[match(f$seed, f$cr$residue)]
    muts$shared <- duplicated(muts$residue)
    for (i in seq_len(nrow(muts))) {
      tc <- layer_colors[[muts$type[i]]]
      sel <- if (muts$shared[i]) m_sel(resi = muts$residue[i], atom = c("N", "CA", "C", "O")) else m_sel(resi = muts$residue[i])
      v <- v |>
        m_set_style(sel = m_sel(resi = muts$residue[i]), style = m_style_cartoon(color = tc, thickness = 1)) |>
        m_add_style(sel = sel, style = m_style_sphere(colorScheme = NULL, color = tc,
                                                      scale = if (muts$type[i] == s$alteration_type) 0.85 else 0.5))
    }
    # ponytail: glycines have no side chain, so a shared glycine shows only the second color
    # Labels: the selected cluster's busiest residues first, then a few of the other types'.
    own <- muts$type == s$alteration_type & s$source == "mutation"
    top <- rbind(head(muts[own, ][order(-muts$events[own]), ], 8), head(muts[!own, ][order(-muts$events[!own]), ], 3))
    top <- top[!duplicated(top$residue), ]
    for (i in seq_len(nrow(top))) {
      prefix <- if (top$type[i] == s$alteration_type) "" else paste0(layer_labels[[top$type[i]]], " ")
      v <- v |> label(sprintf("%s%d (%d)", prefix, top$residue[i], top$events[i]), top$residue[i],
                      border = layer_colors[[top$type[i]]])
    }

    v |> m_zoom_to(sel = m_sel(resi = zoom, expand = 8))
  }

  # The tool opens on the top convergent cluster (PTEN: the HOMDEL exon 1/2/5
  # site), so the first thing on screen is a worked example rather than an empty
  # structure. Only on load; after that the table drives the selection.
  observeEvent(input$cluster_table_rows_current, once = TRUE,
    if (length(input$cluster_table_rows_current)) selectRows(dataTableProxy("cluster_table"), 1))

  # A protein with no convergent clusters (no significant clusters of two types
  # at one site) falls back to All rather than showing an empty table.
  observeEvent(dat(), {
    if (!identical(input$view, "conv")) return()
    cl <- dat()$clusters[dat()$clusters$significant_3d, ]
    any_conv <- any(vapply(seq_len(nrow(cl)),
                           function(i) nrow(colocated(dat(), cl[i, ], dat()$info$types[[1]])) > 0, TRUE))
    if (!any_conv) updateRadioButtons(session, "view", selected = "all")
  })

  observeEvent(input$clear_selection, {
    selected_id(NULL)
    selectRows(dataTableProxy("cluster_table"), NULL)
  })
  output$clear_button <- renderUI(if (!is.null(selected())) tagList(
    tags$button(class = "btn btn-sm btn-outline-secondary", onclick = "grin3dExport()",
                title = "Download this view as an SVG figure with a CSV of the cluster's statistics", "Export"),
    actionButton("clear_selection", "Show all layers", class = "btn-sm btn-outline-secondary")))

  # Export: the browser sends the 3D snapshot, then the zip download is triggered.
  snapshot <- reactiveVal(NULL)
  observeEvent(input$structure_png, {
    snapshot(input$structure_png)
    session$sendCustomMessage("grin3d-download", "export_zip")
  })
  output$export_zip <- downloadHandler(
    filename = function() sprintf("%s_%s.zip", dat()$info$protein, selected()$cluster_id),
    content = function(file) {
      f <- req(focus())
      loc <- cluster_location(f$s)
      stem <- file.path(tempdir(), paste(dat()$info$protein, f$s$cluster_id, sep = "_"))
      fm <- feature_map_svg(nrow(dat()$coords), dat()$ann, dat()$residues,
                            intersect(layers()$mut, mut_types), layer_colors, f$s)
      writeLines(cluster_figure_svg(dat(), f, req(snapshot()), fm, loc), paste0(stem, ".svg"))
      write.csv(cluster_stats_row(dat(), f, loc), paste0(stem, "_stats.csv"), row.names = FALSE)
      utils::zip(file, paste0(stem, c(".svg", "_stats.csv")), flags = "-jq")
    },
    contentType = "application/zip")
  outputOptions(output, "export_zip", suspendWhenHidden = FALSE)  # the link is hidden, but needs its URL

  key_row <- function(color, text) div(class = "key-row", tags$i(class = "sw env", style = sprintf("background:%s", color)), span(text))

  # Legend key for what the structure shows beyond the layer toggles.
  output$context_key <- renderUI({
    f <- focus()
    if (!is.null(f)) {
      types <- unique(c(f$s$alteration_type, f$co$alteration_type))
      return(div(class = "legend-section focus-key",
        div(class = "eyebrow", "At this site"),
        lapply(types, function(t) div(class = "key-row", HTML(type_chip(t)),
          span(class = "count", if (t == f$s$alteration_type) "selected"
               else if (f$s$source == "cna" && !t %in% mut_types) "within exons"))),
        lapply(seq_len(NROW(f$dom)), function(i) key_row(f$dom$color[i], f$dom$feature_name[i])),
        lapply(unique(f$sites$kind), function(k) key_row(site_color, k))))
    }
    if (identical(input$backbone, "domains") && !is.null(dat()$ann)) {
      d <- dat()$ann$domains
      d <- d[!duplicated(d$group), ]
      div(class = "legend-section focus-key", div(class = "eyebrow", "Domains"),
          lapply(seq_len(nrow(d)), function(i) key_row(d$color[i], d$group[i])))
    }
  })

  # Cluster panel: other significant clusters and functional sites at the same 3D site.
  site_section <- function(f) {
    if (!nrow(f$co) && !nrow(f$sites)) return(NULL)
    mut_seed <- f$s$source == "mutation"
    tagList(
      div(class = "eyebrow", style = "margin-top:14px", "At this 3D site"),
      lapply(split(f$co, f$co$cluster_id), function(x) div(class = "ann-item",
        HTML(type_chip(x$alteration_type[1])), span(class = "mono", x$cluster_id[1]),
        span(class = "res", compact_residues(paste(x$residue, collapse = ";"))))),
      lapply(split(f$sites, f$sites$name), function(x) div(class = "ann-item",
        span(class = "dom-sw", style = sprintf("background:%s", site_color)), span(x$name[1]),
        span(class = "res", compact_residues(paste(unique(x$residue), collapse = ";"))))),
      p(class = "scope-note",
        if (mut_seed) sprintf("Significant clusters of other types with residues within %g Å of this cluster (CNA clusters: exons containing it), and functional sites within %g Å.", near_cutoff, near_cutoff)
        else "Significant clusters of other types lying inside this cluster's exons, and functional sites in those exons.",
        "Each cluster is significant on its own; the combined site is not a separate test.")
    )
  }

  output$inspector <- renderUI({
    s <- selected()
    if (is.null(s)) {
      cl <- dat()$clusters
      sig <- table(factor(cl$alteration_type[cl$significant_any], levels = unique(cl$alteration_type)))
      return(div(class = "empty",
        h5("No cluster selected"),
        p("Select a row in the cluster table to focus the structure on that cluster and see its evidence here."),
        div(class = "eyebrow", style = "margin-top:14px", "Significant clusters"),
        tags$ul(class = "sig-list", lapply(names(sig), function(t)
          tags$li(HTML(type_chip(t)), tags$b(sig[[t]])))),
        chembl_section(dat()$ann))
      )
    }
    r <- dat()$residues
    r <- r[r$cluster_id %in% s$cluster_id, ]
    fmt_p <- function(x) formatC(x, digits = 2, format = "g")
    tagList(
      div(style = "display:flex; align-items:center; gap:10px; flex-wrap:wrap",
          HTML(type_chip(s$alteration_type)),
          if (s$significant_any) span(class = "status yes", "✓ Significant") else span(class = "status no", "Not significant")),
      div(class = "mono", style = "margin-top:8px; font-size:15px", s$cluster_id),
      div(style = "color:#6B6A66", s$hotspot_class),
      div(class = "kpis",
        div(class = "kpi", div(class = "v", s$n_subjects), div(class = "l", "Subjects")),
        div(class = "kpi", div(class = "v", s$n_events), div(class = "l", "Events")),
        div(class = "kpi", div(class = "v", round(s$diameter_3d, 1)), div(class = "l", "3D diameter (Å)"))),
      tags$dl(class = "kv",
        tags$dt("Location"), tags$dd(cluster_location(s)),
        tags$dt(if (s$source == "cna") "1D span (exons)" else "1D span (aa)"), tags$dd(round(s$diameter_1d, 1)),
        tags$dt("Residues"), tags$dd(sprintf("%d (aa %d–%d)", s$n_residues, s$residue_min, s$residue_max)),
        tags$dt("Compaction"), tags$dd(if (s$source == "cna") "n/a for exon clusters" else if (s$diameter_3d > 0) sprintf("%.1f× (%d aa span, %.1f Å across)",
          (s$residue_max - s$residue_min) * 3.8 / s$diameter_3d, s$residue_max - s$residue_min, s$diameter_3d) else "single position"),
        tags$dt("Significant in"), tags$dd(c(if (s$significant_1d) "sequence (1D)", if (s$significant_3d) "structure (3D)", "neither")[1:max(1, s$significant_1d + s$significant_3d)] |> paste(collapse = " and ")),
        tags$dt("Mean pLDDT"), tags$dd(s$mean_plddt,
          span(style = "color:#6B6A66", sprintf(" · %d residues < 70", length(unique(r$residue[r$plddt < 70])))))),
      div(class = "eyebrow", style = "margin-top:14px", "Empirical p-values"),
      div(class = "pgrid",
        span(), span(class = "h", "1D"), span(class = "h", "3D"),
        span("Size-specific"), span(fmt_p(s$p_1d_size)), span(fmt_p(s$p_3d_size)),
        span("Joint"), span(fmt_p(s$p_1d_joint)), span(fmt_p(s$p_3d_joint)),
        tags$b("Joint, either"), span(), tags$b(fmt_p(s$p_any_joint))),
      site_section(focus()),
      div(class = "eyebrow", style = "margin-top:14px", "Contributing events"),
      DTOutput("member_table")
    )
  })

  # Annotations for the selected cluster, or the whole protein when none is selected.
  # server = FALSE: these tables change their column set (a selected cluster shows
  # "Overlapping residues", the whole protein shows "Position"; only CNA members
  # have a coding fraction). With server-side processing the browser keeps asking
  # for the previous table's columns and DataTables raises "column name ... is not
  # found in data". Both tables are small, so the client holds all the rows.
  output$annotation_table <- renderDT(server = FALSE, {
    d <- annotation_rows(dat()$ann, selected()$cluster_id)
    if (is.null(d) || !nrow(d)) return(NULL)
    datatable(d, rownames = FALSE, selection = "none", escape = -2, style = "bootstrap5",
              class = "table table-sm table-hover",
              options = list(pageLength = 10, dom = "ftip", scrollX = TRUE, order = list(),
                             language = list(search = "", searchPlaceholder = "Filter annotations")))
  })

  output$annotation_sub <- renderUI({
    ann <- dat()$ann
    if (is.null(ann)) return("no annotations for this protein")
    s <- selected()
    n <- nrow(annotation_rows(ann, s$cluster_id))
    if (is.null(s)) sprintf("%d features on %s · select a cluster to see the ones it overlaps",
                            n, dat()$info$protein)
    else sprintf("%d features overlapping %s", n, s$cluster_id)
  })

  output$annotation_note <- renderUI({
    ann <- dat()$ann
    if (is.null(ann)) return(p(class = "hint", "This protein has no annotation files in data/annotations/."))
    s <- selected()
    if (!is.null(s) && !nrow(annotation_rows(ann, s$cluster_id)))
      return(p(class = "hint", "No annotated features overlap this cluster."))
    p(class = "scope-note", paste0(
      "UniProt, InterPro and PDBe-KB features ",
      if (is.null(s)) "for the whole protein"
      else if (s$source == "cna") "overlapping the residues this cluster's exons encode"
      else "overlapping the cluster's residues",
      ". Annotations add context only; they do not change significance."))
  })

  output$member_table <- renderDT(server = FALSE, {
    s <- req(selected())
    members <- dat()$members
    d <- members[members$cluster_id == s$cluster_id, c("subject_id", "event_type", "detail", "coding_fraction")]
    names(d) <- c("Subject", "Type", "Detail", "Coding fraction")
    # CNA events only; a bar on a fixed 0-1 track compares at a glance, which the
    # raw number does not. Mutation events have no coding fraction.
    frac <- which(names(d) == "Coding fraction") - 1
    if (all(is.na(d$`Coding fraction`))) d$`Coding fraction` <- NULL
    datatable(d, rownames = FALSE, selection = "none", style = "bootstrap5",
              class = "table table-sm",
              options = list(pageLength = 6, dom = "tip", scrollX = TRUE,
                             columnDefs = if (!is.null(d$`Coding fraction`)) list(list(
                               targets = frac, className = "frac-cell",
                               render = JS("function(d, type) { if (type !== 'display') return d;",
                                           "if (d === null) return '';",
                                           sprintf("return '<span class=\"frac\" title=\"coding fraction ' + d.toFixed(2) + ' of 1\"><i style=\"width:' + Math.max(d * 100, 1.5).toFixed(1) + '%%; background:%s\"></i></span>'; }",
                                                   layer_colors[[s$alteration_type]]))))))
  })
}

shinyApp(ui, server)
