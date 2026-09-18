# GRIN3D biological annotation workflow

## Purpose and statistical boundary

This R workflow adds biological context after GRIN3D mutation-hotspot testing.
It does not change cluster membership, p-values, significance flags, or hotspot
classifications. Annotation availability therefore cannot make a candidate
statistically significant.

The reusable implementation is
`development-code/GRIN3D_annotate_protein_clusters.R`. It accepts any GRIN3D
mutation-cluster CSV plus the UniProt accession for the exact analyzed protein
isoform. The TP53 runner is
`examples/TP53_mutations/run_TP53_annotation_example.R`.

## Environment and command

Create or update the repository environment and run the example from the
repository root:

```bash
conda env update -n kids26-team13 -f environment.yml
conda activate kids26-team13
Rscript examples/TP53_mutations/run_TP53_annotation_example.R
```

The only non-base R dependency is `jsonlite`. Base R creates the CSV and
single-file HTML dashboard outputs.

For another protein, source the module and change the paths and accession:

```r
source("development-code/GRIN3D_annotate_protein_clusters.R")

annotation.results <- annotate_grin3d_clusters(
  cluster.file = "path/to/results/mutation_hotspot_clusters.csv",
  uniprot.accession = "P12345",
  output.dir = "path/to/annotations",
  protein = "GENE",
  sources = c("uniprot", "interpro", "pdbe", "chembl")
)

write_cluster_annotation_dashboard(annotation.results)
```

Responses are cached under `OUTPUT_DIR/raw`. Set `offline = TRUE` to require
cached responses or `refresh = TRUE` to deliberately retrieve current records.
Keep `annotation_sources.csv` and the raw JSON with reported results because
source databases change over time.

## How clusters enter the workflow

The unchanged `mutation_hotspot_clusters.csv` from
`run_grin3d_mutation_hotspots()` is the input. It must contain:

- `cluster_id`: the unique GRIN3D candidate identifier;
- `residues`: semicolon- or comma-separated, one-based protein positions.

The module retains the statistical fields—including candidate/tree source,
event and subject counts, p-values, and significance flags—but never uses them
to choose annotations. It rejects duplicate cluster IDs, invalid residue lists,
multiple protein labels, and residues beyond the retrieved UniProt sequence.

### Exon-based CNA candidates

`grin3d_prepare_cna_annotation_input()` bridges a
`CNA_exon_hotspot_clusters.csv` file to the residue-coordinate annotation
schema. For each candidate, it looks up every `exon_id` in the selected
transcript's exon-to-protein map, expands `aa_start:aa_end`, takes the union, and
deduplicates residues shared across splice boundaries. It writes:

- `CNA_annotation_clusters.csv`, the annotation input retaining the CNA
  cluster's statistics; and
- `CNA_cluster_exon_residue_mapping.csv`, the auditable exon-to-residue bridge.

An annotation overlap in this setting means that a feature intersects a region
encoded by a candidate exon. It does not mean that every residue in that exon
was altered, nor does it convert a CNA into a point mutation. The PTEN example
runner is `examples/PTEN_CNA/run_PTEN_CNA_annotation_example.R`.

## Sources and evidence

### UniProtKB REST

The endpoint is `https://rest.uniprot.org/uniprotkb/{ACCESSION}.json`.
The normalized catalog includes domains, regions, motifs, binding and active
sites, PTMs, interaction regions, repeats, and selected structural regions.
Natural variants, mutagenesis records, secondary structure, and sequence
conflicts are excluded because they answer different questions and would
overwhelm the cluster interpretation.

Feature evidence codes and identifiers are preserved. Publications attached to
a feature are recorded as `feature_evidence` with PubMed/DOI metadata when the
UniProt record supplies it.

### InterPro REST

Protein matches come from
`https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/{ACCESSION}/`.
InterPro contributes protein families, domains, homologous superfamilies, and
conserved sites. Each matching entry is also retrieved so its description and
literature are retained.

InterPro papers are labelled `entry_literature`: they describe the database
entry and are not direct experimental evidence for a particular GRIN3D cluster.
UniProt and InterPro boundaries remain separate source records so agreement or
disagreement stays visible.

### PDBe-KB graph API and SIFTS

The adapter queries `uniprot/interface_residues/{ACCESSION}` for experimentally
observed macromolecular interfaces and `uniprot/ligand_sites/{ACCESSION}` for
experimentally observed ligand contacts. SIFTS supplies the mapping from PDB
residues to the analyzed UniProt sequence.

The normalized records retain the partner or PDB chemical-component identifier,
exact UniProt contact positions, and supporting PDB IDs. Interfaces are grouped
as protein-protein, self, nucleic-acid, or other macromolecular interfaces.
Non-solvent ligand contacts are grouped as cofactor, explicitly drug-linked, or
other small-molecule sites. “Non-solvent ligand” does not by itself mean “drug”;
a drug label requires an explicit ChEMBL or DrugBank cross-reference.

These are experimental structure contacts aggregated across structures reported
by PDBe-KB. They are not AlphaFold proximity predictions and do not change the
GRIN3D 3D hotspot statistic.

### ChEMBL REST

The adapter finds ChEMBL targets whose components include the supplied UniProt
accession, then retrieves drug-mechanism records and molecule names/development
phases. A target can be the single protein or a multi-protein complex.

ChEMBL rows are deliberately exported as target-level context and are never
converted into residue features. Only a separately observed PDBe-KB ligand
contact supports a residue-level structural-contact claim. The claims remain in
different tables unless the source records provide an explicit shared ID.

## Coordinate mapping

Every normalized feature has one or more inclusive, one-based UniProt sequence
segments. For cluster positions `C` and annotation positions `A`, the mapper
calculates:

```text
overlap = C intersect A
```

No row is produced when the intersection is empty. A single-position curated
feature is labelled `exact_residue_overlap`; a regional feature is labelled
`interval_overlap`. PDBe-KB contacts use the specific labels
`experimental_interface_contact` and `experimental_ligand_contact`. Each
individual residue-feature match is also exported.

This workflow does not call a residue “annotated” merely because it is nearby
in sequence or in three-dimensional space. A PDBe-KB row requires that the
exact UniProt-indexed residue was reported as a structure-observed contact.

## Classification hierarchy

The workflow uses three explicit levels while preserving the source record:

```text
annotation_class -> feature_subcategory -> source-native type and name
```

Examples include `site -> metal_binding_site -> UniProt Zn(2+)`,
`modification -> phosphorylation -> UniProt Phosphoserine`, and
`domain -> dna_binding_domain -> InterPro p53 DNA-binding domain`.
`classification_rule` records the versioned deterministic rule used. The
hierarchy organizes source evidence; it does not replace or reinterpret the
original `feature_type`, `feature_name`, description, or accession.

## CSV outputs

- `protein_features.csv`: normalized feature catalog with source coordinates,
  evidence, cross-references, source version/date, and retrieval time.
- `annotation_classification_summary.csv`: the broad classes, normalized
  subclasses, versioned rules, source-native feature types, sources, and counts
  present for the analyzed protein.
- `feature_publications.csv`: feature-publication relationships. The
  `evidence_scope` column distinguishes direct feature evidence from general
  entry literature.
- `cluster_feature_overlaps.csv`: one row per cluster-feature intersection,
  retaining the original statistical fields and overlap fractions.
- `cluster_residue_annotations.csv`: one row per cluster residue per feature;
  this is the authoritative explanation of each residue-level assignment.
- `cluster_annotation_summary.csv`: compact one-row-per-cluster review table.
- `annotation_sources.csv`: API endpoints, source record IDs, versions/dates,
  retrieval timestamps, and cached raw-response paths.
- `pdbe_ligand_evidence.csv`: non-solvent PDB ligands, exact contact residues,
  supporting structures, and explicit ChEMBL/DrugBank/cofactor identifiers.
- `chembl_target_evidence.csv`: single-protein and complex ChEMBL targets that
  contain the analyzed UniProt accession.
- `chembl_drug_mechanisms.csv`: ChEMBL molecule, mechanism, action, and maximum
  development-phase fields, labelled as target-level rather than residue-level.

Tables are normalized deliberately. `feature_id` links an overlap or residue to
its complete feature record and its literature without duplicating large
publication metadata into every cluster row.

## Single-file cluster annotation dashboard

`write_cluster_annotation_dashboard()` creates one self-contained file:
`annotations/cluster_annotation_dashboard.html`. It is generated by R, embeds
all data and JavaScript, and can be opened locally without R, Shiny, or a web
server.

The dashboard lets a reviewer search by cluster or residue, filter by
statistical status, hotspot class, annotation class, or annotation subclass,
and move to the previous or next cluster.
For the selected cluster it contains:

- the full protein sequence axis and all cluster residues;
- a separate, dynamically drawn track for each overlapping feature;
- regional feature boundaries and highlighted overlapping residues;
- relationship type (curated exact/interval overlap or an experimental PDBe-KB
  interface/ligand contact);
- the broad class, normalized subclass, and source-native feature;
- clickable UniProt, InterPro, or PDBe-KB source records;
- evidence codes and up to five direct literature links per feature.
- a separate ChEMBL mechanism table clearly labelled as protein/complex-level
  context rather than an annotation of the selected cluster.

This is a review interface over the CSVs, not another analysis. The dashboard
repeats the statistical/annotation boundary and labels InterPro literature
appropriately. Shiny is R's closest analogue to Dash, but the static dashboard
is preferable here because it is one portable artifact and requires no running
server.

## Interpretation rules

1. Statistical significance comes only from the hotspot result table.
2. Exact site overlap, interval overlap, experimental contact, and predicted 3D
   proximity are distinct.
3. A protein interaction or broad interaction region is not automatically an
   experimentally resolved residue-level interface.
4. InterPro entry literature is background for a feature, not proof of a
   cluster’s functional effect.
5. Reports should retain the accession, evidence code, source/version, and
   retrieval date.
6. Licensed resources such as PhosphoSitePlus or ELM should be optional
   adapters and must follow their redistribution terms.

## Adapter boundary and possible extensions

The four core adapters are UniProtKB, InterPro, PDBe-KB/SIFTS, and ChEMBL. The
normalized feature schema separates retrieval from mapping, so an optional
adapter only needs to emit the same coordinates, category, evidence, source,
and publication fields. IntAct is a reasonable later addition for interaction
evidence, but protein-level interaction evidence and residue-level interfaces
must remain distinct.
