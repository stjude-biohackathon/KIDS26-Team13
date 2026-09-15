# GRIN3D BioHackathon Project

GRIN3D is a proof-of-concept framework for identifying mutation and copy-number alteration (CNA) hotspots in three-dimensional protein structures.

The project asks whether genomic lesions that are distant in the linear protein sequence can converge into significant clusters in 3D space.

## Project Profile

- **Project name:** GRIN3D: From Linear Sequence to 3D Protein Hotspots  
- **Main question:** Can we identify statistically significant clusters of mutations and CNAs in protein structures?
- **Data and inputs:** Mutation and CNA data, coding-exon annotations, transcript-to-protein mappings, and AlphaFold protein structures.
- **Example proteins:** TP53, PTEN, and CDKN2A.
- **Expected output:** Reproducible analyses, evaluated null models, documented limitations, and integrated protein visualization.
- **Tools and stack:** R, AlphaFold, Ensembl, UniProt, `bio3d`, `ensembldb`, `AnnotationHub`, `GenomicRanges`, `r3dmol`, and R Shiny.
- **Team lead:** Abdelrahman Elsayed ([@abdel-elsayed87](https://github.com/abdel-elsayed87))
- **Team members and roles:** See [project-management/team.md](project-management/team.md).
- **Communication:** Use GitHub Issues and pull requests for code-related discussions.

## Vision and Mission

- **Vision:** Establish a unified framework for identifying and interpreting different genomic alteration types in protein 3D space.
- **Mission:** Develop, test, and integrate statistical, structural, and visualization methods for detecting mutation and CNA hotspots.

## About

Current protein 3D hotspot tools primarily focus on mutations. Other alterations, particularly CNAs, are difficult to evaluate because they affect genomic intervals and often cannot be assigned to individual amino acids.

GRIN3D evaluates different alteration types at appropriate levels of structural resolution. Mutations are mapped to individual residues, while CNAs are mapped to their affected coding exons.

The proposed framework tests whether mutations or CNA-affected exons form spatial clusters that are more compact than expected by chance. It then evaluates their biological and clinical relevance by mapping them to protein domains, protein–protein interaction interfaces, catalytic sites, post-translational modification sites, and drug-binding regions. The simulations account for subjects, mutation burden, CNA type, event size, and structural eligibility.

## Project Components

### Mutation-Hotspot Analysis

The mutation workflow:

- Maps mutation events to amino-acid residues.
- Builds a 1D sequence tree and a 3D structural tree.
- Combines clusters identified by both trees.
- Simulates mutation positions while preserving subject IDs, event IDs, mutations per subject, and event lengths.
- Rebuilds both trees in every simulation.
- Reports size-specific, protein-wide, and joint 1D–3D empirical p-values.

TP53 is included as an example with a dense mutation pattern. Additional proteins will be used to test sparse mutation patterns.

### CNA-Hotspot Analysis

The CNA workflow:

- Preserves each CNA as one intact genomic interval.
- Maps each CNA to the coding exons it overlaps.
- Separates all coding overlaps from partial, structurally localizable events.
- Retains broad and whole-gene CNAs as gene-level evidence.
- Builds a 1D exon-order tree and a 3D exon-proximity tree.
- Analyzes HOMDEL, HETDEL, GAIN, and AMP events separately and as a combined CNA group.
- Simulates CNA placement while preserving subjects, CNA types, exon spans, and approximate coding-residue spans.

The current 3D exon-distance metric is the mean of the five smallest reliable cross-exon Cα distances. Alternative settings, including 10 residue pairs, will be tested.

PTEN is used to test multi-exon structural convergence. CDKN2A is used to examine limitations in proteins with few coding exons.

### Integrated Visualization

The visualization component will:

- Display mutation and CNA hotspots on the same protein.
- Color-code mutations, HOMDELs, HETDELs, gains, and amplifications.
- Allow users to show or hide different alteration layers.
- Connect each cluster to its statistical and structural evidence.
- Add domains, motifs, interaction sites, and other functional annotations when available.

## Repository Contents

- **`development-code/`** — GRIN3D mutation, CNA, and input-preparation scripts.
- **`examples/`** — Example runners, input files, and results for TP53, PTEN, and CDKN2A.
- **`datasets/`** — Public mutation and CNA datasets for additional testing.
- **`docs/`** — GitHub guidance and troubleshooting information.
- **`project-management/`** — Team information, project plans, meeting notes, and checklists.

## Getting Started

1. Clone or download the repository.
2. Review the scripts in `development-code/`.
3. Select an analysis from the `examples/` folder.
4. Open the example runner and update `analysis.dir`.
5. Confirm that the required files are in `input_files/`.
6. Start with a small number of simulations to test the workflow.
7. Use at least 10,000 simulations for more stable comparisons.
8. Review the output files in the example's `results/` folder.

## Roadmap and Milestones

| When | Focus | Expected outcome |
| --- | --- | --- |
| Day 1 | Reproduce the examples and review the statistical logic | Confirmed inputs, outputs, null-model assumptions, and assigned tasks |
| Day 2 | Test statistical, structural, and visualization improvements | Updated methods and a working visualization draft |
| Day 3 | Integrate results and document limitations | Reproducible demonstration, summary slides, and clear next steps |

## BioHackathon Priorities

### 1. Evaluate the Statistical Framework

- Confirm what each null model preserves.
- Verify mutation and CNA simulations.
- Review size-specific, protein-wide, and joint p-values.
- Identify events with insufficient alternative placements.

### 2. Improve the Structural Framework

- Evaluate alternative distance and linkage methods.
- Compare CNA exon distances based on the 5 versus 10 closest Cα pairs.
- Use pLDDT, PAE, and IDR annotations to assess structural reliability.
- Retain 1D evidence when reliable 3D analysis is not possible.

### 3. Test Robustness

- Compare dense and sparse mutation patterns.
- Evaluate PTEN, CDKN2A, TP53, and additional proteins.
- Test whether clusters remain stable across reasonable settings.

### 4. Build an Integrated Visualization

- Display mutation and CNA hotspots together.
- Color-code different alteration types.
- Connect each cluster to its statistical and biological evidence.

## Interpretation and Limitations

GRIN3D is development code and should be treated as a proof of concept.

Important limitations include:

- AlphaFold confidence varies across protein regions.
- PAE may indicate uncertainty in the relative positions of protein domains.
- Intrinsically disordered regions may not support reliable 3D conclusions.
- CNAs provide exon-level rather than residue-specific localization.
- Sparse lesion patterns may provide limited statistical power.
- Transcript, exon, protein, and UniProt mappings must be consistent.
- CNA coverage p-values are conditional on a CNA already overlapping the selected gene. They do not test genome-wide gene recurrence.

The BioHackathon goal is to demonstrate feasibility, identify limitations, and define the next development steps. The current framework is not fully validated or intended for clinical use.

## Data

The `datasets/` folder contains public mutation and CNA data from a published T-ALL study.

## Contributing

1. Create a branch from the current `main` branch.
2. Make changes related to one defined task.
3. Test the affected example runner.
4. Describe the changes in the commit and pull request.
5. Submit a pull request for team review.
6. Merge the pull request after approval.

Keep experimental results, assumptions, and known limitations clearly documented so other team members can reproduce the work.

This code is under active development. Please contact the team lead before redistributing it, presenting it as completed work, or using it outside this project. Do not claim ownership of the code or project materials.


