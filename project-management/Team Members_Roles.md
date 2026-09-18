# Team and Roles

- **Team name:** GRIN3D / KIDS26 Team 13
- **Team lead:** Abdelrahman Elsayed ([@abdel-elsayed87](https://github.com/abdel-elsayed87))
- **Communication channel:** GitHub Issues and pull requests
- **Project question/problem:** Can genomic mutations and copy-number alterations (CNAs) that are distant in protein sequence converge into meaningful three-dimensional protein clusters?
- **Expected output:** Working mutation and CNA hotspot prototypes, a visualization example, and documented methods, results, and limitations to guide further development.
- **Tools and stack:** R, AlphaFold structures, Ensembl annotations, public mutation and CNA data, and GitHub.

## Roles

| Person | Role | Main responsibility | Backup or support needed |
| --- | --- | --- | --- |
| Abdelrahman Elsayed | Team lead and integration | Prepare the starting prototype code, assign tasks, coordinate the team, test the workflows across approximately 150 proteins, and gather recommendations for a future production-ready tool | Team review of code, results, and final recommendations |
| Yonghui | Statistical framework | Review the mutation and CNA null models, simulations, and p-value calculations | Input from prototype developers on the intended behavior of each test |
| Candice | CNA hotspot framework | Refine CNA eligibility and localization rules; test the framework using PTEN and CDKN2A | Statistical and structural review of proposed rules |
| Pramesh | Mutation hotspots and interpretation | Test the mutation prototype and explore functional annotation of candidate clusters | Input on protein annotations and suitable comparison proteins |
| Eliijah | Tree construction and cluster selection | Review the mutation and CNA trees, cluster selection, and alternatives | Statistical and structural input on distance definitions |
| Andrew | Integrated visualization | Develop an interactive view of mutation hotspots and CNA-affected exons | Standard result tables and a suitable PTEN example |

## Team Responsibilities and Deliverables

### 1. Validate the Statistical Framework

**Lead: Yonghui**

- Review the mutation and CNA null models separately.
- Confirm what each simulation preserves, including subjects, mutation burden, CNA type, event size, and structural eligibility.
- Verify the empirical p-value calculations and the split between calibration and testing simulations.
- Review size-specific, protein-wide, and joint 1D–3D corrections.
- Identify tests with too few alternative placements or too little variation in the simulated results.
- Document assumptions, potential problems, and unresolved questions.

**Deliverables:** A concise description of each null model and p-value calculation, with confirmed behavior, identified concerns, and prioritized recommendations.

### 2. Refine and Test the CNA Hotspot Framework

**Lead: Candice**

- Define which CNAs are eligible for exon-level and 3D analysis.
- Distinguish partially localizable, broad, and whole-gene CNAs.
- Evaluate coding coverage, exon overlap, and structural eligibility thresholds.
- Test whether results remain stable across CNA types and reasonable analysis settings.
- Use PTEN and CDKN2A as initial contrasting examples.
- Document events that cannot be reliably localized or tested.

**Deliverables:** Proposed CNA eligibility and localization criteria, initial sensitivity results, and limitations that need further investigation.

### 3. Test and Interpret Mutation Hotspots

**Lead: Pramesh**

- Test the mutation prototype on TP53 and proteins with sparser mutation patterns.
- Check whether detected clusters remain stable across reasonable analysis settings.
- Map candidate clusters to protein domains, motifs, catalytic sites, post-translational modification sites, protein–protein interaction interfaces, and drug-binding regions.
- Identify reliable sources and tools for these annotations.
- Keep biological annotations separate from statistical significance.

**Deliverables:** A documented evaluation of the mutation prototype, examples of annotated clusters, and recommendations for adding functional annotations to GRIN3D.

### 4. Evaluate Tree Construction and Cluster Selection

**Lead: Eliijah**

- Review how the mutation prototype builds amino-acid 1D and 3D trees.
- Review how the CNA prototype builds exon-order and exon-proximity trees.
- Evaluate complete linkage, branch selection, and duplicate-cluster handling.
- Identify alternative clustering and distance methods worth testing.
- Compare selected alternatives when feasible.

**Deliverables:** An assessment of the current approach, its limitations, and prioritized improvements, supported by a small sensitivity comparison where feasible.

### 5. Develop Integrated Mutation–CNA Visualization

**Lead: Andrew**

- Show mutation hotspots and CNA-affected exon regions on the same protein structure.
- Use consistent colors for mutations, HOMDELs, HETDELs, gains, and amplifications.
- Let users show or hide individual alteration layers.
- Highlight a selected cluster and its contributing residues or exons.
- Show cluster details, including subjects, events, p-values, and structural confidence.
- Define a standard input table connecting analysis results to protein residues.
- Consider 3Dmol with R Shiny as a starting option, and assess alternatives if needed.

**Deliverables:** A reproducible interactive viewer for PTEN or another suitable protein, its source code, and a defined input format for results from both prototypes.

### 6. Coordinate and Integrate the Project

**Lead: Abdelrahman Elsayed**

- Prepare and share the starting mutation and CNA prototype code.
- Assign tasks and follow up with team members on progress and questions.
- Run the workflows on a larger test dataset of approximately 150 proteins, recording successes, failures, and limitations.
- Review and integrate contributions from the team.
- Gather the team's findings and recommendations into a plan for a future production-ready GRIN3D tool.

**Deliverables:** Shared starting code, a multi-protein test summary, and a consolidated set of recommendations and next steps.
