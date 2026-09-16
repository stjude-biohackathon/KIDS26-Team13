# Team and Roles

- **Team name:** GRIN3D / KIDS26 Team 13
- **Team lead:** Abdelrahman Elsayed ([@abdel-elsayed87](https://github.com/abdel-elsayed87))
- **Communication channel:** GitHub Issues and pull requests
- **Project question/problem:** Can genomic mutations and copy-number alterations (CNAs) that are distant in protein sequence converge into meaningful three-dimensional protein clusters?
- **Expected output:** A reproducible framework for mutation and CNA hotspot analyses, with documented null models, structural eligibility criteria, results, and limitations.
- **Tools and stack:** R, AlphaFold structures, Ensembl annotations, public mutation/CNA data, and GitHub.

## Roles

| Person | Role | Main responsibility | Backup or support needed |
| --- | --- | --- | --- |
| Abdelrahman Elsayed | Team lead | Overall project direction, integration, and coordination | Team review of methods and deliverables |
| Ramzi Alsallaq | CNA analysis and structural eligibility | Developed the CNA eligibility/localization scan; evaluated PTEN and CDKN2A CNA events across coding coverage, exon mapping, structural eligibility, and sensitivity thresholds; produced the analysis result tables and documented untestable whole-gene events | Review of eligibility assumptions and integration with the permutation workflow |

## Current CNA Analysis Contribution

The CNA eligibility scan separates events that can support exon-level and 3D hotspot testing from events that are only interpretable as broad or whole-gene evidence. The current results include:

- Per-event coding coverage and eligibility classifications for PTEN and CDKN2A.
- Event-to-exon mapping and exon structural-eligibility tables.
- A sensitivity analysis across pLDDT, exon-overlap, and affected-exon eligibility thresholds.
- A cross-gene overview of testable and whole-gene-only CNA events.
- A table of untestable events with the reason for exclusion from permutation testing.

Result files are in `02_results/`.

### PTEN CNA sensitivity analyses

Multiple sensitivity analyses were run to test whether the PTEN CNA findings depended on a single modeling choice. The analyses varied: (1) the structural-confidence filter (pLDDT 60, 70, and 80), (2) the number of closest cross-exon C-alpha distances used in the 3D metric (5 versus 10), (3) the minimum sequence separation for cross-exon residue pairs (0, 5, 10, and 20 residues), (4) the coding-coverage and exon-overlap rules used to define localizable events, (5) the effective-number-of-placement-windows diagnostic threshold, and (6) the mixture of uniform and span-matched CNA placement in the null model.

Across these reasonable parameter settings, the exact set and size of reported clusters changed—as expected when different CNA events or residue pairs became eligible—but the principal combined-CNA signal remained. The all-CNA analysis repeatedly produced highly significant joint results (minimum empirical `p_any_joint = 0.002` in the tested runs), with the strongest results consistently involving overlapping groups among exons 1–5. The pLDDT sweep retained 17 significant clusters at each threshold; sequence-separation and cross-exon-distance settings retained 12–15 and 12 significant clusters, respectively. These checks support that the observed PTEN exon-level pattern is not an artifact of one specific structural-distance or sequence-separation setting. The findings remain exploratory and conditional on a CNA already overlapping PTEN; they do not establish genome-wide gene recurrence.
