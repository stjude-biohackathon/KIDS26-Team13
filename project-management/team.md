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
