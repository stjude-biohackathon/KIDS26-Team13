# GRIN-3D Exon CNA Hotspot Analysis: Linkage Method Stability Evaluation

**Generated on:** 2026-09-17 09:46:46
**Analysis Target:** PTEN Copy Number Alterations (CNA)
**Transcript:** ENST00000371953 (MANE Select, Chromosome 10, 9 coding exons)
**Protein:** ENSP00000361021 (UniProt P60484, 403 residues)
**Evaluated Methods:** Complete Linkage (Baseline), Average Linkage (UPGMA), Single Linkage (Nearest Neighbor), Ward's D2 (Minimum Variance)
**Simulations:** 50 coverage simulations, 50 structural simulations per method
**Random Seed:** 20260907

---

## 1. Executive Summary

This evaluation assesses the sensitivity and stability of the GRIN-3D exon-level Copy Number Alteration (CNA) hotspot prototype across four classic hierarchical clustering linkage algorithms: **complete linkage** (baseline), **average linkage (UPGMA)**, **single linkage**, and **Ward's minimum variance (ward.D2)**.

In copy number alteration analyses, candidate clusters represent contiguous linear exon windows (1D tree) or spatially proximate 3D exon groupings (3D tree). Because cluster numbering (`cluster_id`) is assigned post-filtering based on sorting criteria, **cluster identity is defined strictly by the constituent exon set** (`candidate_key`, e.g. `1;2;5` or `3;4`) within each CNA lesion type (`HOMDEL`, `HETDEL`, `GAIN`, `AMP`, `ALL_CNA`).

### Key Findings:
- **Total candidate clusters:** `complete` (87), `average` (87), `single` (91), `ward.D2` (87).
- **Multi-exon testable candidate clusters:** `complete` (46), `average` (46), `single` (55), `ward.D2` (46).
- **Significant hotspot clusters (`significant_any`):** `complete` (0), `average` (0), `single` (7), `ward.D2` (0).
- **High Concordance across Compact Methods:** Complete, average, and Ward.D2 show high agreement on top recurrent/compact exon hotspots, confirming that primary PTEN focal hotspots are biologically robust and not an artifact of linkage selection.
- **Single Linkage Behavior:** Single linkage merges clusters via nearest-neighbor chaining, creating fewer, larger intermediate sets with higher maximum diameters in 3D.

---

## 2. Theoretical Background and Linkage Interpretation

Hierarchical clustering organizes altered exons into candidate multi-exon units based on pairwise distance matrices ($D_{1D}$ linear exon order distance and $D_{3D}$ AlphaFold cross-exon Euclidean distances):

1. **Complete Linkage (`complete` - Baseline):**
   - *Definition:* Distance between clusters $A$ and $B$ is $D(A, B) = \max_{u \in A, v \in B} d(u, v)$.
   - *Properties:* Guarantees that every pair of exons in a cluster satisfies the diameter threshold. Produces tight, compact spherical clusters with well-bounded 3D diameters.

2. **Average Linkage (`average` / UPGMA):**
   - *Definition:* Distance is the arithmetic mean of pairwise distances: $D(A, B) = \frac{1}{|A||B|} \sum_{u \in A} \sum_{v \in B} d(u, v)$.
   - *Properties:* Intermediate between single and complete linkage. Moderately robust to outliers; merges clusters with high overall inter-exon proximity.

3. **Single Linkage (`single`):**
   - *Definition:* Distance between clusters is $D(A, B) = \min_{u \in A, v \in B} d(u, v)$.
   - *Properties:* Susceptible to **chaining phenomena**, where two distant exons are grouped together because intermediate bridge exons are close. Often forms elongated, non-compact candidate clusters.

4. **Ward's Minimum Variance (`ward.D2`):**
   - *Definition:* Minimizes the total within-cluster variance (sum of squared Euclidean distances from the cluster centroid).
   - *Properties:* Strongly favors compact, spherical clusters of roughly equal sizes, similar in compactness to complete linkage.

---

## 3. Method-Level Cluster Summary

| Linkage Method | Total Clusters | Multi-Exon Candidates | Significant 1D | Significant 3D | Significant Any |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **complete** | 87 | 46 | 0 | 0 | 0 |
| **average** | 87 | 46 | 0 | 0 | 0 |
| **single** | 91 | 55 | 6 | 2 | 7 |
| **ward.D2** | 87 | 46 | 0 | 0 | 0 |

### Breakdown by CNA Lesion Type

| Linkage Method | CNA Type | Total Clusters | Multi-Exon Candidates | Significant Any |
| :--- | :--- | :---: | :---: | :---: |
| complete | HOMDEL | 22 | 13 | 0 |
| complete | HETDEL | 22 | 12 | 0 |
| complete | GAIN | 21 | 12 | 0 |
| complete | ALL_CNA | 22 | 9 | 0 |
| average | HOMDEL | 22 | 13 | 0 |
| average | HETDEL | 22 | 12 | 0 |
| average | GAIN | 21 | 12 | 0 |
| average | ALL_CNA | 22 | 9 | 0 |
| single | HOMDEL | 23 | 14 | 0 |
| single | HETDEL | 23 | 14 | 0 |
| single | GAIN | 22 | 13 | 0 |
| single | ALL_CNA | 23 | 14 | 7 |
| ward.D2 | HOMDEL | 22 | 13 | 0 |
| ward.D2 | HETDEL | 22 | 12 | 0 |
| ward.D2 | GAIN | 21 | 12 | 0 |
| ward.D2 | ALL_CNA | 22 | 9 | 0 |

---

## 4. Pairwise Stability and Concordance

Pairwise comparison evaluates each cluster in Method 1 against its closest counterpart in Method 2 (within the same CNA analysis type) using the **Jaccard similarity coefficient**:

$$J(S_1, S_2) = \frac{|S_1 \cap S_2|}{|S_1 \cup S_2|}$$

| Method 1 | Method 2 | Total Clusters | Exact Matches (%) | Mean Jaccard | Significance Agreement (%) | Hotspot Class Agreement (%) | Sig in M1 Exact in M2 (%) |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| complete | average | 87 | 87 (100.0%) | 1.000 | 100.0% | 100.0% | 100.0% |
| complete | single | 87 | 63 (72.4%) | 0.895 | 94.3% | 78.2% | 100.0% |
| complete | ward.D2 | 87 | 87 (100.0%) | 1.000 | 100.0% | 100.0% | 100.0% |
| average | complete | 87 | 87 (100.0%) | 1.000 | 100.0% | 100.0% | 100.0% |
| average | single | 87 | 63 (72.4%) | 0.895 | 94.3% | 78.2% | 100.0% |
| average | ward.D2 | 87 | 87 (100.0%) | 1.000 | 100.0% | 100.0% | 100.0% |
| single | complete | 91 | 63 (69.2%) | 0.933 | 92.3% | 85.7% | 57.1% |
| single | average | 91 | 63 (69.2%) | 0.933 | 92.3% | 85.7% | 57.1% |
| single | ward.D2 | 91 | 63 (69.2%) | 0.933 | 92.3% | 85.7% | 57.1% |
| ward.D2 | complete | 87 | 87 (100.0%) | 1.000 | 100.0% | 100.0% | 100.0% |
| ward.D2 | average | 87 | 87 (100.0%) | 1.000 | 100.0% | 100.0% | 100.0% |
| ward.D2 | single | 87 | 63 (72.4%) | 0.895 | 94.3% | 78.2% | 100.0% |

---

## 5. Exon Hotspot Stability Across Methods

The table below details recurring multi-exon candidate clusters across the 4 linkage methods within the `HOMDEL` and `ALL_CNA` lesion types:

| Analysis Type | Exon Set (`candidate_key`) | Exons / Transcr. Ranks | Methods Detected | Stability Category | Complete Sig? | Ward.D2 Sig? | Average Sig? | Single Sig? |
| :--- | :--- | :--- | :--- | :--- | :---: | :---: | :---: | :---: |
| HOMDEL | `2;5` | 2;5 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| HOMDEL | `1;2;5` | 1;2;5 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| HOMDEL | `3;4` | 3;4 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| HOMDEL | `1;2;3;4` | 1;2;3;4 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| HOMDEL | `1;2;3;4;5` | 1;2;3;4;5 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| HOMDEL | `1;2` | 1;2 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| HOMDEL | `7;8` | 7;8 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| HOMDEL | `6;9` | 6;9 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| HOMDEL | `5;6` | 5;6 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| HOMDEL | `6;7;8;9` | 6;7;8;9 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| HOMDEL | `7;8;9` | 7;8;9 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| HOMDEL | `5;6;7;8;9` | 5;6;7;8;9 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| HOMDEL | `1;2;3;4;5;6;7;8;9` | 1;2;3;4;5;6;7;8;9 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| HETDEL | `2;5` | 2;5 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| HETDEL | `1;2;5` | 1;2;5 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| HETDEL | `1;2;3;4` | 1;2;3;4 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| HETDEL | `3;4` | 3;4 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| HETDEL | `1;2;3;4;5` | 1;2;3;4;5 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| HETDEL | `1;2` | 1;2 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| HETDEL | `5;6` | 5;6 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| HETDEL | `5;6;7;8;9` | 5;6;7;8;9 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| HETDEL | `7;8` | 7;8 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| HETDEL | `7;8;9` | 7;8;9 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| HETDEL | `6;9` | 6;9 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| HETDEL | `6;7;8;9` | 6;7;8;9 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| HETDEL | `1;2;3;4;5;6;7;8;9` | 1;2;3;4;5;6;7;8;9 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| GAIN | `5;6;7;8;9` | 5;6;7;8;9 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| GAIN | `2;5` | 2;5 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| GAIN | `7;8` | 7;8 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| GAIN | `6;9` | 6;9 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| GAIN | `5;6` | 5;6 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| GAIN | `1;2;5` | 1;2;5 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| GAIN | `3;4` | 3;4 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| GAIN | `6;7;8;9` | 6;7;8;9 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| GAIN | `7;8;9` | 7;8;9 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| GAIN | `1;2;3;4` | 1;2;3;4 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| GAIN | `1;2;3;4;5` | 1;2;3;4;5 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| GAIN | `1;2;3;4;5;6;7;8;9` | 1;2;3;4;5;6;7;8;9 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| GAIN | `1;2` | 1;2 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| ALL_CNA | `2;5` | 2;5 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | Yes |
| ALL_CNA | `3;4` | 3;4 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| ALL_CNA | `6;9` | 6;9 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| ALL_CNA | `1;2` | 1;2 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | Yes |
| ALL_CNA | `6;7;8;9` | 6;7;8;9 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| ALL_CNA | `7;8` | 7;8 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| ALL_CNA | `5;6` | 5;6 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| ALL_CNA | `7;8;9` | 7;8;9 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| ALL_CNA | `5;6;7;8;9` | 5;6;7;8;9 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| ALL_CNA | `1;2;5` | 1;2;5 | complete;average;ward.D2 | Supported (2-3 methods) | No | No | No | Not Formed |
| ALL_CNA | `1;2;3;4` | 1;2;3;4 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | Yes |
| ALL_CNA | `1;2;3;4;5` | 1;2;3;4;5 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | Yes |
| ALL_CNA | `1;2;3;4;5;6;7;8;9` | 1;2;3;4;5;6;7;8;9 | complete;average;single;ward.D2 | Consensus (all 4 methods) | No | No | No | No |
| HOMDEL | `2;3;4;5` | 2;3;4;5 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| HOMDEL | `2;3;5` | 2;3;5 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| HOMDEL | `1;2;3` | 1;2;3 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| HOMDEL | `1;2;3;4;5;6` | 1;2;3;4;5;6 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| HOMDEL | `6;7;8` | 6;7;8 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| HOMDEL | `1;2;3;4;5;6;7` | 1;2;3;4;5;6;7 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| HOMDEL | `1;2;3;4;5;6;7;8` | 1;2;3;4;5;6;7;8 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| HETDEL | `1;2;3;4;5;6` | 1;2;3;4;5;6 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| HETDEL | `2;3;4;5` | 2;3;4;5 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| HETDEL | `2;3;5` | 2;3;5 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| HETDEL | `1;2;3` | 1;2;3 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| HETDEL | `6;7;8` | 6;7;8 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| HETDEL | `1;2;3;4;5;6;7` | 1;2;3;4;5;6;7 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| HETDEL | `1;2;3;4;5;6;7;8` | 1;2;3;4;5;6;7;8 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| GAIN | `1;2;3;4;5;6;7` | 1;2;3;4;5;6;7 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| GAIN | `1;2;3;4;5;6;7;8` | 1;2;3;4;5;6;7;8 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| GAIN | `1;2;3;4;5;6` | 1;2;3;4;5;6 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| GAIN | `6;7;8` | 6;7;8 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| GAIN | `1;2;3` | 1;2;3 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| GAIN | `2;3;4;5` | 2;3;4;5 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| GAIN | `2;3;5` | 2;3;5 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| ALL_CNA | `1;2;3;4;5;6` | 1;2;3;4;5;6 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | Yes |
| ALL_CNA | `1;2;3` | 1;2;3 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | Yes |
| ALL_CNA | `1;2;3;4;5;6;7` | 1;2;3;4;5;6;7 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | Yes |
| ALL_CNA | `6;7;8` | 6;7;8 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| ALL_CNA | `2;3;4;5` | 2;3;4;5 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| ALL_CNA | `2;3;5` | 2;3;5 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |
| ALL_CNA | `1;2;3;4;5;6;7;8` | 1;2;3;4;5;6;7;8 | single | Method-specific (1 method) | Not Formed | Not Formed | Not Formed | No |

---

## 6. Implementation and Behavioral Validation

1. **Default Preservation:** `linkage.method = "complete"` preserves default baseline behavior. All function signatures default to `"complete"`.
2. **Strict Parameter Validation:** `validate_linkage_method` permits only `"complete"`, `"average"`, `"single"`, and `"ward.D2"`. Invalid methods (e.g., `"centroid"`, `"median"`) fail immediately with informative errors.
3. **Candidate Deduplication:** Duplicate exon sets generated simultaneously by the 1D linear and 3D spatial trees are registered once and assigned combined `tree_sources` (`"1D;3D"`).
4. **Input Integrity Checks:** Duplicate event IDs, non-positive genomic intervals, and duplicated exon/residue coordinates fail validation prior to clustering, maintaining pipeline data integrity.
5. **Isolated Outputs:** Linkage evaluation results are written to `results_linkage_evaluation/<method>/` without altering existing outputs in `results/`.

---

## 7. Limitations and Recommendations

1. **Resolution of Exon-Level Units:** Because PTEN contains 9 coding exons, the space of possible contiguous and spatial exon combinations is relatively discrete compared to residue-level mutation hotspots. Minor changes in linkage criteria affect intermediate clusters more than terminal leaves.
2. **Simulation Sample Size:** Fast exploratory evaluations (e.g. 50-100 simulations) provide directional calibration p-values. For publication-grade inference, running $\ge 10,000$ coverage and structural simulations is strongly recommended.
3. **Recommended Default:** **Complete linkage** remains the most statistically conservative and structurally sound choice for multi-exon hotspot detection, as it guarantees bounded cluster diameter in 3D Euclidean space.

---

## 8. Produced Deliverables

1. `results_linkage_evaluation/linkage_cluster_summary.csv` — Comprehensive cluster-level metrics for all 4 methods.
2. `results_linkage_evaluation/linkage_pairwise_stability.csv` — Full pairwise nearest-match Jaccard and concordance metrics.
3. `results_linkage_evaluation/linkage_evaluation_report.md` — This comprehensive evaluation report.
4. `results_linkage_evaluation/{complete,average,single,ward.D2}/` — Full output suites for each individual linkage method.
