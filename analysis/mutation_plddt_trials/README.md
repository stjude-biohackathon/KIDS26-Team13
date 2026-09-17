# Mutation hotspot sensitivity to trials and pLDDT

This analysis changes only two inputs to the GRIN3D mutation-hotspot model:

1. `n.sim`, the number of null-simulation trials; and
2. `min.confidence`, the AlphaFold pLDDT cutoff used to retain coordinates.

All other model settings are fixed in `01_run_sensitivity.R`: seed, mutation
events, minimum subject/event support, calibration fraction, significance
threshold, overlap rules, and the uniform positional null.

TP53 is included for trial-count sensitivity, but its current coordinate CSV
does not contain pLDDT. Its pLDDT threshold is therefore recorded as `NA` and
is not varied. SUZ12 and EZH2 contain residue-level pLDDT and are evaluated on
the full factorial grid.

## Run

From the repository root:

```bash
conda run -n kids26-team13 Rscript \
  analysis/mutation_plddt_trials/01_run_sensitivity.R quick

conda run -n kids26-team13 Rscript \
  analysis/mutation_plddt_trials/02_summarize_visualize.R

conda run -n kids26-team13 Rscript \
  analysis/mutation_plddt_trials/03_validate_analysis.R

conda run -n kids26-team13 Rscript \
  analysis/mutation_plddt_trials/04_analyze_trial_convergence.R

conda run -n kids26-team13 Rscript \
  analysis/mutation_plddt_trials/05_build_pdf_report.R
```

The `quick` profile uses 100 and 200 trials with pLDDT cutoffs 0, 70, and 90.
It is intended for workflow validation, not final inference. The `standard`
profile uses 200, 500, and 1,000 trials with cutoffs 0, 50, 70, and 90.
The `dense` profile uses every 100 trials from 100 through 1,000 with cutoffs
0, 70, and 90, matching the included extended results.
Custom comma-separated values can be supplied through environment variables:

```bash
GRIN3D_TRIALS=200,500 GRIN3D_PLDDT=0,70,90 \
  conda run -n kids26-team13 Rscript \
  analysis/mutation_plddt_trials/01_run_sensitivity.R custom
```

Completed run directories are reused, so an interrupted grid can be resumed.
Delete or move one run directory only when that configuration must be rerun.

## Reference and interpretation

Within each protein, the reference is the run with the largest requested trial
count and the least restrictive structural filter (`pLDDT >= 0`, or no filter
for TP53). Trial-count changes affect empirical p-value resolution but should
not change observed cluster membership at a fixed pLDDT cutoff. A pLDDT change
can alter both the eligible structure and the mutation events that remain fully
mapped; mapping retention is therefore reported beside cluster similarity.

The summary stage produces:

- `run_metrics.csv`: structural/mutation retention, cluster counts, runtime,
  significance, and recovery relative to the reference;
- `reference_cluster_matches.csv`: the best match for every reference cluster,
  with Jaccard, Dice, containment, precision/recall, boundary changes, p-value
  changes, and significance flips;
- `cluster_overlap_edges.csv`: every reference/comparison pair with nonzero
  residue overlap, preserving evidence of cluster splits and merges;
- `reference_cluster_stability.csv`: cluster-level reproducibility and an
  exploratory priority score;
- `validation_checks.csv`: automated checks that fixed settings stayed fixed
  and trial count alone did not change observed membership; and
- PNG/PDF figures showing global sensitivity, cluster similarity, empirical
  p-value trajectories, residue-level membership, and overlap transitions.

The `*_all_trials_cluster_overlap_atlas.png` figures place all trial counts in
one view per protein and print Jaccard directly in every cluster/configuration
cell. The `*_configuration_overlap_matrix.png` figures summarize similarity
between the complete cluster landscapes for every pair of runs.

`04_analyze_trial_convergence.R` then compares p-values and significance calls
to the 1,000-trial result at the same pLDDT. It reports p-value error, rank
correlation, significance agreement, and the earliest trial count after which
the exact significant-cluster set remains equal to the 1,000-trial set.

`05_build_pdf_report.R` creates a curated, self-contained PDF in `report/` from
the saved summaries. It does not rerun hotspot simulations.

The priority score is descriptive, not a statistical test. Biological meaning
and functional annotation remain separate from hotspot significance.
