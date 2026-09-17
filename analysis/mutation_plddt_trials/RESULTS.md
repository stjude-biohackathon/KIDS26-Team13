# Dense trial and pLDDT sensitivity results

A self-contained, executable R notebook is available as
`mutation_plddt_trials_report.Rmd`, with a portable rendered handoff at
`mutation_plddt_trials_report.nb.html`. It organizes each experiment into
separate narrative, calculation, figure, and conclusion sections without the
fixed-page limitations of PDF. The older curated PDF remains available at
`report/mutation_hotspot_stability_report.pdf`. Neither report reruns the
hotspot simulations.

The completed dense grid contains 150 runs. It evaluates every 100 null trials
from 100 through 1,000 for TP53, PTEN, SUZ12, EZH2, and LEF1, crossed with pLDDT
cutoffs of 0, 70, and 90. TP53 now uses the canonical AlphaFold DB v6 model for
UniProt P04637. Its 393 XYZ coordinates are exactly identical to the legacy
TP53 coordinate table; the new table adds the model's residue-level pLDDT
values, so the TP53 rerun isolates confidence filtering rather than a structure
change.

The reference for each protein is its 1,000-trial, least-restrictive run.

## TP53 pLDDT provenance

- Source model: AlphaFold DB `AF-P04637-F1-model_v6.pdb`, canonical UniProt
  P04637 entry, downloaded 2026-09-17 from
  <https://alphafold.ebi.ac.uk/entry/P04637>.
- The accompanying API response is preserved as
  `examples/TP53_mutations/input_files/alphafold_P04637_metadata.json`.
- `grin3d_alphafold_pdb_to_coordinates()` extracted each alpha-carbon XYZ
  coordinate and read pLDDT from the PDB B-factor field, producing
  `tp53_alphafold_v6_coordinates.csv`.
- The result contains 393 unique, contiguous residues (1--393). Its XYZ values
  have maximum absolute difference 0 from the legacy `tp53_coordinates.csv`.
  Only the pLDDT column is new.
- AlphaFold's pLDDT interpretation and PDB encoding are documented in the
  [AlphaFold DB FAQ](https://alphafold.ebi.ac.uk/faq).

## Main observations

- Trial count never changed the observed residue membership of a cluster at a
  fixed structural threshold. Every within-threshold comparison from 100 to
  1,000 trials had Jaccard 1.0. This is expected because trial count changes
  null calibration, not the deterministic observed clustering tree.
- Trial count did change empirical p-values and significance calls. TP53 had
  no significant clusters through 700 trials and 11 at 800--1,000 trials.
  This step behavior reflects increasing Monte Carlo resolution and clusters
  crossing the corrected significance threshold; it does not reflect cluster
  membership appearing at 800 trials.
- For TP53, pLDDT >= 70 retained 235/393 coordinates and 26/30 events. Mean
  best-match Jaccard against the unfiltered reference was 0.880, and 87.8% of
  reference clusters were recovered at Jaccard >= 0.8. Six clusters were
  significant at 1,000 trials, with the exact final significance set stable
  from 700 trials onward.
- At pLDDT >= 90, TP53 retained 207/393 coordinates and 25/30 events. Mean
  Jaccard remained 0.873, with 85.4% high-Jaccard recovery. Six clusters were
  significant at 1,000 trials; the 800-trial run had one additional call, and
  the exact final set was stable from 900 trials onward.
- PTEN retained 313/403 coordinates and 309/313 mutation events at pLDDT >= 70.
  Mean best-match Jaccard was 0.973, with 96.9% of unfiltered reference
  clusters recovered at Jaccard >= 0.8.
- At pLDDT >= 90, PTEN retained 277/403 coordinates and 303/313 events. Mean
  Jaccard remained 0.922, with 90.1% high-Jaccard recovery. No corrected
  significant PTEN clusters occurred at any evaluated threshold or trial count.
- At pLDDT >= 70, SUZ12 retained 453/739 coordinates and 49/72 events. Mean
  best-match Jaccard against the unfiltered reference was 0.736, with 61.5% of
  reference clusters recovered at Jaccard >= 0.8.
- At pLDDT >= 90, SUZ12 retained 279/739 coordinates and 32/72 events. Mean
  Jaccard fell to 0.455 and high-Jaccard recovery to 32.3%.
- At pLDDT >= 70, EZH2 retained 513/746 coordinates and 77/95 events. Mean
  Jaccard was 0.775, with 64.5% of reference clusters recovered at Jaccard >=
  0.8.
- At pLDDT >= 90, EZH2 retained 326/746 coordinates and 42/95 events. Mean
  Jaccard fell to 0.380 and high-Jaccard recovery to 20.0%.
- At pLDDT >= 70, LEF1 retained 84/399 coordinates and 27/50 events. Mean
  best-match Jaccard was 0.511, with 42.0% of unfiltered reference clusters
  recovered at Jaccard >= 0.8.
- At pLDDT >= 90, LEF1 retained 63/399 coordinates and 26/50 events. Mean
  Jaccard was 0.484, with 42.0% high-Jaccard recovery. Its significance sets
  at all three thresholds matched the 1,000-trial result from 900 trials onward.
- At 1,000 trials, corrected significance counts at pLDDT 0, 70, and 90 were
  11, 6, and 6 for TP53; 0, 0, and 0 for PTEN; 0, 0, and 2 for SUZ12; and 0,
  25, and 9 for EZH2; and 9, 6, and 1 for LEF1. The filtered results operate on
  different retained event cohorts, so these counts must not be interpreted as
  stronger evidence caused by higher confidence alone.

With 1,000 total trials and a 50/50 calibration split, the minimum evaluation
empirical p-value is approximately 0.002 before protein-wide and joint search
correction. This grid is appropriate for sensitivity assessment but remains
below the prototype's recommended final-inference simulation count.

## One-view overlap figures

The most direct figures for comparing all trials are:

- `results/figures/TP53_all_trials_cluster_overlap_atlas.png`
- `results/figures/PTEN_all_trials_cluster_overlap_atlas.png`
- `results/figures/SUZ12_all_trials_cluster_overlap_atlas.png`
- `results/figures/EZH2_all_trials_cluster_overlap_atlas.png`
- `results/figures/LEF1_all_trials_cluster_overlap_atlas.png`

Each row is a reference cluster, each column is a trial count, and each cell is
labelled with its Jaccard score. All five proteins have separate pLDDT panels in
the same image. Constant values across a panel show that trial count did not
alter residue membership; differences between panels show pLDDT-driven
membership changes.

The `*_configuration_overlap_matrix.png` figures provide a second overview.
They compare the complete candidate-cluster landscapes between every pair of
configurations using symmetric mean best-match Jaccard. The corresponding
numeric matrices are saved in `results/summary/`.

For more detailed changes:

- `*_cluster_residue_maps.pdf` shows residue additions and losses;
- `*_cluster_overlap_transitions.pdf` preserves one-to-many edges so splits and
  merges remain visible; and
- `*_pvalue_trajectories.pdf` separates p-value behavior from membership.

## Tables

- `results/summary/run_metrics.csv`: one row per configuration;
- `results/summary/reference_cluster_matches.csv`: best reference match and
  Jaccard, Dice, containment, precision/recall, boundary, p-value, and
  significance-change metrics;
- `results/summary/cluster_overlap_edges.csv`: every nonzero residue overlap,
  including split/merge relationships;
- `results/summary/reference_cluster_stability.csv`: descriptive cluster-level
  reproducibility ranking;
- `results/summary/*_configuration_overlap_matrix.csv`: all-configuration
  landscape similarity; and
- `results/summary/validation_checks.csv`: reproducibility and configuration
checks. All checks passed for all 150 runs.

The cluster priority score is only a sorting aid. It is not a p-value or
biological evidence, and functional annotation remains separate.

## Trial convergence relative to 1,000 trials

The follow-up convergence analysis compares each run with the 1,000-trial run
at the same pLDDT, so structural filtering is held constant within every
comparison. It records mean, median, root-mean-square, and maximum absolute
p-value error; p-value rank correlation; significance agreement; and exact
significance-set matches.

The earliest trial counts whose significant-cluster sets equal the 1,000-trial
set and remain equal at every subsequent evaluated count were:

| Protein | pLDDT | Earliest persistent match | Significant clusters at 1,000 |
|---|---:|---:|---:|
| TP53 | 0 | 800 | 11 |
| TP53 | 70 | 700 | 6 |
| TP53 | 90 | 900 | 6 |
| PTEN | 0 | 100 | 0 |
| PTEN | 70 | 100 | 0 |
| PTEN | 90 | 100 | 0 |
| SUZ12 | 0 | 100 | 0 |
| SUZ12 | 70 | 100 | 0 |
| SUZ12 | 90 | 1,000 | 2 |
| EZH2 | 0 | 100 | 0 |
| EZH2 | 70 | 1,000 | 25 |
| EZH2 | 90 | 600 | 9 |
| LEF1 | 0 | 900 | 9 |
| LEF1 | 70 | 900 | 6 |
| LEF1 | 90 | 900 | 1 |

The 100-trial matches with zero significant reference clusters are trivial:
every evaluated run had an empty significance set. A first match at exactly
1,000 trials is not evidence of convergence because there is no later trial
count to confirm persistence. SUZ12 pLDDT 90 and EZH2 pLDDT 70 therefore need
more trials before their significance sets should be treated as stable.

The convergence figures are `*_trial_convergence_overview.png` and
`*_significance_by_trial.png`. Numeric results are in
`trial_convergence_metrics.csv`, `cluster_pvalue_convergence.csv`, and
`earliest_significance_convergence.csv`.

All trial series use the same random seed, which makes increasing runs nested
and useful for resolution diagnostics. Separate-seed runs would be required to
measure Monte Carlo variability independently.
