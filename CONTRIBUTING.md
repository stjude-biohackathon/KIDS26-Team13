# Contributing to KIDS26 Team 13

## Branch workflow

- `main` is the production branch. It should contain only reviewed, releasable work.
- `development` is the integration branch. Feature work is merged here first.
- Create short-lived branches from `development`, using names such as
  `feature/mutation-stability`, `fix/exon-mapping`, or `docs/annotation-workflow`.

Typical workflow:

```bash
git fetch origin
git switch development
git pull --ff-only origin development
git switch -c feature/short-description
```

Commit and publish the feature branch:

```bash
git add <specific-files>
git commit -m "Describe the focused change"
git push -u origin feature/short-description
```

Open a pull request from the feature branch into `development`. After integrated
work has been reviewed and validated together, open a separate release pull
request from `development` into `main`.

Do not push directly to `main` or `development`, rewrite their history, or mix
unrelated changes in one pull request.

## Before requesting review

Run the relevant checks from the repository root:

```bash
Rscript tests/test_GRIN3D_annotate_protein_clusters.R
Rscript tests/test_GRIN3D_mutation_stability.R
```

For statistical analyses, record the random seed, simulation count, input
version, coordinate source, and run settings. Do not interpret development runs
with small simulation counts as final inference.

Keep biological annotation separate from statistical significance. Generated
reports should retain source provenance and should not silently replace primary
analysis outputs.

## Pull-request destinations

- Normal change: `feature/*` → `development`
- Production release: `development` → `main`
- Urgent production fix: `fix/*` → `main`, followed by a synchronization pull
  request back into `development`

