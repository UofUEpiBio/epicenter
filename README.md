# Epicenter

This repository contains the analysis workflow for the supplemental materials to
the paper on bipartite contact networks between health care workers (HCWs) and
residents in long-term care facilities in the United States.

The current publication target is:

```bash
quarto render models/05-supplemental-materials.qmd
```

The supplement is intentionally built from prepared network objects and saved
model fits. Recreating every upstream ERGM and bootstrap object is possible, but
expensive; for routine manuscript work, the saved `.RData` and `.rds` inputs are
the practical starting point.

## Publication Inputs

`models/05-supplemental-materials.qmd` directly uses these files:

- `models/05-supplemental-materials.qmd`: final supplemental-materials source.
- `models/helpers.R`: shared ERGM setup, plotting, tabulation, and cache helpers.
- `data/network93_f2.RData`: prepared list of 93 bipartite facility networks.
- `models/04d-boot-ergms-results.rds`: final bootstrap result object; provides the final fitted model and bootstrap uncertainty estimates.
- `models/05-supplemental-materials-gof.rds`: cached GOF object. If absent, the supplement recomputes it with `gofN(naive_fit)` and saves it.
- `models/03-pooled-ergms-rds/res4_gwdsp_var_bnodematch.rds`: baseline comparison model.
- `models/03-pooled-ergms-rds/res6_full_and_facility_mcmle.rds`: state-intercept comparison model.
- `models/03-pooled-ergms-rds/res6_full_and_logn_gwdsp.84mcmle.rds`: final log-size GWDSP comparison model.
- `models/04b-boot-ergms-results.rds`, `models/04c-boot-ergms-results.rds`, and `models/04d-boot-ergms-results.rds`: bootstrap summaries used in the alternative-specification table.

The render creates manuscript outputs and figures such as
`models/05-supplemental-materials.docx`,
`models/05-supplemental-materials.html`,
`models/05-supplemental-materials.md`, and
`models/05-supplemental-materials_files/`. These are generated artifacts, not
source files.

## Requirements

The preferred environment is the devcontainer defined by:

- `.devcontainer/Containerfile`
- `.devcontainer/devcontainer.json`

The container installs R 4.5, Quarto 1.9.35, `libglpk-dev`, and the main R
packages used by the analysis:

- `ergm.multi`, `ergm`, `ergmito`, `network`, and `sna`
- `data.table`
- `ggplot2`, `ggrepel`, `ggridges`, `ggExtra`, `gridExtra`, and `patchwork`
- `texreg`
- `netplot`
- `quarto` and `knitr`
- `tabulergm`, installed from GitHub

Build and enter the container from the repository root with:

```bash
make container_build
make container_run
```

Then render the supplement:

```bash
quarto render models/05-supplemental-materials.qmd
```

or equivalently:

```bash
make render
```

## Execution Order

For normal manuscript edits, run only the final document after confirming the
saved inputs listed above are present:

```bash
quarto render models/05-supplemental-materials.qmd
```

To regenerate upstream objects from the prepared networks, use this order:

1. `quarto render models/03-pooled-ergms.qmd`
2. `quarto render models/04-boot-ergms.qmd`
3. `quarto render models/04b-boot-ergms.qmd`
4. `quarto render models/04c-boot-ergms.qmd`
5. `quarto render models/04d-boot-ergms.qmd`
6. `quarto render models/05-supplemental-materials.qmd`

Optional context reports can be rendered independently:

- `models/01-descriptive.qmd`
- `models/02-individual_ergms.qmd`
- `models/20260520-ergm-fit-diagnostics.qmd`
- `models/ergms_filtering_nets.qmd`

These optional reports help with provenance, diagnostics, and model selection,
but they are not direct dependencies of the final supplemental-materials file.

## Project Layout

- `models/`: current Quarto analysis workflow and saved model objects.
- `models/helpers.R`: shared code used by the active workflow.
- `data/`: prepared analysis inputs. The raw interaction CSV referenced by
  `data/model_data.Rmd` is not present, so the repository is not currently a
  raw-data-to-paper reconstruction.
- `.devcontainer/`: active reproducible development environment.
- `Makefile`: helper targets for rendering, container builds, and Singularity.
- `render.slurm`: CHPC/SLURM entry point that uses Singularity and `make render`.
- `abm/`: separate agent-based outbreak simulation subproject. It is useful for
  related work, but it does not feed `models/05-supplemental-materials.qmd`.
- `fig/`: ERGM term illustrations used for documentation or presentation.
- `review/` and `20250707-socnet-submission/`: manuscript review/submission
  records, useful for provenance but not part of the supplement build.
- `chong/` and `LTCF_final_report/`: earlier analysis and report archive from
  prior model-development cycles.

## Generated And Cached Files

The repository contains many rendered reports and saved computations. They are
useful for reproducibility and audit trails, but they should be treated as
generated artifacts unless the paper archive intentionally includes them.

Common generated outputs:

- `models/*.html`, `models/*.pdf`, `models/*.docx`, and `models/*.md`
- `models/*.rmarkdown`
- `README.html` and `README_files/`
- `models/*_files/`
- `models/*_cache/`
- `models/02-individual_ergms_fits/`
- `models/03-pooled-ergms-rds/`
- `models/04*-boot-ergms-cache/`
- `models/05-supplemental-materials-*.png`
- `models/05-supplemental-materials-*.pdf`

Common saved model/data objects:

- `*.RData`, `*.rda`, `*.rds`, and `*.RDS`
- `*.csv`, `*.tsv`, `*.dat`, and `*.data`
- `*.zip` and `*.sif`

These patterns are ignored by `.gitignore` for new files. Existing local copies
may still be needed to render the supplement, especially the prepared network
and model objects listed in "Publication Inputs".

## Cleanup And Archive Notes

Files and directories that do not contribute to
`models/05-supplemental-materials.qmd` and should be archived or removed from a
publication source release:

- `chong/`: old 2023 analysis scripts and notebooks.
- `LTCF_final_report/`: old fitted models, GOF objects, scripts, and PDFs.
- `20250707-socnet-submission/` and `review/`: submission/review records.
- Top-level rendered or packaged artifacts: `Epicenter.pdf`, `Rplots.pdf`,
  `errs*.pdf`, `gof*.pdf`, `mcmc_diagnostics.pdf`, `4d.zip`,
  `final-report.zip`, and `response_letter_edited.docx`.
- `README.html` and `README_files/`: generated from the README.
- `models/05-supplemental-materials-sims.rds`,
  `models/05-supplemental-materials-example1.svg`, and
  `models/example1.svg`: remnants of an older simulated-network visualization
  section that is no longer present in the current supplement.
- `data/big_net_mat_f2.RData`, `data/dd_nodal_attr.RData`,
  `data/gofN_model12_93.rds`, `data/model_may2026*.rds`, and
  `data/res12_*`: useful for older diagnostics and model-selection history, but
  not read by the final supplement.
- `models/04-boot-ergms-results.rds` and `models/04-boot-ergms-cache/`: earlier
  bootstrap outputs. They are not direct supplement inputs, but the
  `04b` report reads the `04-boot` cache if you regenerate reports from scratch.
- `models/*.rmarkdown`: older/intermediate R Markdown copies of Quarto sources.
- `models/.gitignore`: redundant local ignore file; the root `.gitignore`
  covers Quarto outputs and generated artifacts.
- `ContainerFile`: older top-level container recipe. The active recipe is
  `.devcontainer/Containerfile`.
- `epicenter.sif` and `.singularity/`: local Singularity image/cache artifacts.
- `summary-stats.R` and `abm/descriptive-stats.R`: small legacy utility scripts
  that load archived `LTCF_final_report` model objects.
- `ergm-terms.svg` and `fig/`: documentation/presentation assets; keep only if
  the manuscript or a related explainer needs them.
- `partial.csv` at the repository root: appears to be an intermediate data file;
  it is not read by the final supplement.
- `data/facilities.csv`: facility metadata, but not read directly by the final
  supplement.
- `data/data`, `models/data`, and `models/models`: symlinks back into the
  repository. They may be convenient locally, but they are easy to confuse with
  real source directories in a publication archive.
- `.vscode/settings.json`: duplicate local editor settings; the devcontainer
  also mounts `.devcontainer/.vscode/settings.json`.

The `abm/` directory is best handled as a separate subproject. In particular,
`abm/networks.R` references `../models/2022-04-25-bipartite-ergms.rds`, which is
not present in this repository, so that workflow is not currently reproducible
from the publication inputs alone.

## Notes

- Several Quarto documents use parallel workers and can be heavy on a laptop or
  shared machine. Reduce the worker count before rerunning long model fits.
- The saved model objects use legacy names such as `n_hcp` and `n_patients`.
  The manuscript language now treats those as HCWs and residents, respectively.
- The prepared data/model binaries are ignored for new commits by default. If a
  formal reproducibility archive needs to include them, add them explicitly with
  a documented exception rather than loosening the general ignore rules.
