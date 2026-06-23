# Epicenter

This repository contains the current analysis workflow for studying bipartite contact networks between health care workers (HCWs) and residents in long-term care facilities in the United States. The main statistical work fits pooled and individual ERGMs across 93 facility-level networks, with supporting descriptive analyses, diagnostics, bootstrap uncertainty estimation, and a separate agent-based modeling (ABM) subproject for outbreak simulations.

A large share of the repository consists of rendered reports, fitted model objects, and cached intermediate results. Many of the scripts intentionally cache expensive computations to disk so repeated runs are much faster.

## Project structure

### Top-level files

- `README.md`: This documentation.
- `README.html`: Rendered HTML copy of the README.
- `README_files/`: Supporting assets used by `README.html`.
- `Makefile`: Top-level helper targets for container builds and Singularity/CHPC execution. The container targets are still useful; the `render` target currently points to an older Quarto path and should be updated before relying on it.
- `ContainerFile`: Older standalone container recipe based on `rocker/r-ver:4.4.0`.
- `render.slurm`: SLURM batch script that loads Singularity and runs the top-level Make target on CHPC.
- `epicenter.sif`: Singularity image used for HPC execution.
- `epicenter.Rproj`: RStudio project file.
- `summary-stats.R`: Small utility script that loads a fitted model object and prints observed sufficient statistics used in the manuscript.
- `partial.csv`: Small intermediate data file kept at the repository root.
- `ergm-terms.svg`: Figure illustrating ERGM terms.
- `4d.zip`, `final-report.zip`: Bundled exports/archive artifacts.
- `Rplots.pdf`, `errs.pdf`, `errs-new.pdf`, `gof12_93-new.pdf`, `gof12_93-obs-vs-fitted.pdf`, `gofN_res2_b2factor_terms.pdf`, `mcmc_diagnostics.pdf`: Standalone diagnostic or figure outputs produced during model development.

### Environment and editor configuration

- `.devcontainer/Containerfile`: Current development container recipe. Installs R 4.5, Quarto, and the main analysis packages.
- `.devcontainer/devcontainer.json`: VS Code devcontainer definition.
- `.vscode/settings.json`: Workspace editor settings.
- `.gitignore`: Git ignore rules.

### Data inputs and derived analysis objects

The `data/` directory contains the analysis-ready objects used by the current Quarto documents.

- `data/model_data.Rmd`: Older data-preparation notebook showing how the raw interaction data were cleaned and converted into network objects.
- `data/network93_f2.RData`: Main analysis-ready list of 93 bipartite networks used throughout the current ERGM workflow.
- `data/big_net_mat_f2.RData`: Combined network object used for summaries and exploratory work.
- `data/dd_nodal_attr.RData`: Nodal attribute table used by older scripts and legacy analyses.
- `data/facilities.csv`: Facility-level metadata.
- `data/gofN_model12_93.rds`: Saved GOF object.
- `data/model_may2026.rds`, `data/model_may2026_log_n_simpler.rds`, `data/model_may2026_mean_centered.rds`: Saved model objects from the 2026 model-selection cycle.
- `data/res12_92.rda`, `data/res12_92.rds`, `data/res12_93.rda`, `data/res12_93.rds`: Saved fitted-model objects in both `.rda` and `.rds` formats.

Note: `data/model_data.Rmd` references a raw CSV (`data/Matched Interaction and Resident.csv`) that is not present in the repository, so the repo is currently set up for analysis from prepared network objects rather than full raw-data reconstruction.

### Main analysis workflow

The `models/` directory contains the active Quarto-based analysis pipeline plus many saved fits and rendered outputs.

#### Source files you are most likely to run

- `models/helpers.R`: Shared helper functions used across the modeling notebooks. This is the central utility file for preparing networks, tabulating fits, plotting GOF residuals, and supporting bootstrap/model-comparison work.
- `models/01-descriptive.qmd`: Descriptive analysis of network size, degree distributions, and resident feature distributions.
- `models/02-individual_ergms.qmd`: Fits simple ERGMs separately to each network and summarizes variability across network-specific fits.
- `models/03-pooled-ergms.qmd`: Main pooled multi-network ERGM analysis and model-selection document.
- `models/04-boot-ergms.qmd`: First bootstrap workflow for the pooled ERGM.
- `models/04b-boot-ergms.qmd`: Bootstrap variant with different stochastic approximation controls.
- `models/04c-boot-ergms.qmd`: Bootstrap variant that adds state-level effects and stores results in its own cache.
- `models/04d-boot-ergms.qmd`: Current bootstrap/specification variant used by the supplemental materials.
- `models/05-supplemental-materials.qmd`: Supplemental document that compiles diagnostics, observed statistics, bootstrap summaries, and additional checks around the final model.
- `models/20260520-ergm-fit-diagnostics.qmd`: Focused diagnostic report for model fit and residual checks.
- `models/ergms_filtering_nets.qmd`: Notes/checks related to how networks are filtered for specific terms.
- `models/03-pooled-ergms.rmarkdown`: Older R Markdown version of the pooled ERGM analysis retained for reference.

#### Rendered outputs and saved results

- `models/*.html`, `models/*.pdf`, `models/*.docx`, `models/*.md`: Rendered outputs of the Quarto/R Markdown source files.
- `models/03-pooled-ergms-final-results.rds`: Saved consolidated results object for the pooled ERGM analysis.
- `models/04-boot-ergms-results.rds`, `models/04b-boot-ergms-results.rds`, `models/04c-boot-ergms-results.rds`, `models/04d-boot-ergms-results.rds`: Saved bootstrap result objects.
- `models/04c-boot-ergms-initial.rds`: Initial object used in the `04c` bootstrap workflow.
- `models/05-supplemental-materials-gof.rds`, `models/05-supplemental-materials-sims.rds`: Saved objects used by the supplemental materials document.
- `models/02-individual_ergms_fits/`: Per-network cached ERGM fits written by `02-individual_ergms.qmd`.
- `models/03-pooled-ergms-rds/`: Cached intermediate pooled-model fits from the model-selection process.
- `models/03-pooled-ergms_cache/`: Cache directory for rendered pooled-model analyses.
- `models/04-boot-ergms-cache/`: Cached bootstrap replicates for `04-boot-ergms.qmd`.
- `models/03-pooled-ergms_files/`, `models/01-supplemental-materials_files/`: Figure/resource directories created during rendering.
- `models/example1.svg`, `models/size-distribution.png`, `models/errs.pdf`, `models/Rplots.pdf`, `models/05-supplemental-materials-network-feats.png`: Generated figures and exploratory graphics.

### ABM subproject

The `abm/` directory contains a separate outbreak simulation workflow built around Epiworld and network variants derived from the ERGM analysis.

- `abm/Makefile`: Build and execution helper for the ABM code.
- `abm/simoutbreak.cpp`: Main C++ simulation program.
- `abm/epiworld.hpp`: Vendored Epiworld header used by the simulation.
- `abm/simoutbreak-params.txt`: Parameter file consumed by the outbreak simulator.
- `abm/networks.R`: Generates network and attribute files for the simulator from ERGM outputs.
- `abm/actor_attributes.txt`: Exported actor/resident attribute file used by the simulator.
- `abm/abm.Rmd`: R Markdown report summarizing the simulation output.
- `abm/descriptive-stats.R`: Small descriptive script for degree summaries on a fitted model object.

### Legacy and archival analysis

These directories are useful for provenance and older manuscript iterations, but they are not the main current entry points.

- `chong/`: Earlier R and R Markdown analysis files from the 2023 model-development cycle.
  - `chong/2023-08-16-bipartite-ergms_multi.Rmd`, `chong/2023-08-17-bipartite-ergms_multi.Rmd`: Earlier manuscript/report drafts.
  - `chong/rcode.R`, `chong/rcode_chpc.R`: Older analysis scripts, including CHPC-oriented versions.
  - `chong/Final_models_assessment.R`: Earlier model assessment script.
  - `chong/school.rmd`: Additional legacy notebook.
- `LTCF_final_report/`: Archived fitted models, GOF objects, PDFs, and scripts associated with earlier manuscript versions.
  - Includes `res*.RData`, `gofN_*.RDS`, diagnostic PDFs, and `rcode_chpc.R`.
- `review/`: Manuscript-review materials.
  - `review/review_email.md`: Editorial decision and reviewer comments captured as Markdown.
  - `review/review.pdf`: PDF review document.
- `20250707-socnet-submission/SON-D-25-00410.pdf`: Archived manuscript submission PDF.
- `fig/`: PNG illustrations of ERGM terms used in documentation/presentation.

## How to run the project

### 1. Recommended local environment

The cleanest setup is the devcontainer or an equivalent R environment with Quarto.

Key software:

- R 4.5 or close equivalent
- Quarto
- Main R packages: `ergm.multi`, `ergm`, `data.table`, `netplot`, `ggplot2`, `gridExtra`, `texreg`, `ggrepel`, `patchwork`, `ggridges`, `ggExtra`
- For some legacy scripts: `Rglpk`, `tidyverse`, `intergraph`
- For the ABM subproject: a C++17 compiler

If you use the repository container setup:

```bash
make container_build
make container_run
```

These targets use `.devcontainer/Containerfile`.

### 2. Run the main Quarto analyses

From the repository root, the most reliable way is to render documents directly:

```bash
quarto render models/01-descriptive.qmd
quarto render models/02-individual_ergms.qmd
quarto render models/03-pooled-ergms.qmd
quarto render models/04d-boot-ergms.qmd
quarto render models/05-supplemental-materials.qmd
quarto render models/20260520-ergm-fit-diagnostics.qmd
```

Suggested order:

1. `models/01-descriptive.qmd`
2. `models/02-individual_ergms.qmd`
3. `models/03-pooled-ergms.qmd`
4. One of the bootstrap documents (`04*`)
5. `models/05-supplemental-materials.qmd`
6. `models/20260520-ergm-fit-diagnostics.qmd`

Important notes:

- These documents assume the prepared `.RData` and `.rds` files in `data/` and `models/` already exist.
- Several documents use `parallel::makeForkCluster(...)` with high core counts. You may want to reduce the requested number of workers on a laptop or shared machine.
- The top-level `make render` target currently references an older Quarto file path, so direct `quarto render ...` commands are the safer choice unless that Make target is updated.

### 3. Run the ABM workflow

From `abm/`:

```bash
make simoutbreak.o
make sims
```

Useful targets:

- `make networks`: Regenerates network/attribute exports from the ERGM-derived objects.
- `make simoutbreak.o`: Compiles the C++ simulator.
- `make sims`: Runs several simulation scenarios and then renders the ABM report.
- `make abm.pdf`: Renders `abm/abm.Rmd`.

The ABM workflow expects additional output directories and simulation result files to exist or be created during execution.

### 4. Run on CHPC / Singularity

For the HPC path, the repository includes both a Singularity image and a SLURM script:

```bash
sbatch render.slurm
```

This loads Singularity and runs the Make target `singularity_render_chpc`. As noted above, the container/HPC wiring is present, but the render target itself still points to an older document path and may need a small update before use in the current repository state.

## Caching and saved results

Many scripts in this repository automatically save intermediate or final results to disk to avoid recomputing expensive ERGM or bootstrap jobs.

Examples:

- `models/02-individual_ergms.qmd` writes per-network fits to `models/02-individual_ergms_fits/`.
- `models/03-pooled-ergms.qmd` restores and reuses saved model objects in `models/03-pooled-ergms-rds/`.
- `models/04*` bootstrap documents write replicate-level caches and final `.rds` summaries.
- Quarto itself also creates figure/resource folders and cache directories for some reports.

This behavior is intentional and improves computational efficiency substantially, especially for model fitting, diagnostics, and bootstrap replication.

## Practical guidance

- If you want to read the current main analysis, start with `models/03-pooled-ergms.qmd` and `models/05-supplemental-materials.qmd`.
- If you want to inspect shared logic, start with `models/helpers.R`.
- If you want to understand how the prepared network objects were originally constructed, read `data/model_data.Rmd`.
- If you want to reproduce the outbreak simulations, work inside `abm/`.
- If you are cleaning up or modernizing the repo further, the first candidate for maintenance is the top-level `Makefile` render target so it points at the current Quarto sources.
