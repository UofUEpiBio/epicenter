# Epicenter

Analysis workflow for the paper on bipartite contact networks between health
care providers (HCPs) and residents in long-term care facilities in the United
States.

The publication target is the supplemental materials document:

```bash
quarto render models/05-supplemental-materials.qmd
```

`models/05-supplemental-materials.qmd` is the only document needed to reproduce
the published supplement. It reads prepared network objects and saved model
fits rather than refitting anything: the upstream ERGM and bootstrap steps take
days of cluster time. The other Quarto documents are the intermediate steps that
produced those saved objects, and are kept for provenance.

## Analysis documents

| File | Purpose | Key output |
| --- | --- | --- |
| `data/model_data.Rmd` | Builds the 93 bipartite facility-unit networks from the raw interaction records. The raw CSV is not distributed with this repository. | `data/network93_f2.RData` |
| `models/01-descriptive.qmd` | Descriptive analysis of resident care-need attributes, degree, and network size. | figures only |
| `models/02-individual_ergms.qmd` | Exploratory single-network ERGMs, one fit per network, used to check which terms are estimable. | `models/02-individual_ergms_fits/` |
| `models/ergms_filtering_nets.qmd` | Early exploration of which networks to keep; documents dropping the 7 networks with a single HCP. | none (report) |
| `models/03-pooled-ergms.qmd` | Main model selection. Fits the pooled multi-network ERGMs term by term with `ergm.multi::N()`. | `models/03-pooled-ergms-final-results.rds`, `models/03-pooled-ergms-rds/` |
| `models/20260520-ergm-fit-diagnostics.qmd` | Investigation of convergence failures and between- vs within-state heterogeneity, motivating the final specification. | none (report) |
| `models/04-boot-ergms.qmd` | Bootstrap of the baseline model using MC-MLE. Its cache is read by `04b` for the SA-vs-MC-MLE comparison. | `models/04-boot-ergms-cache/` |
| `models/04b-boot-ergms.qmd` | Block bootstrap of the baseline model using stochastic approximation, compared against the MC-MLE run above. | `models/04b-boot-ergms-results.rds` |
| `models/04c-boot-ergms.qmd` | Block bootstrap of the model with state-level intercepts (`facility_state`). | `models/04c-boot-ergms-results.rds` |
| `models/04d-boot-ergms.qmd` | Block bootstrap of the final GWDSP × log(size) model. Supplies the standard errors reported in the paper. | `models/04d-boot-ergms-results.rds` |
| `models/05-supplemental-materials.qmd` | **Final supplement.** Observed statistics, GOF, MCMC diagnostics, bootstrap standard errors, and the alternative-specification comparison. | `models/05-supplemental-materials.docx` |
| `models/helpers.R` | Shared code: `set_ergm()` network preparation, `boot_ergm()` bootstrap driver, `facility_cluster_id()` blocking variable, GOF plotting, `texreg` wrappers, and result caching. | — |

The `04*` documents share one bootstrap driver, `boot_ergm()` in
`models/helpers.R`. Each replicate is cached as an individual `.rds` under
`models/04*-boot-ergms-cache/`, so an interrupted run resumes where it stopped.

`04b`, `04c`, and `04d` run a **block bootstrap**: they pass `cluster_id` to
`boot_ergm()`, which switches the driver from resampling individual networks to
resampling whole facilities, so every network belonging to a drawn facility
enters the replicate together and within-facility correlation is respected. The
blocking variable is built once by `facility_cluster_id()` in
`models/helpers.R`, which blocks on `state` + `Fac_Name` and maps the 93
networks onto 24 facilities. `04` is the older MC-MLE run and still resamples
individual networks; it is kept only as the comparison column in `04b`.

Because the resampling scheme is part of what each cached replicate encodes, a
cache written by an earlier naive run is not reusable: clear
`models/04*-boot-ergms-cache/` before re-running these documents, or the stale
replicates will be restored instead of recomputed.

## Inputs required by the supplement

These files must be present to render `models/05-supplemental-materials.qmd`:

- `data/network93_f2.RData` — the 93 prepared networks.
- `models/04b-boot-ergms-results.rds`, `models/04c-boot-ergms-results.rds`, and
  `models/04d-boot-ergms-results.rds` — block-bootstrap summaries; `04d` also
  supplies the final fitted model.
- `models/03-pooled-ergms-rds/res4_gwdsp_var_bnodematch.rds`,
  `res6_full_and_facility_mcmle.rds`, and `res6_full_and_logn_gwdsp.84mcmle.rds`
  — the three specifications in the comparison table.
- `models/05-supplemental-materials-gof.rds` — cached GOF object; recomputed with
  `gofN()` and saved if absent.

Model objects and prepared data are `.gitignore`d, so they travel outside git.

## Environment

The analysis runs in the container defined by `.devcontainer/Containerfile`
(R 4.5, Quarto 1.9.35, `ergm.multi`, `data.table`, `texreg`, `ggplot2` and
friends, and `tabulergm` from GitHub).

```bash
make container_build
make container_run
make render            # quarto render models/05-supplemental-materials.qmd
```

On CHPC, `render.slurm` runs the same target through Singularity. The `04*`
bootstrap documents use `parallel::makeForkCluster()`; lower the worker count
before running them anywhere other than a compute node.

## Repository layout

```
data/       prepared networks and the data-preparation source
models/     analysis documents, helpers.R, and saved model objects
.devcontainer/  container definition
Makefile    render, container, and Singularity targets
render.slurm    CHPC entry point
```

Superseded analyses and side projects (`abm/`, `chong/`, `fig/`,
`LTCF_final_report/`, `summary-stats.R`) are no longer tracked. They may still
exist in a local working copy, and remain in the git history.

## Note on terminology

Saved model objects and coefficient names use the legacy labels `n_patients` and
`n_hcp`. Throughout the manuscript these mean residents and health care
providers, respectively.
