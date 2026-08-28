# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

An R codebase for simulating cane toad (*Rhinella marina*) spread across waterpoints in north-western Australia. It supports Dunlop *et al.* (submitted 2025), which extends the point-process spread model of Tingley *et al.* (2013, *J. Appl. Ecol* 50:129-137) to:

1. Estimate arrival time of toads in the Pilbara / at the Toad Containment Zone (TCZ) under a do-nothing scenario.
2. Optimise investment in the TCZ (fencing/converting artificial waterpoints, accounting for natural waterpoint drying probability) to test containment scenarios.

This is a research script repository, not a package — there is no `DESCRIPTION`/`renv.lock`, no test suite, and no build/lint tooling. Everything runs interactively from RStudio (`TCZ-sims.Rproj`) or via `Rscript`, and scripts are meant to be `source()`d in a specific order (see below), populating the global environment as they go.

## Running things

There is no build, lint, or test command. To run a simulation:

```r
# from repo root, in an R session
source("src/tcz-sims/tcz-sims.R")      # full TCZ implementation scenario
source("src/tcz-sims/forecast.R")      # do-nothing forecast from current invasion front
```

Both of these `source("src/tcz-sims/load-data.R")` internally, which pulls in raw spatial data. **`load-data.R` requires the `DATA_PATH` environment variable** to be set to a local directory containing the Fulcrum infrastructure audit export, GIS layers, and invasion-front parameters (`Toads/TCZ/infrastructure/...`, `GIS - General/GIS_layers_read_only/...`, `invasion-front-parameters.Rdata`) — none of this lives in the repo. `dea-drying-probs.R` and `rainfall-drying-model.R` additionally need network access (DEA WFS, SILO — the latter via `SILO_API_KEY`) and are meant to be run once, with their output objects (`living.waters.dea`, `annual_metrics`, `fit_drying`, `rain_scale`) persisted into the session/posteriors rather than re-run every time.

For the HPC (cluster) workflow, `src/cluster/a-setup-environment.sh` loads the required environment modules (`r/4.3.0`, `gdal`, `geos`, `proj`, `udunits`) before invoking scripts with `commandArgs`.

The Shiny app (`shiny/app.R`) is a separate deployable unit: run `shiny/prep_app_data.R` after a `forecast.R` run to regenerate `shiny/data/app_data.rds`, then launch the app from the `shiny/` directory (it has its own `.Rproj`). The app is deliberately thin — it `readRDS()`s `app_data.rds` once and only filters/redraws Leaflet markers by year; all spatial work (`sf`, Albers reprojection, per-year colonisation probabilities) happens up-front in `prep_app_data.R`. That means it can be shipped as a **static shinylive site** (`Rscript shiny/build_site.R` → `shiny/site/`, gitignored) that runs entirely in the browser via webR and needs no Shiny Server — the preferred way to host it for public / high-traffic use. `prep_app_data.R` runs both inside a `forecast.R` session (boundary sf objects already in scope) and standalone (it then reads `tcz.boundary` from `DATA_PATH` and `pastoral.boundaries` from `out/pastoral-boundaries.shp`). `build_site.R`'s first run downloads ~400 MB of shinylive web assets into a user cache.

Manuscript output (`ms/tcz-optimisation-methods.qmd`) renders via Quarto (`_quarto.yml` sets `output-dir: ms/output`); it references figures written to `out/` by the simulation scripts, so simulations must be run before rendering.

## Architecture

### Core simulation engine: `src/pprocess_functions.R`

Every scenario script sources this file. It is the single workhorse module with no other internal dependencies within `src/`:

- **`setup()`** — builds the `spread.table` (matrix of waterpoint ID/X/Y/target/kernel-params/presence/age/natural-flag) and the sparse pairwise dispersal matrix (`pairs`), from a point dataset, kernel parameters (`dat/Kernel-fits_truncated.RData`), and posterior draws (`dat/Posteriors*.RData`). Both outputs are placed directly into the global environment via `list2env`. Handles spatial-duplicate removal (`remove_spatial_duplicates`) and optional TCZ-area filtering. Also unconditionally adds `occ.years`/`first.col.gen` bookkeeping columns, and (only when the corresponding optional arguments are supplied) `control.type`/`r.eff`/`fence.units`/`p.breach`/`surv.prob`/`control.fail.prob` columns that drive the optional control step in `spread.pilb()` — see below. All of these default to no-op values (`CONTROL_NONE`, `NA` `r.eff`, `surv.prob=1`, etc.) when the relevant `setup()` argument is omitted, so every scenario script that doesn't pass them is unaffected.
- **`spread.pilb()`** — advances the population one generation at a time: draws Poisson propagule counts from occupied sites, spreads them through the dispersal kernel matrix, draws realised colonisation via a multinomial, and stops early once any/all target site(s) are colonised. Optional args (all default off): `control.step=TRUE` turns on a within-timestep "control step" after colonisation each generation — natural points survive via `surv.prob` (a static 1 = never dries, until a real drying model is wired in), tank/fence points are re-tested every generation they're occupied via `control.fail.prob` (`1-control.fail.prob` = extinction probability; defaults 0.01 tank / 0.05 fence per `ms/tcz-optimisation-methods.qmd`), and fenced points additionally get their arrival counts thinned by `p.breach` (from fence length + gates) at the colonisation step, using an inflated per-point detection radius `r.eff` in place of the global scalar `r`. `stop.on.target=FALSE` disables the early-exit-on-breach behavior (needed for full-horizon point-years accounting). `rollup=TRUE` and `control.step=TRUE` together raise an error — rollup's propagule-emission shortcut is irrelevant at TCZ scale (a few hundred points) and would silently hide extinction-eligible sites from testing. See `src/tcz-sims/control-step-demo.R` for a synthetic worked example of all of this against made-up tank/fence/natural point data.
- **`run_sims()`** — repeatedly calls `spread.pilb()` (default 100 reps), resampling `lambda` from a log-normal posterior each rep. Passes through `control.step`/`stop.on.target`/`drying.matrix` to `spread.pilb()`.
- **`save_outputs()`** — writes `.RData`/`.csv` summaries (mean arrival year per point, per-year colonisation-probability tables for animation) to `out/`. `extinction.aware=TRUE` switches the arrival-year calculation from `age`-based inference (which assumes permanent occupation) to the direct `first.col.gen` record, needed once `control.step` extinction/recolonisation is in play.
- **`make_plots()`** — builds static (`tmap`) arrival-year maps and an animated GIF (`magick`) from the dynamic-map CSVs; cleans up intermediate frame files afterward.
- Dispersal kernel: `dcncross()` / `dcncross.trunc()` implement a truncated Cauchy-normal (2D t-like) kernel; `rain_to_days()` converts rainy days to days of toad movement.

Every scenario script follows the same shape: source dependencies → load point data → call `setup(...)` with dataset-specific column names → `run_sims()` → `save_outputs()` → `make_plots()`.

### Scenario/data pipeline: `src/tcz-sims/`

This is the current, actively-developed pipeline (superseding `src/tingley-sims/`, `src/pilbara-impact-sims/timing-to-tcz.R`, and other legacy folders — see below).

- **`load-data.R`** — the shared data-assembly step. Reads the Fulcrum TCZ infrastructure audit (geojson) and costing spreadsheet, joins costs to sites by fuzzy name-matching (`normalize_name()` + a manual `recode()` table for known mismatches), loads legacy LaGrange points and Living Waters natural-point layers, extracts rainfall raster values at each point, calls `score_colonised()` to flag which points are already colonised, reprojects to Albers (EPSG:3577) for metric distance work, and produces the diagnostic figures `out/time-to-tcz.png` / `out/tcz-waterpoints.pdf`.
- **`get-invasion-front.R`** — `score_colonised()`: reconstructs the current invasion-front line (fit elsewhere, stored in `invasion-front-parameters.Rdata`) in Albers coordinates, and classifies each point as colonised/not based on position relative to that line and to fixed geographic anchor points (Fitzroy River mouth, a southern boundary point).
- **`tcz-sims.R`** — the "what if the TCZ is fully implemented" scenario: buffers the TCZ boundary by 80km, defines north-of-TCZ as already colonised and south-of-TCZ as target points, and runs the model with `TCZ = TRUE` (i.e., artificial waterpoints inside the TCZ are excluded per the containment design).
- **`forecast.R`** — the "do nothing" scenario from the real current invasion front (`TCZ = FALSE`); this is the one that feeds the Shiny app (`source("shiny/prep_app_data.R")` at the end).
- **`dea-drying-probs.R`** — queries the Digital Earth Australia Waterbodies WFS/API to derive empirical dry-season drying probabilities for natural waterpoints (`living.waters`), matching by nearest polygon within `MAX_DIST_M`. Must run with `living.waters` already in session (from `load-data.R`).
- **`rainfall-drying-model.R`** — fits a `glmer` (binomial, random intercept per waterbody) relating wet-season SILO rainfall (current + lag-1) to the empirical drying classification from `dea-drying-probs.R`; exposes `predict_drying()` for scoring new rainfall records. Requires `annual_metrics` / `living.waters.dea` in session.
- **`control-step-demo.R`** — self-contained synthetic demo/regression script (no `DATA_PATH` dependency) for the control-step mechanics in `pprocess_functions.R` (see above): builds made-up tank/fence/natural point data, checks that `control.step=FALSE` reproduces legacy permanent-occupation behavior, checks tank/fence extinction rates and fence breach-thinning behave as designed, and sanity-checks the `occ.years`/`first.col.gen` bookkeeping columns. Real per-site `control.type`/`fence.area`/`fence.units` data is not yet wired into `load-data.R`/`tcz-sims.R` (see `data-requirements.md`).

Data flow through this pipeline: `load-data.R` → (`tcz-sims.R` | `forecast.R`) → `save_outputs()`/`make_plots()` → (for forecast only) `shiny/prep_app_data.R` → `shiny/app.R`. Drying-probability outputs (`dea-drying-probs.R` → `rainfall-drying-model.R`) feed into the cost/optimisation side of the model described in `ms/tcz-optimisation-methods.qmd` but are not yet wired into `setup()`/`spread.pilb()` — see `data-requirements.md` for the current gaps (fence length, non-matched-waterpoint drying probability, infrastructure cost merge) before extending this further.

### Legacy/adjacent folders

- **`src/tingley-sims/`** — original Pilbara barrier scenarios (natural-only, with/without artificial waterpoints). Superseded by `src/tcz-sims/` but kept for reference/reproducibility.
- **`src/pilbara-impact-sims/`** — `timing-to-pilbara.R`/`timing-to-tcz.R` (the latter explicitly marked deprecated in favour of `tcz-sims/forecast.R`), `timing-figure*.R` for plotting, and `spread-rate-rainfall.R` implementing **Method 2** from the paper (a rainfall-proportional spread-rate extrapolation, independent of the point-process model — see `Readme.md`).
- **`src/paruku-impact-sims/`** — same model pattern applied to the Paruku waterpoint dataset.
- **`src/ABC/`** and **`src/convolutions/`** — model-fitting/validation code: ABC rejection sampling against Kimberley/VRD invasion observations to estimate `lambda`/`r` posteriors, and the telemetry-resampling procedure (correlated random walk, `cor_angle()`/`nday()`) used to fit dispersal kernel parameters and truncation distances (see `ms/change_log.Rmd` for the rationale). These are largely one-off calibration scripts, not part of the regular run pipeline, and predate the `dat/Posteriors*.RData` / `dat/Kernel-fits_truncated.RData` artifacts they produced.
- **`src/cluster/`** — HPC environment setup for running large simulations off a local machine.

### Manuscript

`ms/tcz-optimisation-methods.qmd` is the methods writeup (Quarto → PDF) referencing the paper's three arrival-time estimation methods (see `Readme.md`); `ms/change_log.Rmd` documents the historical evolution of the model since the original Tingley et al. code (kernel truncation switch, VRD→Kimberley refit, code restructuring into `setup()`/`spread.pilb()`).

## Data conventions worth knowing

- Coordinates are carried in **Australian Albers (EPSG:3577)** internally (metric, needed for distance-based kernel calculations); raw spatial layers are loaded in WGS84 (EPSG:4326) or their native CRS and reprojected.
- `present.id`/`colonised` columns use `0` = uncolonised, `1` = already colonised, `2` = target site (converted to `0` internally by `setup()`, tracked separately via the `target` column).
- `artificial.natural.id` must resolve to a logical/0-1 flag for "natural" (`TRUE`/`1` = natural, matching `nats` in `spread.table`); scripts differ in how they derive this (`as.numeric(as.factor(...))-1`, `origin_des == "Natural"`, etc.) — check the source column's encoding before reusing a script as a template for a new scenario.
- `dat/` and `out/` are gitignored (data and generated outputs respectively); everything under them is either fetched from `DATA_PATH`/external APIs or regenerated by the scripts above.
