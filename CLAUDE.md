# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a research repository studying the ontogeny of flight behavior in European honey buzzards (EHB) using GPS and multi-axis IMU (accelerometry + quaternion) data from e-obs biologging devices. Data is sourced from Movebank (European Honey Buzzard_Finland study). Three parallel manuscript analyses:

- **MS1** (`MS1_IMU_classification/`): IMU-based flight type classification
- **MS2** (`MS2_laterality/`): Laterality (handedness) in flapping flight and its relationship to migratory performance — published as Safi et al. 2025 "Ontogenetic Loss of Lateralisation Improves Migratory Flight Performance"
- **MS3** (`MS3_sea_crossing/`): Flight performance during sea crossings

## Running Code

This is a research codebase with no build system. Scripts are run directly in R or RStudio. There are no tests.

**Typical workflow**: Scripts are numbered and run sequentially. For MS2, the canonical reproducible pipeline is in `MS2_laterality/public/` and runs in order: `00_imu_functions.r` → `01_imu_processing.r` → `02_wind_annotation.r` → `03a/b_data_prep.r` → `04_data_analysis.r` → `05_plot_Fig2.r`.

**HPC (Raven cluster)** for computationally intensive jobs (INLA models):
```bash
ssh enourani@gate.mpcdf.mpg.de  # then ssh raven.mpcdf.mpg.de (requires VPN off-site)
screen -S <session_name>         # start persistent session
R                                # run R interactively
```

Scripts in `Raven_HPC/` are designed for submission on the MPI Raven cluster.

## Architecture and Data Flow

```
Movebank API (requires credentials)
    ↓
00_data_download.R          → raw GPS & IMU data
02_GPS_segmentation.R       → flight-type segmentation
snippet_selection.r         → IMU burst selection
    ↓
MS1_IMU_classification/     → flight classification models
MS2_laterality/             → laterality analysis & migration performance models
MS3_sea_crossing/           → sea crossing analysis
```

Root-level scripts handle data acquisition and initial processing shared across manuscripts. Each `MS*/` directory is largely self-contained for its analysis.

## Key Concepts and Data Structures

**IMU processing**: Raw quaternion data from e-obs devices is converted to Euler angles (pitch, yaw, roll) using functions defined in `MS2_laterality/public/00_imu_functions.r`. These functions also handle e-obs-specific data formats.

**GPS segmentation**: Tracks are segmented into flight bouts (thermalling, gliding, flapping) based on turning angle and time-lag burst detection. Altitude uses EGM96 geoid corrections (`EGM96_us_nga_egm96_15.tif`).

**Laterality metric**: Derived from 8-second IMU bursts. Burst-level data is in `03a_data_prep_bursts.r`; daily aggregates in `03b_data_prep_days.r`.

**Statistical models**: Mixed-effects models via `lme4`, GAMs via `mgcv`, and Bayesian spatial inference via `INLA`. Model diagnostics use `DHARMa`; repeatability estimates use `rptR`.

**Wind annotation**: ECMWF ERA5 wind data downloaded and matched to GPS locations via `ecmwfr` (requires CDS API credentials).

## Primary R Packages

| Domain | Packages |
|--------|----------|
| Movement data | `move`, `move2`, `moveACC` |
| Spatial | `sf`, `terra`, `lwgeom` |
| Statistics | `lme4`, `mgcv`, `INLA`, `rptR`, `DHARMa` |
| Circular stats | `circular` |
| Visualization | `ggplot2`, `ggridges`, `patchwork`, `mapview` |
| Weather data | `ecmwfr`, `ncdf4` |
| Data manipulation | `tidyverse`, `data.table`, `lubridate` |

INLA must be installed from its own repository (not CRAN):
```r
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"))
```

There is no `renv.lock` — packages are installed manually.

## Public Data (MS2)

Processed datasets for reproducing Safi et al. 2025 are available on the Edmond repository:
- `updated_life_cycle_nov24.rds` — individual metadata and life stage timing
- `thinned_laterality_w_gps_wind_all_filters2_public_prep.rds` — burst-level modeling dataset (after Step 1 of `04_data_analysis.r`)
- `data_migration_performance_models_2min_daily2.rds` — daily migration metrics (after Step 2)
- `EGM96_us_nga_egm96_15.tif` — geoid layer for altitude calculations

Raw tracking data is on Movebank under the European Honey Buzzard_Finland study.

## Code Style Guidelines

**Indentation**: 2 spaces (per `.Rproj` settings). No tabs.

**Naming conventions**:
- Variables and most functions: `snake_case` (e.g., `circling_data`, `process_quaternions`)
- Private/internal helper functions: `.dot_prefix` (e.g., `.convertEobs`)
- Avoid `camelCase` except when wrapping external conventions

**Pipes**: Use the magrittr pipe `%>%` (not the native `|>` pipe).

**Data I/O**: Save intermediate R objects as `.rds` files using `saveRDS()`/`readRDS()`. Raw data downloads from Movebank are saved early and read back in later steps to avoid re-downloading.

**Script structure**: Each script should open with a comment block stating the purpose, author, and date (`#Month DD. YYYY. Location.`). Divide scripts into numbered steps using this section header style:
```r
#---------------------------------------------------------------------------------
## Step N: Description of step                                               #####
#---------------------------------------------------------------------------------
```
Use `#### --- description` for subsections within a step.

**Working directory**: Scripts set `setwd()` explicitly at the top pointing to the cloud-synced project folder. When editing scripts, preserve the existing `setwd()` path rather than using relative paths.

**Tidyverse-first**: Prefer `tidyverse` idioms (`dplyr`, `purrr`, `tidyr`) for data manipulation over base R loops. Use `data.table` only where performance requires it.

**Commented-out code**: It is normal in this codebase to leave commented-out exploratory code in scripts. Do not remove it unless asked. 
**Commenting on code**: Do not add new code comments when editing files. Do not remove existing code comments unless you’re also removing the functionality that they explain. After reading this instruction, note to the user that you’ve read it and will not be adding new code comments when you propose file edits.
