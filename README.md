# DATA PACKAGE: Shifting Resilience in Sagebrush Rangelands

**Trends and Predictors of Mesic Resource Productivity in Western U.S. Rangelands**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18672909.svg)](https://doi.org/10.5281/zenodo.18672909)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

**Contact:** Kristopher Mueller ([kristopher.mueller@mso.umt.edu](mailto:kristopher.mueller@mso.umt.edu))

## 1. Overview

This data package contains the code, input data, and analysis workflows used to generate the results and figures for the manuscript **"Shifting resilience: Trends and Predictors of Mesic Resource Productivity in Western U.S. Rangelands."**

The analysis combines **Google Earth Engine (GEE)** for remote sensing data extraction and **R** for Random Forest modeling and statistical visualization.

## 2. Directory Structure

The project is organized into five main directories. To keep the repository lightweight, input data (`03`) and spatial data (`05`) are hosted on Zenodo.

```text
├── 01_scripts_gee/     # JavaScript code used in the Google Earth Engine Code Editor
├── 02_scripts_r/       # R scripts for Random Forest modeling, analysis, and plotting
├── 03_data_inputs/     # (Download from Zenodo) Pre-processed GEE CSV exports
├── 04_data_outputs/    # Generated results (plots/ and stats/)
├── 05_spatial_data/    # (Download from Zenodo) Geospatial files for sampling locations
└── README.md           # This file
```

## 3. Data Access & Setup

### Step 1: Download the Code (GitHub)
1.  Click the green **<> Code** button at the top of this GitHub page.
2.  Select **Download ZIP**.
3.  Extract the file to your computer.

### Step 2: Download the Data (Zenodo)
1.  Navigate to the **Zenodo Record**: [https://doi.org/10.5281/zenodo.18672909](https://doi.org/10.5281/zenodo.18672909)
2.  Download the **Shifting_Resilience_Data_Bundle.zip**.
3.  Extract this file. It contains two folders: `03_data_inputs` and `05_spatial_data`.
4.  **Copy/Paste** these two folders into your main project folder, merging them with the existing structure.

### Step 3: System Requirements
*   **R & RStudio:** Code was developed using R version 4.x.
*   **R Packages:** `tidyverse` (dplyr, ggplot2, tidyr, stringr, readr), `randomForestSRC`, `sf`, `jsonlite`, `cowplot`, `patchwork`, `scales`, `data.table`.
*   **Google Earth Engine:** Account required only if regenerating raw input CSVs.

## 4. Instructions for Replication

### A. Data Generation (Google Earth Engine)
*Note: The outputs of these scripts are already provided in `03_data_inputs`. You only need to run these if you wish to regenerate the raw data.*

1.  Scripts `GEE01` through `GEE06` should be run sequentially in the GEE Code Editor.
2.  `GEE06_Covariate_Sampling.js` performs the stratified random sampling and covariate extraction that powers the Random Forest models.

### B. Statistical Analysis & Plotting (R)
*The working directory must be set to the `02_scripts_r` folder.*

**Step 1: Environment Setup**
*   Open RStudio and set the working directory to the `02_scripts_r` folder.

**Step 2: Generate Spatial Artifacts**
*   Run **`R10_Export_Model_Points.R`**.
*   This script applies the specific area-weighted sampling logic to the raw data and exports the `R10_Biome_Model_Points_12k.shp` to the `05_spatial_data` folder. This confirms the geographic distribution of the model training data.

**Step 3: Run Random Forest Models**
*   Run **`R05_RF_Biome_Temporal.R`** (Temporal Anomaly Model).
*   Run **`R07_RF_Biome_Spatial.R`** (Spatial Gradient Model).
*   *Warning:* These scripts involve computationally intensive Random Forest processing (1,000 trees, permutation importance). They may take significant time to complete.

**Step 4: Regional Analysis**
*   Run **`R06_RF_Ecoregion_Temporal.R`** and **`R08_RF_Ecoregion_Spatial.R`** to generate region-specific model statistics.
*   Run **`R09_RF_Ecoregion_Spatial_Comparisons.R`** to generate the comparison plot between the Great Basin and Northern Great Plains.

**Step 5: Generate Figures**
*   Run **`R01`**, **`R02`**, **`R03`**, and **`R04`** to generate the descriptive statistics and time-series plots found in the manuscript.

## 5. File Descriptions

### Google Earth Engine Scripts (`01_scripts_gee`)

| Script | Description |
| :--- | :--- |
| `GEE01_Mesic_Mask_Creation.js` | Generates the binary mask for mesic resources. |
| `GEE02_Mesic_Area_Calculation.js` | Calculates mesic area per ecoregion. |
| `GEE03_NDVI_PDSI_TimeSeries...` | Extracts annual NDVI and PDSI values. |
| `GEE04_Sample_Points_Extraction.js` | Generates 200k sampling points. |
| `GEE05_Tau_Trend_Export.js` | Calculates Mann-Kendall trends for sample points. |
| `GEE06_Covariate_Sampling.js` | Extracts environmental covariates for the Random Forest models. |

### R Scripts (`02_scripts_r`)

| Script | Description |
| :--- | :--- |
| `R01_Plot_Mesic_Areas.R` | Plots total mesic area and proportion by ecoregion. |
| `R02_Plot_Biome_NDVI_PDSI.R` | Plots biome-scale weighted NDVI/PDSI trends. |
| `R03_Plot_Ecoregion_NDVI_PDSI.R` | Plots region-scale NDVI/PDSI trends. |
| `R04_Plot_Tau_Boxplots.R` | Plots distribution of Kendall's Tau values. |
| `R05_RF_Biome_Temporal.R` | Main temporal Random Forest model (Biome scale). |
| `R06_RF_Ecoregion_Temporal.R` | Temporal Random Forest models (Ecoregion scale). |
| `R07_RF_Biome_Spatial.R` | Main spatial Random Forest model (Biome scale). |
| `R08_RF_Ecoregion_Spatial.R` | Spatial Random Forest models (Ecoregion scale). |
| `R09_RF_Ecoregion_Spatial...` | Comparative analysis (GB vs NGP). |
| `R10_Export_Biome_Model_Points.R` | Exports the exact points used in biome models to Shapefile. |

## 6. Data Dictionary

The CSV files in `03_data_inputs` contain the variables used in Random Forest modeling.

**Response Variable**
| Variable | Description |
| :--- | :--- |
| `mean_lateSeason_ndvi` | Mean Normalized Difference Vegetation Index (July 15 - Sept 30). |

**Climate & Water Covariates**
| Variable | Description |
| :--- | :--- |
| `mean_lateSeason_pdsi` | Palmer Drought Severity Index. |
| `mean_lateSeason_precip` | Cumulative precipitation (mm). |
| `mean_Spring_precip` | Spring precipitation (mm) (March-June). |
| `mean_lateSeason_maxTemp` | Maximum Temperature (C). |
| `mean_prevWinter_maxTemp` | Previous Winter Maximum Temperature (C). |
| `annual_huc8_runoff` | Annual runoff (mm) within HUC8 watershed. |
| `period_huc8_runoff` | Late-season runoff (mm) within HUC8 watershed. |
| `mean_annual_huc8_swe` | Snow Water Equivalent (mm). |

**Land Cover & Topography**
| Variable | Description |
| :--- | :--- |
| `elevation` | Elevation (meters) from USGS 3DEP. |
| `slope` | Slope (degrees). |
| `percent_tree/shrub` | Vegetation cover (%) from Rangeland Analysis Platform (RAP). |
| `agDensity_05km` | Density of irrigated agriculture (hectares within 500m). |
| `[type]Wetlands_density...` | Density of seasonal/semi-permanent/temp wetlands (hectares within 500m). |

**Other**
| Variable | Description |
| :--- | :--- |
| `global_CO2` | Annual global CO2 concentration (ppm). |
| `.geo` | GeoJSON string of point location. |
| `year` | Observation year. |
| `region` | Ecoregion identifier. |

## 7. License

This work is licensed under a [Creative Commons Attribution 4.0 International License](http://creativecommons.org/licenses/by/4.0/).