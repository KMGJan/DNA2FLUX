# DNA2FLUX

The project contains code and data.

# **code/**
This folder contains all the code needed to reproduce this study.

-   `InstallPacakges.R` contains all the needed libraries, it is called in CleanAnalyses.R
-   `SHARKphysical.sh` and `SHARKplankton.sh` download SharkWeb data from the API to data/imported/sharkweb/
-   `GetMonitoringData.R` calls the two bash scripts, merge the data and save them to data/processed/shark/
-   `InterpolateWeekly.R` interpolates the processed shark data and the fish raw data (data/raw/fish_parameters.csv) and save them to data/processed/interpolation
-   `CombineMetabarcodingData.R` combines all metabarcoding data into data/processed/predator_selectivity.csv
-   `CalculateFluxes.R` contains helper functions to calculate the energy fluxes between each predator-prey interactions
-   `ProcessFluxes.R` calculates the energy fluxes throughout the timeseries for all bootstrap iterations and summarise the fluxes with 95% confidence intervals
-   `HelpMyAnalyses.R` contains helper functions to analyse the data
-   `CleanAnalyses.R` contains the data analyses and visualisation
-   `PlotDataMap.R` downloads the ICES data and plot the sampling map



## **data/raw/**


### Metabarcoding

Data files are named using the format: `predator_barcode_[asv|metadata].csv`

#### ASV Tables

- **`fish_18s_asv.csv`** and **`fish_coi_asv.csv`** contain ASV tables from [Jan et al. (2025)](https://doi.org/10.1093/icesjms/fsaf122)[^1].  
  - Each **row** represents a unique ASV.  
  - Columns **1–9** (`18s`) or **1–8** (`coi`) contain taxonomic annotations.  
  - Column **10** (`18s`) or **9** (`coi`) contains the ASV identifier.  
  - All subsequent columns represent unique samples, with cell values corresponding to the number of reads for each ASV in that sample.

- **`zooplankton_16s_asv.csv`** and **`zooplankton_18s_asv.csv`** contain ASV tables from [Zamora-Terol et al. (2020)](https://doi.org/10.1111/mec.15555)[^2], [Novotny et al. (2021)](https://doi.org/10.1098/rspb.2021.0908)[^3], and [Serandour et al. (2023)](https://doi.org/10.1093/plankt/fbad007)[^4].  
  - Each **row** represents a unique ASV.  
  - Columns **1–8** (`18s`) or **1–10** (`16s`) contain taxonomic annotations.  
  - All subsequent columns represent unique samples, with cell values corresponding to read counts for each ASV.

#### Sample Metadata

Each corresponding metadata file contains one row per sample.  
The **row name** matches the **column names** in the respective ASV table, ensuring direct linkage between sample data and metadata.

| Field Name  | Description |
|---|-------------|
| `library_ID` | Unique identifier for each sequencing library (matches ASV table column). |
| `title` | Short descriptive title for the sample. |
| `organism`   | Taxonomic name or environmental sample descriptor (e.g., "Clupea harengus"). |
| `collection_date` | Date when the sample was collected (YYYY-MM-DD). |
| `geo_loc_name` | Geographic location, correspond to the ICES statistical rectangle or to the monitoring station name. |
| `depth` | Sampling depth interval in meters. |
| `samp_size` | Size or volume of the sample collected (e.g., 1 L, 1 individual). |
| `size_frac` | Size fraction (length for fish species, mesh size for WP2 samples, filter size for water samples) |
| `lat_lon` | Latitude and longitude in decimal degrees. |
| `design description` | Short description of the sampling design or context. |
| `env_broad_scale` | Broad environmental context (e.g., "Pelagic Baltic Sea"). |
| `env_local_scale` | More specific local context (e.g., "ICES statistical rectangle 44G7"). |
| `env_medium` | Type of environmental material sampled (e.g., "Seawater"). |
| `collection_method` | Description or name of the sampling protocol or method. |
| `samp_collection_device` | Equipment or device used to collect the sample (e.g., "Niskin bottle", "Pelagic trawl", "WP2 , 90um"). |
| `samp_mat_process` |  Description of any material processing steps (e.g., "Bulk DNA extraction, QIAmp Micro Kit"). |
| `source_material_ID` | Identifier linking to a parent or source material, if applicable. |
| `sample_name` | Unique sample name used in the study. |
| `sample_accession` | ENA BioSample accession number. |
| `study_accession` | Associated ENA BioProject or study accession number. |

### Fish parameters

Biomass, abundance, bodymass for each year and ICES statistical rectangle in `fish_parameters.csv`, this dataset is a cleaned version of the data in `BIAS Herr Sprat` and `BIAS Stickleback`

| Field Name  | Description |
|----|-------------|
| `node_name` | Fish genus name. |
| `ICES_rect` | ICES statisitcal rectangle. |
| `abundance` | Abundance of the fish in ind/m2. |
| `biomass` | Biomass of the fish in g/m2. |
| `bodymass` | Average bodymass of the fish in g/ind. |
| `year` | Year of the survey (YYYY). |

### Zooplankton bodymass

For each taxa based on season, station, sex and development stage according to the HELCOM Combine manual in `zooplankton_bodymass.csv`

| Field Name  | Description |
|----|-------------|
| `node_name` | Zooplankton genus name. |
| `taxon_genus` | Zooplankton genus. |
| `taxon_species` | Zooplankton species. |
| `dev_stage_code` | Zooplankton development stage code. |
| `sex_code` | Zooplankton sex code. |
| `station_name` | Sampling station name. |
| `season` | Sampling season (e.g., "spring", "fall", "summer", "winter"). |
| `bodymass` | Average bodymass of the zooplankton. |
| `bodymass_unit` | Zooplankton bodymass unit (g). |

### Zooplankton count data

From May 2022 ([Jan et al. 2025](https://doi.org/10.1093/icesjms/fsaf122)) are in `count_spras2022.csv`

| Field Name  | Description |
|----|-------------|
| `trawl_id` | trawl id corresponding to the DNA metabarcoding metadata. |
| `node_name` | Zooplankton genus name. |
| `sample_week` | Date of the Monday of the sampling week (YYYY-MM-DD). |
| `abundance` | Zooplankton abundance in ind/L. |

### Nodes information

This dataset contains taxonomy and parameters for estimating metabolic rates in `node_data.csv`

| Field Name  | Description |
|----|-------------|
| `node_name` | Node taxa. |
| `type` | Correspond to the trophic level (e.g., "phytoplankton", "zooplankton" or "fish"). |
| `16s_rrna_name` | Corresponding taxa name in the 16s metabarcoding datasets. |
| `18s_rrna_name` | Corresponding taxa name in the 18s metabarcoding datasets. |
| `dyntaxa_name` | Corresponding taxa name in the monitoring datasets. |
| `tax_level` | Taxonomy resolution of the node. |
| `efficiencies` | Assimilation efficiencies used for flux estimations. |
| `intercept` | Intercept used for metabolic rate calculation. |
| `slope` | Slope used for metabolic rate calculation. |
| `energy_density` | Node energy density, this is not used. |
| `trophic_level` | Trophic level in numeric. |
| `horizontal_position` | Helper numeric variable to plot the trophic network. |
| `color` | Node color used throughout the analyses. |

[^1]: Jan KMG, Hentati-Sundberg J, Larson N, Winder M. 2025. *Limited resource use overlaps among small pelagic fish species in the central Baltic Sea*. ICES Journal of Marine Science. https://doi.org/10.1093/icesjms/fsaf12
[^2]: Zamora-Terol S, Novotny A, Winder M. 2020. *Reconstructing marine plankton food web interactions using DNA metabarcoding*. Molecular Ecology. https://doi.org/10.1111/mec.15555
[^3]: Novotny A, Zamora-Terol S, Winder M. 2021. *DNA metabarcoding reveals trophic niche diversity of micro and mesozooplankton species*. Proceedings of the Royal Society B. https://doi.org/10.1098/rspb.2021.0908
[^4]: Serandour B, Jan KMG, Novotny A, Winder M. 2023. *Opportunistic vs selective feeding strategies of zooplankton under changing environmental conditions*. Journal of Plankton Research. https://doi.org/10.1093/plankt/fbad007