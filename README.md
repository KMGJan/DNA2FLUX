# DNA2FLUX

The project contains code and data.

The **code/** folder contains all the code needed to reproduce this study.

-   *SHARKphysical.sh* and *SHARKplankton.sh* download SharkWeb data from the API to data/imported/sharkweb/
-   *GetMonitoringData.R* call the two bash scripts, merge the data and save them to data/processed/shark/
-   *InterpolateWeekly.R* interpolate the processed shark data and the fish raw data (data/raw/fish_parameters.csv) and save them to data/processed/interpolation

The **data/raw/** folder contains:

- Metabarcoding data have this form <predator>_<barcode>_[asv,metadata]>.csv
- Fish parameters (biomass, abundance, bodymass) for each year and ICES statistical rectangle in fish_parameters.csv