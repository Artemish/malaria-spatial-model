# Model Data Files

This directory contains the processed data files used in the climate-malaria spatiotemporal modeling analysis. All files are in R data format (.rda) and represent aggregated, district-level data prepared for the Bayesian spatiotemporal model.

## Data Files Overview

### `full_m_district_PRCP.rda`
**Precipitation Data**
- **Source**: NASA GPM 3IMERGEDE_06 satellite data
- **Description**: Monthly precipitation measurements averaged per district using administrative shapefiles
- **Temporal Coverage**: 60 months of data
- **Spatial Coverage**: All 130 districts in Mozambique
- **Units**: mm of precipitation per district per month
- **Variables**: 
  - `PRCP_timetotal_spatialmean_district`: Spatial mean of temporal total precipitation
  - `PRCP_timemean_spatialtotal_district`: Temporal mean of spatial total precipitation

### `temperature_district_final.rda`
**Temperature Data**
- **Source**: MODIS MOD11A1 satellite Land Surface Temperature (LST) data
- **Description**: District-level temperature summary statistics for each month
- **Variables**: 
  - `lst_day_c_max`: Maximum daily land surface temperature (°C)
  - `lst_day_c_min`: Minimum daily land surface temperature (°C)
  - `lst_day_c_mean`: Mean daily land surface temperature (°C)
- **Spatial Coverage**: All 130 districts in Mozambique
- **Temporal Coverage**: 60 months of data

### `malaria_related.rda`
**Malaria Epidemiological Data**
- **Source**: Mozambique National Malaria Control Programme (NMCP)
- **Description**: Monthly confirmed malaria case counts by district
- **Variables**:
  - `marlaria_dst_case_final_new`: Malaria case counts per district per month
  - `marlaria_dst_population_final_new`: District population data for rate calculations
- **Temporal Coverage**: 60 months of data
- **Spatial Coverage**: All 130 districts in Mozambique

### `storm_district_final.rda`
**Storm Event Data**
- **Description**: Binary indicators of extreme precipitation events ("storm events") by district and month
- **Source**: Derived from precipitation data using threshold-based definitions
- **Format**: District-month pairs indicating presence of storm events
- **Use**: For modeling interaction effects between storms and precipitation-malaria relationships

### `map_mat.rda`
**Spatial Adjacency Matrix**
- **Description**: Binary adjacency matrix defining spatial neighborhood structure for the 130 districts
- **Format**: 130 x 130 matrix where entries indicate spatial adjacency between districts
- **Use**: Enables spatial autocorrelation modeling in the CAR (Conditional Autoregressive) framework
- **Construction**: Based on administrative district shapefiles with shared boundary definitions

### `id.rda`
**Regional Stratification Indices**
- **Description**: District identification vectors for regional analysis
- **Variables**:
  - `South_id`: Indices for southern region districts
  - `Central_id`: Indices for central region districts  
  - `North_id`: Indices for northern region districts
- **Use**: Enables region-specific model fitting and analysis (current analysis focuses on southern districts)

## Data Processing Notes

- All climate variables are centered (mean-centered) for improved model convergence
- Temporal lags (1-4 months) are constructed for both precipitation and temperature variables
- Population data serves as offset terms for modeling incidence rates
- Storm events are integrated as binary interaction terms with precipitation
- Data spans 60 months with the first 4 months used for lag construction, resulting in 56 months of complete data for analysis

[gsod]: https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/ec4af8ce752a4898b308687c8c18351f/html
[sqlite_ref]: https://www.sqlitetutorial.net/sqlite-cheat-sheet/
[ibtracs]: https://www.ncei.noaa.gov/products/international-best-track-archive
[moz_shapefiles_drive]: https://drive.google.com/drive/u/2/folders/1Zu3GPjpWQZaPwGI4ZQyA6gT4qafAjFbB
[malaria_case_data_drive]: https://docs.google.com/spreadsheets/d/1ir3fL3tCK7ICLZ2iUIDsfUI72UN7R0fzE9NAziYjSJI/edit
[modis_source]: https://lpdaac.usgs.gov/products/mod11a1v061/
[GPM_3IMERGDE_06]: https://disc.gsfc.nasa.gov/datasets/GPM_3IMERGDE_06/summary?keywords=IMERG