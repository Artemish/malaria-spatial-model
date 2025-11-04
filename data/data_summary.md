This document provides summarizes the tables and data structures present in the
climate database. The database is organized along dimensions, which are:

- Mozambique Provinces
- Date
- Location
- Storms
- Weather Stations

And contains data on:

- Precipitation
- Temperature
- Malaria Data
- District Populations
- Storm Data
- Station (point) Data

The data can be accessed inside of a sqlite3 database on the MSI server. You
can copy and paste this into your terminal on MSI to get a prompt:

```sh
cd /home/ksearle/shared/downloaded-data/
sqlite3 combined_dataset.sqlite
```

Sqlite is similar to SQL, you can find a quick primer here: [sqlite cheat sheet][sqlite_ref].

## Dimensions

### Mozambique Provinces

A list of all the Provinces in Mozambique.

- **Table Name:** `mozambique_provinces`

| Field | Type | Description |
|-------|------|-------------|
| id    | INTEGER | Unique identifier for the province |
| name  | TEXT | Name of the province |

### Date

This is a simple dimension table holding all of our dates. It might come in
handy in the future if we want to assign variables to dates, etc.

- **Table Name:** `date`

| Field | Type | Description |
|-------|------|-------------|
| ymd   | INTEGER | Combined Year-Month-Day identifier |
| year  | INTEGER | Year part of the date |
| month | INTEGER | Month part of the date |
| day   | INTEGER | Day part of the date |

### Location

This is a dimension table containing one entry for every grid point in the
range from latitude (-8, -28) to longitude (28, 42). Annotated with the
Province, District, and Town that contains the center of the grid point.

Shapefiles for Mozambique's administrative districts are available on MSI,
as well as [google drive][moz_shapefiles_drive].

- **Table Name:** `location`

| Field      | Type  | Description |
|------------|-------|-------------|
| id         | INTEGER | Unique identifier for the location |
| latitude   | REAL  | Latitude coordinate |
| longitude  | REAL  | Longitude coordinate |
| province   | TEXT  | Name of the province |
| district   | TEXT  | Name of the district |
| town       | TEXT  | Name of the town |

### Storms

This table is a list of all the named and unnamed storms in the
[International Best Track Archive for Climate Stewardship (IBTrACS)][ibtracs]
dataset.

- **Table Name:** `storms`

| Field  | Type | Description |
|--------|------|-------------|
| sid    | TEXT | Unique identifier for the storm |
| numobs | INTEGER | Number of observations related to the storm |
| season | INTEGER | Year of the storm season |
| number | INTEGER | Number assignment for the storm in the season |
| name   | TEXT | Name of the storm |

### GSOD Stations

List of stations associated with the Global Summary of the Day ([GSOD][gsod]) data
collection instrument.

- **Table Name:** `gsod_stations`

| Field     | Type  | Description |
|-----------|-------|-------------|
| id        | TEXT  | Unique identifier for the station |
| name      | TEXT  | Name of the station |
| latitude  | FLOAT | Latitude coordinate of the station |
| longitude | FLOAT | Longitude coordinate of the station |
| elevation | FLOAT | Elevation of the station (in meters) |

## Data Tables

### Precipitation

The following is a summary of the data format for the NASA daily precipitation measurements
provided under [1 Day early run precipitation][GPM_3IMERGDE_06]

- **Table Name:** `precipitation`

| Field                          | Type    | Description                                                               | Units |
|--------------------------------|---------|---------------------------------------------------------------------------|-------|
| loc_id                         | INTEGER | Foreign key linking to location                                           | -     |
| date                           | TEXT    | Foreign key linking to date                                               | -     |
| HQprecipitation                | FLOAT   | Daily accumulated High Quality precipitation from all available SSMI/S... | mm    |
| HQprecipitation_cnt            | FLOAT   | Count of all valid half-hourly HQprecipitation retrievals for the day    | count |
| precipitationCal               | FLOAT   | Daily accumulated precipitation (combined microwave-IR) estimates         | mm    |
| precipitationCal_cnt           | FLOAT   | Count of all valid half-hourly precipitationCal retrievals for the day    | count |
| precipitationCal_cnt_cond      | FLOAT   | Count of valid half-hourly precipitationCal retrievals for the day        | count |
| randomError                    | FLOAT   | Daily total error of precipitation estimate                               | mm    |
| randomError_cnt                | FLOAT   | Count of valid half-hourly randomError retrievals for the day             | count |

### Temperature

- **Table Name:** `temperature`
- **Description:** Temperature-related data points.

Sourced from the [MODIS MOD11A1][modis_source] product provided by NASA. 

We currently sample this dataset to produce temperature aggregates aligned to the
10km (0.1 lat/long) grid, by windowing the Land Surface Temperature variables
for Day and Night.

The Temperature dataset originates from the Terra MODIS Land Surface
Temperature/Emissivity Daily (MOD11A1) Version 6.1 product. It offers daily
readings of Land Surface Temperature and Emissivity on a 1 km spatial
resolution, spread over a 1,200 by 1,200 km grid. The dataset also incorporates
daytime and nighttime temperature bands, complemented by quality assessments,
observation timings, view angles, clear-sky metrics, and emissivities from
bands 31 and 32 tied to specific land covers.

| Field             | Type  | Description                                                  | Units |
|-------------------|-------|--------------------------------------------------------------|-------|
| loc_id            | INTEGER | Foreign key linking to location                            | -     |
| date              | TEXT  | Foreign key linking to date                                  | -     |
| Clear_day_cov     | REAL  | Day clear-sky coverage                                       | -     |
| Clear_night_cov   | REAL  | Night clear-sky coverage                                     | -     |
| Day_view_angl     | REAL  | View zenith angle of daytime Land-surface Temperature        | deg   |
| Day_view_time     | REAL  | Time of daytime Land-surface Temperature observation         | hrs   |
| Emis_31           | REAL  | Band 31 emissivity                                           | -     |
| Emis_32           | REAL  | Band 32 emissivity                                           | -     |
| LST_Day_1km       | REAL  | Daily daytime 1km grid Land-surface Temperature              | K     |
| LST_Night_1km     | REAL  | Daily nighttime 1km grid Land-surface Temperature            | K     |
| Night_view_angl   | REAL  | View zenith angle of nighttime Land-surface Temperature      | deg   |
| Night_view_time   | REAL  | Time of nighttime Land-surface Temperature observation       | hrs   |
| QC_Day            | REAL  | Quality control for daytime LST and emissivity               | -     |
| QC_Night          | REAL  | Quality control for nighttime LST and emissivity             | -     |

### Malaria Data

This data was provided to us by the Mozambique Ministry of Health, and
summarizes Malaria case data from 2015 to present. We keep a CSV copy
of this data in Google Drive: [Malaria Case Data][malaria_case_data_drive]

- **Table Name:** `malaria_data`

| Field          | Type  | Description |
|----------------|-------|-------------|
| province      | TEXT  | Foreign key linking to mozambique_provinces |
| year          | INTEGER | Year of the data |
| month         | INTEGER | Month of the data |
| cases         | INTEGER | Number of malaria cases |
| hospitalizations | INTEGER | Number of hospitalizations due to malaria |
| deaths        | INTEGER | Number of deaths due to malaria |

### District Populations

The following data was summarized from the Malaria Dataset mentioned earlier.

- **Table Name:** `province_populations`
- **Description:** Population data for each province by year.

| Field      | Type  | Description |
|------------|-------|-------------|
| province   | TEXT  | Foreign key linking to mozambique_provinces |
| year       | INTEGER | Year of the data |
| population | INTEGER | Population count |

### Storm Data

Storm observations are instantaneous records of wind, pressure, and maximum
gusts, organized by latitude/longitude, time, and storm ID. 

- **Table Name:** `storm_data`
- **Description:** Data points related to storm observations.


| Field                 | Type  | Long Name                                                  | Description                                                                                              | Units     | Content Type                    |
|-----------------------|-------|------------------------------------------------------------|----------------------------------------------------------------------------------------------------------|-----------|-----------------------------------|
| obs_id                | INTEGER | -                                                         | Unique identifier for the storm observation                                                            | -         | -                               |
| sid                   | TEXT  | -                                                         | Foreign key linking to storms                                                                            | -         | -                               |
| date                  | INTEGER | -                                                         | Foreign key linking to date                                                                            | -         | -                               |
| lat                   | REAL  | -                                                         | Latitude of the observation                                                                              | -         | -                               |
| lon                   | REAL  | -                                                         | Longitude of the observation                                                                             | -         | -                               |
| time                  | TEXT  | Time (ISO)                                                | Nominally, time steps are 3 hourly, but can be more often since some agencies include extra position (e.g., times near landfall, maximum intensity, etc.) | -         | physicalMeasurement             |
| wmo_wind              | REAL  | Maximum sustained wind speed from Official WMO agency     | -                                                                                                        | kts       | physicalMeasurement             |
| wmo_pres              | REAL  | Minimum central pressure from Official WMO agency         | -                                                                                                        | mb        | physicalMeasurement             |
| wmo_agency            | TEXT  | Official WMO agency                                       | -                                                                                                        | -         | thematicClassification          |
| storm_speed           | REAL  | Storm translation speed                                   | -                                                                                                        | kts       | physicalMeasurement             |
| storm_dir             | REAL  | Storm translation direction                               | Unit is degrees east of north                                                                            | degrees   | physicalMeasurement             |
| combined_wind         | REAL  | Combined wind speed                                       | Combined estimate of wind speed across all sources                                                       | kts      |                           |
| combined_pres         | REAL  | Combined atmospheric pressure                             | Combined estimate of atmospheric pressure across all sources                                             | mb      |                           |
| combined_rmw          | REAL  | Combined radius of max winds                              | Combined estimate of radius of max winds across all sources                                              | km      |                           |
| combined_gust         | REAL  | Combined max gust speed                                   | Combined estimate of max gust speed across all sources                                                   | kts      |                           |


### GSOD Measurements

Contains daily weather measurement data from [GSOD][gsod] stations.

- **Table Name:** `gsod_measurements`

| Field    | Type  | Description |
|----------|-------|-------------|
| station  | TEXT  | Foreign key linking to `gsod_stations` by its `id` |
| date     | TEXT  | Date of the measurement |
| temp     | REAL  | Temperature measurement for the day |
| dewp     | REAL  | Dew point measurement |
| slp      | REAL  | Sea level pressure measurement |
| stp      | REAL  | Station pressure measurement |
| visib    | REAL  | Visibility measurement |
| wdsp     | REAL  | Wind speed measurement |
| mxspd    | REAL  | Maximum wind speed recorded for the day |
| gust     | REAL  | Gust speed measurement |
| max      | REAL  | Maximum temperature recorded for the day |
| min      | REAL  | Minimum temperature recorded for the day |
| prcp     | REAL  | Precipitation measurement |


[gsod]: https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/ec4af8ce752a4898b308687c8c18351f/html
[sqlite_ref]: https://www.sqlitetutorial.net/sqlite-cheat-sheet/
[ibtracs]: https://www.ncei.noaa.gov/products/international-best-track-archive
[moz_shapefiles_drive]: https://drive.google.com/drive/u/2/folders/1Zu3GPjpWQZaPwGI4ZQyA6gT4qafAjFbB
[malaria_case_data_drive]: https://docs.google.com/spreadsheets/d/1ir3fL3tCK7ICLZ2iUIDsfUI72UN7R0fzE9NAziYjSJI/edit
[modis_source]: https://lpdaac.usgs.gov/products/mod11a1v061/
[GPM_3IMERGDE_06]: https://disc.gsfc.nasa.gov/datasets/GPM_3IMERGDE_06/summary?keywords=IMERG
