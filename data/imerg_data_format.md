## NASA IMERG Daily Precipitation data format

The following is a summary of the data format for the NASA daily precipitation measurements
provided under [1 Day early run precipitation](https://disc.gsfc.nasa.gov/datasets/GPM_3IMERGDE_06/summary?keywords=IMERG)

### GridHeader (Metadata):
GridHeader contains metadata defining the grids in the grid structure. See Metadata for
GPM Products for details.

### lon (4-byte float, array size: lon):
Longitude of the center of the grid box. Values range from -180 to 180 degrees east.
Special values are defined as:
-9999.9 Missing value

### lat (4-byte float, array size: lat):
Latitude of the center of the grid box. Values range from -90 to 90 degrees north. Special
values are defined as:
-9999.9 Missing value

### precipitation (4-byte float, array size: lat x lon):
Precipitation estimate using gauge calibration over land. Values range from 0 to 1000
mm/hr. Special values are defined as:
-9999.9 Missing value

### randomError (4-byte float, array size: lat x lon):
Random error estimate of precipitation. Values range from 0 to 1000 mm/hr. Special
values are defined as:
-9999.9 Missing value

### gaugeRelativeWeighting (2-byte integer, array size: lat x lon):
Surface gauge weighting in the final precipitation estimate. The values range from 0 to
100, where 0 is no gauge weighting and 100 is entirely based on gauge. Values range from
0 to 100 percent. Special values are defined as:
-9999 Missing value

### probabilityLiquidPrecipitation (2-byte integer, array size: lat x lon):
Probability of liquid precipitation. 0=definitely frozen. 100=definitely liquid. 50=equal
probability frozen or liquid. This field is globally complete and provided irrespective of
the presence of precipitation. Values range from 0 to 100 percent.

### precipitationQualityIndex (4-byte float, array size: lat x lon):
Estimated quality of precipitationCal where 0 is worse and 100 is better. Values range
from 0 to 100. Special values are defined as:
-9999.9 Missing value