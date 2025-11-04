We're going with this datasource for NASA Precipitation data:

GPM IMERG Early Precipitation L3 1 day 0.1 degree x 0.1 degree V06 (GPM\_3IMERGDE)
[1 Day early run precipitation](https://disc.gsfc.nasa.gov/datasets/GPM_3IMERGDE_06/summary?keywords=IMERG)

Some notes about the [data format](./imerg_data_format.md)

Here are our bounds:

Time: 2012 to 2023 (11 years, or ~4020 days)

Space:
  * 10km grid of Mozambique,
  * Latitudes: -29.00 to -8.00
  * (or, 900+285=1185 to 900-274=626)
  * 626:1:1185
  * Longitudes: 30.00 to 41.00
  * (or, 1800+420=2220 to 1800-88=1712)
  * 1712:1:2220

Data can be segmented through OpenNDAP, we will use the
[xarray](https://pypi.org/project/xarray/) python library for programmatic
access. It will be visualized through the [holoviews](https://holoviews.org/)
library

Some notes on the Data Format [here](data_format.md).

You need a user account with NASA's EarthData platform in order to access this
data. You can follow the steps outlined [here](https://disc.gsfc.nasa.gov/information/documents?title=Data%20Access)
and put your username/password in the `netrc_urs` file in order to download the data.