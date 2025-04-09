# swiss-histmetgrids
This repository contains the code for creating daily high-resolution (1x1 km) historical grids for Switzerland from 1763 to 2020 for the following variables:
Daily maximum temperature, daily minimum temperature, relative sunshine duration, daily mean and minimum relative humidity, uv-wind at 10m. 
The reconstruction is based on the analogue resampling method and data assimilation. Some sample data is provided, but for running the full code, you rely on downloading data from the respective sources stated in the reference. 

The code is organized in three folders:

Folder #1 code_txtnsd: This folder contains all the scripts needed to reconstruct and evaluate maximum and minimum temperature, and relative sunshine duration.

Folder #2 code_wsrh: This folder contains all the scripts needed to reconstruct and evaluate uv-wind compenents and relative humidity. 

Folder #3 code_fwi: For the calculation of the Canadian Forest Fire Weather Index.


References:
Imfeld N., Brönnimann S.: A daily gridded high-resolution meteorological data set for historical impact studies in Switzerland since 1763, EGUsphere (in review), 2025.

Data set:
Imfeld N, A daily gridded high-resolution meteorological data set for Switzerland since 1763, BORIS Portal [data set], 2025.
