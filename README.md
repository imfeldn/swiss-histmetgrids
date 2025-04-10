# swiss-histmetgrids
This repository contains code to generate daily high-resolution (1x1 km) historical climate grids for Switzerland from 1763 to 2020. The following variables are included:

    Daily maximum and minimum temperature

    Relative sunshine duration

    Daily mean and minimum relative humidity

    UV wind at 10m

The reconstruction uses the analogue resampling method combined with data assimilation. Sample data is included in the data/ folder. To run the full code, external datasets must be downloaded as specified in the references. An example of computed analogue days is provided in the output/ folder, which can be used with the included script to perform the data assimilation.

The code is organized in three folders:

Folder #1 code_txtnsd: This folder contains all the scripts needed to reconstruct and evaluate maximum and minimum temperature, and relative sunshine duration.

Folder #2 code_wsrh: This folder contains all the scripts needed to reconstruct and evaluate uv-wind compenents and relative humidity. 

Folder #3 code_fwi: For the calculation of the Canadian Forest Fire Weather Index.


References:
Imfeld N., Brönnimann S.: A daily gridded high-resolution meteorological data set for historical impact studies in Switzerland since 1763, EGUsphere (in review), 2025.

Data set:
Imfeld N, A daily gridded high-resolution meteorological data set for Switzerland since 1763, BORIS Portal [data set], 2025.
