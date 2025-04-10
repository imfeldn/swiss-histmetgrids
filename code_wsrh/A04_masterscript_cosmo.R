###  Analogues masterscript

rm(list=ls())
RhpcBLASctl::blas_set_num_threads(1)

suppressWarnings(library(StatMatch))
suppressWarnings(library(hydroGOF))  
suppressMessages(library(ncdf4))     
suppressMessages(library(lubridate))  
suppressMessages(library(doParallel))
suppressWarnings(library(purrr))      
suppressWarnings(library(XLConnect))  

print("packages loaded")

date = "2025-04-02"
start = "1961-01-01"
end = "2020-12-31"
scen = "scen_hist2_wind"
window = 80
dist = "Gower"
vars = c("ta","p","rr0","spm","sp")
stnval = "all"
calonly = "F"
climoff = "T"
startCal = start
endCal = end


args <- commandArgs(TRUE)
print(args)
scen <- args[1]
stnval <- args[2]
dist <- args[3]
startCal <- args[4]
endCal <- args[5]
calonly <- args[6]
climoff <- args[7]
date <- args[8]


###################################################################################################
### ANALOGUE METHOD
###################################################################################################

analogs <- T
source("code_wsrh/A02_scenarios.R")

# Calibration of the analog method (list)
# pass a selection of dates because cosmo is missing for 2019-09 to 2019-12
nc <- nc_open("data/cosmo_20162020_daymean_temp_LV95.nc")
selDates <- as.Date(ncvar_get(nc,"time")/86400, origin = "1970-01-01")

specifications <- list(
    variables = vars,     # "temp", "press", "prec", "precipitation occurrence"
    scenario = scen,                        # predefine a network in 011_scenarios.R based on variables and station for which the script should be run. 
    distance = dist,                        # "Euclid","RMSE","Mahalanobis","Gower". RMSE and weighted Euclid are the same!
    spatpress = TRUE,                       # use spatial pressure differences
    nAnalog = 50,                           # number of analogs
    startRef = "2016-04-01",                # start of reference period
    endRef = "2020-09-30",                  # end of reference period
    selDates = selDates,                    # in case of missing data, specify the dates
    startCal = startCal,                    # calibration start year
    endCal = endCal,                        # calibration end year == this is end of historical period, i.e.scenario, not end of entire calibration
    calonly = ifelse(calonly=="T",T,F),     # do not include predictor dates in final analog files
    endHist = "1863-12-31",                 # historical period end year
    dayselection = window,                  # analog dates only +-dayselection days from target date  
    standardize = F,                        # standardize data?
    trend = T,                              # correct temperature trend. This needs a file with the temperature trends already calculated.  
    offset = ifelse(climoff=="T",T,F),      # remove an climate offset of station data in reference period? Needs a file with the offset calculated.
    seasonality = T,                        # correct temperature seasonality
    sep_ref_hist = F,                       # if seasonality = T, should seasonality be removed separately for ref and hist?
    window=5,                               # number of days in reference period not to consider in analog pool for cross-validation
    rr0crit=""                              # should the preselection of analogs include as well the precipitation occurrence. If yes, which would be the threshold then? eg.0.5
)              


### load helpfunction
source("code_wsrh/A03_helpfun.R")

###################################################################################################
### ANALOGUE RESAMPLING METHOD
###################################################################################################

### execute preprocessing according to specifications
source("code_wsrh/A05_preprocessing_cosmo_fulldetrending.R")

## calculate analogues
source("code_wsrh/A06_analogs.R")

