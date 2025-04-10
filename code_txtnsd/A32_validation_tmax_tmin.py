#!/usr/bin/env python
# coding: utf-8
# validation of tmax

import numpy as np
import pandas as pd
import xarray as xr

import glob
import sys
import time
import warnings  # to omit warnings for spatial correlation (as division by zero (no rain days) is not possible)
import os
from dask.distributed import Client

np.seterr(invalid='ignore')

## source functions
exec(open("A30_validation_functions.py").read())

### define evaluation path in console
temp = "tmax"
scen = "scen1_45"
date = ""
date2 = ""
arm = False
EnKF = "yourENKFdefs"

print(arm)

if temp == "tmax":
    varnam = "TmaxD"
else:
    varnam = "TminD"

if arm:
    tmax_recfolder = "output/ARM_run_" + date + "/" + scen + "/ARM_" + varnam + "_" + date +"/"
    varnam2 = varnam

else:
    tmax_recfolder = "output/ARM_run_" + date + "/" + scen + "/EnKF_" + temp + "_" + EnKF + "_" + date2 + "/"
    varnam2 = temp

#%% omit warnings (see above)
warnings.filterwarnings("ignore")

# get start time
stm = time.time()

# create validation output folder
indir = os.path.dirname(tmax_recfolder)
outdir = indir + "/validation_" + time.strftime("%Y-%m-%d") + "/"
print(outdir)
if not os.path.isdir(outdir):
    try:
        os.mkdir(outdir)
    except OSError:
        print("")
    else:
        print("")
 
# --------------------------------------------------------------------------------------------------------------------
# READ DATA
# --------------------------------------------------------------------------------------------------------------------
print("read observation data")
tmax_obsfolder = "data/" + temp +"/" ## tmax data from MeteoSwiss
tmax_obsfiles = glob.glob(f'{tmax_obsfolder}*{varnam}*detrended_deseas*1971*.nc')
tmax_obs = xr.open_mfdataset(tmax_obsfiles).chunk({"time": -1, "E": "auto", "N":"auto"})
tmax_obs["time"] = pd.to_datetime(tmax_obs.time.dt.strftime("%Y-%m-%d").values)
tmax_obs = tmax_obs[varnam]

# get files
print("read reconstruction data")
tmax_recfiles = glob.glob(f'{tmax_recfolder}/*_{temp}*deseas*.nc')

if arm:
    tmax_rec = xr.open_mfdataset(tmax_recfiles).chunk({"time": -1, "E": "auto", "N":"auto"})
    tmax_rec["time"] = pd.to_datetime(tmax_rec.time.dt.strftime("%Y-%m-%d").values)
    tmax_rec = tmax_rec[varnam2]
else:     
    tmax_rec = xr.open_mfdataset(tmax_recfiles).chunk({"time": -1, "easting": "auto", "northing":"auto"})
    tmax_rec["time"] = pd.to_datetime(tmax_rec.time.dt.strftime("%Y-%m-%d").values)
    tmax_rec = tmax_rec.rename({'easting': 'E','northing': 'N'})
    tmax_rec = tmax_rec.reindex(N = tmax_obs["N"])
    tmax_rec = tmax_rec[varnam2]

#%% select the same dates
if len(tmax_obs.time.values) > len(tmax_rec.time.values):
    tmax_obs = tmax_obs.reindex(time = tmax_rec["time"])
else:
    tmax_rec = tmax_rec.reindex(time = tmax_obs["time"]) 

#%% define seasons
seasons = {"MAM": np.arange(3,6),"SON":np.arange(9,12),"JJA":np.arange(6,9),"DJF": [12,1,2]}
print(seasons)

for name,value in seasons.items():
    print(value)
    print(name)
    
    obs_seas = tmax_obs.sel(time=np.isin(tmax_obs.time.dt.month,value))
    rec_seas = tmax_rec.sel(time=np.isin(tmax_rec.time.dt.month,value))
    print("pcorr")
    pcorr_out = pearson_correlation(obs_seas, rec_seas, "time").compute()
    print("rmse")
    rmse_out = rmse_calculation(obs_seas, rec_seas, "time").compute()
    print("msess")
    msess_out = msess_calculation(obs_seas, rec_seas, "time").compute()
    print("bias")
    bias_out = bias_calculation(obs_seas, rec_seas, "time").compute()
         
    metrics = xr.concat([pcorr_out, rmse_out,msess_out,bias_out], dim ="metrics", fill_value = np.nan)
    newcoords = ["corr_" + name, "rmse_"  + name,"msess_" + name ,"bias_"  + name]
    metrics = metrics.assign_coords(metrics=newcoords)

    print("write netcdf")
    metrics.to_netcdf(outdir + time.strftime("%Y-%m-%d") + '_validation_maps_' + temp + '_' + name +'.nc')
    del metrics, newcoords

print("eval done for " + temp)
    
    
