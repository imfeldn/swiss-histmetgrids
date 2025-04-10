#!/usr/bin/env python
# coding: utf-8
# -----------------------------------------------
# validation of wind speed
# 08.04.2024
# -----------------------------------------------

# import packages 
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from natsort import natsorted
import matplotlib.patches as mpatches
import warnings  # to omit warnings for spatial correlation (as division by zero (no rain days) is not possible)

import glob
import time
import os

#%% define new functions because of nan values
rmse  = lambda x,y: np.sqrt(np.nanmean((x-y)**2, axis = -1))
bias  = lambda x,y: np.nanmean(y,axis = -1) - np.nanmean(x, axis = -1) 

def rmse_calculation(x, y, dim):
    return xr.apply_ufunc(
        rmse, x, y,
        input_core_dims=[[dim], [dim]],
        dask='allowed',
        vectorize=True,
        output_dtypes=[float])

def bias_calculation(x, y, dim):
    return xr.apply_ufunc(
        bias, x, y,
        input_core_dims=[[dim], [dim]],
        dask='allowed',
        vectorize=True,
        output_dtypes=[float])


#%% define input for evaluations
scen = "scen1"
window = "45"
arm = False
date = ""
date2 = ""
EnKF = "yourEnKFdefs"
add = ""
temp = "tmin" # define tmin or tmax

if temp == "tmax":
    varnam = "TmaxD"
else:
    varnam = "TminD"

if arm:
    temp_recfolder = "output/ARM_run_" + date + "/" + scen + "_" + window  + add + "/ARM_" + varnam + "_" + date2 + "/"
    varnam2 = varnam

else:
    temp_recfolder = "output/ARM_run_" + date + "/" + scen + "_" + window + add + "/EnKF_" + temp + "_" + EnKF + "_" + date2 + "/"
    varnam2 = temp


print("start")
# omit warnings (see above)
warnings.filterwarnings("ignore")

# create validation output folder
outdir = temp_recfolder + "/validation_" + time.strftime("%Y-%m-%d") + "/"

print(outdir)
if not os.path.isdir(outdir):
    try:
        os.mkdir(outdir)
    except OSError:
        print("")
    else:
        print("")

#%% orig data
print("read observation data")
tmax_obsfolder = "/data/" + temp +"/" ## tmax data from MeteoSwiss detrended with ERA5 and deseasonlized
tmax_obsfiles = glob.glob(f'{tmax_obsfolder}*{varnam}*detrended_deseas*1971*.nc')
tmax_obs = xr.open_mfdataset(tmax_obsfiles).chunk({"time": -1, "E": "auto", "N":"auto"})
tmax_obs["time"] = pd.to_datetime(tmax_obs.time.dt.strftime("%Y-%m-%d").values)
tmax_obs = tmax_obs[varnam]

#%% get files
print("read reconstruction data")
tmax_recfiles = glob.glob(f'{temp_recfolder}/*_{temp}*deseas*200*.nc')

if arm:
    tmax_rec = xr.open_mfdataset(tmax_recfiles).chunk({"time": -1, "E": "auto", "N":"auto"})
    tmax_rec["time"] = pd.to_datetime(tmax_rec.time.dt.strftime("%Y-%m-%d").values)
    tmax_rec = tmax_rec.reindex(time = tmax_obs["time"])[varnam2]
else:     
    tmax_rec = xr.open_mfdataset(tmax_recfiles).chunk({"time": -1, "easting": "auto", "northing":"auto"})
    tmax_rec["time"] = pd.to_datetime(tmax_rec.time.dt.strftime("%Y-%m-%d").values)
    tmax_rec = tmax_rec.rename({'easting': 'E','northing': 'N'})
    tmax_rec = tmax_rec.reindex(N = tmax_obs["N"])
    tmax_rec = tmax_rec[varnam2]

tmax_obs = tmax_obs.reindex(time = tmax_rec["time"])

#%% define seasons
seasons = {"JJA":np.arange(6,9),"DJF": [12,1,2],"MAM": np.arange(3,6),"SON":np.arange(9,12)}

print(seasons)

for name,value in seasons.items():
    print(value)
    print(name)
    
    obs_seas = tmax_obs.sel(time=np.isin(tmax_obs.time.dt.month,value))
    rec_seas = tmax_rec.sel(time=np.isin(tmax_rec.time.dt.month,value))
    
    qt_dims = ("time")
    qt_values = np.arange(0.1,1.1,0.1)

    ds_qt = obs_seas.quantile(qt_values, dim=qt_dims)
    ## reorder to vecotr and remove NA values 
    rmse_out = []
    bias_out = []

    for qt in np.arange(0,len(qt_values)):
        print(qt)
        if qt == 0:
            selqs_obs = obs_seas.where(obs_seas < ds_qt[qt,:,:]).stack(points = ("E","N"))
            selqs_rec = rec_seas.where(obs_seas < ds_qt[qt,:,:]).stack(points = ("E","N"))
        else:
            selqs_obs = obs_seas.where((obs_seas < ds_qt[qt,:,:]) & (obs_seas > ds_qt[qt - 1,:,:])).stack(points = ("E","N"))
            selqs_rec = rec_seas.where((obs_seas < ds_qt[qt,:,:]) & (obs_seas > ds_qt[qt - 1,:,:])).stack(points = ("E","N"))
        
        print("rmse")
        rmse_out.append(rmse_calculation(selqs_obs, selqs_rec, "points").compute().to_numpy())   
        print("bias")
        bias_out.append(bias_calculation(selqs_obs, selqs_rec, "points").compute().to_numpy())
       
    bias_mat = np.array(bias_out)
    rmse_mat = np.array(rmse_out)
    
    np.savetxt(outdir + "bias_byquantile_" + name + ".csv",  bias_mat,  delimiter = ",")
    np.savetxt(outdir + "rmse_byquantile_" + name + ".csv",  rmse_mat,  delimiter = ",")

