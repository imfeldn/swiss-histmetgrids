#!/usr/bin/env python
# coding: utf-8
# validation of srel

# import packages 
import numpy as np
import pandas as pd
import xarray as xr

import glob
import sys
import time
import warnings  # to omit warnings for spatial correlation (as division by zero (no rain days) is not possible)
import os

## define simple functions
rmse  = lambda x,y: np.sqrt(np.mean((x-y)**2, axis = -1))
bias  = lambda x,y: y.mean(axis = -1) - x.mean(axis = -1) 
msess2 = lambda x,y,clim: 1 - np.sum((x-y)**2, axis = -1)/np.sum((clim - x)**2, axis = -1)

def covariance_gufunc(x, y):
    return ((x - x.mean(axis = -1, keepdims = True))
            * (y - y.mean(axis = -1, keepdims = True))).mean(axis = -1)

def pearson_correlation_gufunc(x, y):
    return covariance_gufunc(x, y) / (x.std(axis=-1) * y.std(axis=-1))

## define wrapper functions for apply_ufunc
def pearson_correlation(x, y, dim):
    return xr.apply_ufunc(
        pearson_correlation_gufunc, x, y,
        input_core_dims=[[dim], [dim]],
        dask='allowed',
        vectorize=True,
        output_dtypes=[float])

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
        dask='parallelized',
        vectorize=True,
        output_dtypes=[float])


def msess_calculation(x, y, clim, dim):
    return xr.apply_ufunc(
        msess2, x, y, clim,
        input_core_dims=[[dim], [dim], [dim]],
        dask='allowed',
        vectorize=True,
        output_dtypes=[float])

#%% arguments
scen = sys.argv[1]
date = sys.argv[2]
date2 = sys.argv[3]

#%% folder
srel_recfolder = "output/ARM_run_" + date + "/" + scen + "/ARM_SrelD_" + date  + "/"

print("start")
# omit warnings (see above)
warnings.filterwarnings("ignore")

# get start time
stm = time.time()

# create validation output folder
indir = os.path.dirname(srel_recfolder)
outdir = indir + "/validation_" + time.strftime("%Y-%m-%d") + "/"
print(outdir)
if not os.path.isdir(outdir):
    try:
        os.mkdir(outdir)
    except OSError:
        print("")
    else:
        print("")


#%% 
# --------------------------------------------------------------------------------------------------------------------
# READ DATA
# --------------------------------------------------------------------------------------------------------------------
print("read observation data")
srel_obsfolder = "data/sunshine/" ## relative sunshine duration data from MeteoSwiss
srel_obsfiles = glob.glob(f'{srel_obsfolder}*SrelD*1971*.nc')
srel_obs = xr.open_mfdataset(srel_obsfiles).chunk({"time": -1, "E": "auto", "N":"auto"})
srel_obs["time"] = pd.to_datetime(srel_obs.time.dt.strftime("%Y-%m-%d").values)
srel_obs = srel_obs["SrelD"]

#%%
print("read reconstruction data")
srel_recfiles = glob.glob(f'{srel_recfolder}*srel_analogs_*.nc')
srel_rec = xr.open_mfdataset(srel_recfiles).chunk({"time": -1, "E": "auto", "N":"auto"})
srel_rec["time"] = pd.to_datetime(srel_rec.time.dt.strftime("%Y-%m-%d").values)
srel_rec
srel_rec = srel_rec.reindex(time = srel_obs["time"])["SrelD"]

seasons = {"MAM": np.arange(3,6),"SON":np.arange(9,12),"JJA":np.arange(6,9),"DJF": [12,1,2]}
print(seasons)

for name,value in seasons.items():
    print(value)
    print(name)

    obs_seas = srel_obs.sel(time=np.isin(srel_obs.time.dt.month,value))
    rec_seas = srel_rec.sel(time=np.isin(srel_rec.time.dt.month,value))
    clim = obs_seas.mean("time").expand_dims({'time': obs_seas['time']}).transpose(*obs_seas.dims)
      
    print("msess")
    msess_out = msess_calculation(obs_seas, rec_seas, clim, "time").compute()
    print("pcorr")
    pcorr_out = pearson_correlation(obs_seas, rec_seas, "time").compute()
    print("rmse")
    rmse_out = rmse_calculation(obs_seas, rec_seas, "time").compute()
    print("bias")
    bias_out = bias_calculation(obs_seas, rec_seas, "time").compute()
         
    metrics = xr.concat([pcorr_out, rmse_out,bias_out, msess_out], dim ="metrics", fill_value = np.nan)
    newcoords = ["corr_" + name, "rmse_"  + name,"bias_"  + name,"msess_" + name]
    metrics = metrics.assign_coords(metrics=newcoords)

    print("write netcdf")
    metrics.to_netcdf(outdir + time.strftime("%Y-%m-%d") + '_validation_maps_srel_' + name +'.nc')
    del metrics, newcoords

print("eval done")
    
    
