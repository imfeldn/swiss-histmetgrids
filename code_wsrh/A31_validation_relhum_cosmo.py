#!/usr/bin/env python
# coding: utf-8
# -----------------------------------------------
# validation of relhum
# -----------------------------------------------

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
msess = lambda x,y,clim: 1 - np.sum((x-y)**2, axis = -1)/np.sum((clim - x)**2, axis = -1)

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
        msess, x, y, clim,
        input_core_dims=[[dim], [dim], [dim]],
        dask='allowed',
        vectorize=True,
        output_dtypes=[float])
        
## source functions
#exec(open("A2_validation_functions.py").read())
# ### the script is called from the console as follows
# ### for line based run indicate filepath and -name in line 25 and uncomment line
#sys.argv = ["","scen3_80","2024-10-31","min"]

#%%
scen = sys.argv[1]
date = sys.argv[2]
rhtype = sys.argv[3]

#%%
relhum_recfolder = "/output/ARM_run_" + date + "/" + scen + "/ARM_"+ rhtype +"relhum_" + date + "/"

if rhtype == "min":
  rh = "daymin"
else:
  rh = "daily"

print("start")
print(scen)
# omit warnings (see above)
warnings.filterwarnings("ignore")

# get start time
stm = time.time()

# create validation output folder
indir = os.path.dirname(relhum_recfolder)
outdir = indir + "/validation_" + time.strftime("%Y-%m-%d") + "/"
print(outdir)
if not os.path.isdir(outdir):
    try:
        os.mkdir(outdir)
    except OSError:
        print("")
    else:
        print("")


#
# --------------------------------------------------------------------------------------------------------------------
# READ DATA
# --------------------------------------------------------------------------------------------------------------------
print("read observation data")
relhum_obsfolder = "/scratch3/nimfeld/wear/swiss_grids/data/humidity/cosmo/"
relhum_obsfiles = glob.glob(f'{relhum_obsfolder}*{rh}_RELHUM_2M_*LV95.nc')
relhum_obs = xr.open_mfdataset(relhum_obsfiles).chunk({"time": -1, "easting": "auto", "northing":"auto"})
relhum_obs["time"] = pd.to_datetime(relhum_obs.time.dt.strftime("%Y-%m-%d").values)
relhum_obs = relhum_obs["RH"]

print("read reconstruction data")
relhum_recfiles = glob.glob(f'{relhum_recfolder}*relhum_analogs_cosmo_*.nc')
relhum_rec = xr.open_mfdataset(relhum_recfiles).chunk({"time": -1, "easting": "auto", "northing":"auto"})
relhum_rec["time"] = pd.to_datetime(relhum_rec.time.dt.strftime("%Y-%m-%d").values)
relhum_rec

relhum_rec = relhum_rec.reindex(time = relhum_obs["time"])["RH"]

seasons = {"MAM": np.arange(3,6),"SON":np.arange(9,12),"JJA":np.arange(6,9),"DJF": [12,1,2]}
print(seasons)

for name,value in seasons.items():
    print(value)
    print(name)
    
    obs_seas = relhum_obs.sel(time=np.isin(relhum_obs.time.dt.month,value))
    rec_seas = relhum_rec.sel(time=np.isin(relhum_rec.time.dt.month,value))
    clim = obs_seas.mean("time").expand_dims({'time': obs_seas['time']}).transpose(*obs_seas.dims)
   
    print("pcorr")
    pcorr_out = pearson_correlation(obs_seas, rec_seas, "time").compute()
    print("rmse")
    rmse_out = rmse_calculation(obs_seas, rec_seas, "time").compute()
    print("bias")
    bias_out = bias_calculation(obs_seas, rec_seas, "time").compute()
    print("msess")
    msess_out = msess_calculation(obs_seas, rec_seas, clim, "time").compute()

    metrics = xr.concat([pcorr_out, rmse_out,bias_out, msess_out], dim ="metrics", fill_value = np.nan)
    newcoords = ["corr_" + name, "rmse_"  + name,"bias_"  + name,"msess_"  + name]
    metrics = metrics.assign_coords(metrics=newcoords)

    print("write netcdf")
    metrics.to_netcdf(outdir + time.strftime("%Y-%m-%d") + '_validation_maps_' + rhtype + 'relhum_' + name +'.nc')
    del metrics, newcoords

print("eval done")
    
    
