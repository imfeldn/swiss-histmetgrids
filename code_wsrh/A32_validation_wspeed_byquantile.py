#!/usr/bin/env python
# coding: utf-8
# -----------------------------------------------
# validation of wind speed
# -----------------------------------------------

# import packages 
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from natsort import natsorted
import matplotlib.patches as mpatches

import glob
import time
import warnings 
import os

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


#%% In[4]:
scen = "scen5_testMIL23"
window = "80"
arm = False
date = "2024-10-31"
date2 = "2024-11-20"
EnKF = "ta_p_sp_spm_anom_blend0.5_Rhist2"
add = ""

if arm:
    wind_recfolder = "/output/ARM_run_" + date + "/" + scen + "_" + window  + add + "/"
else:
    wind_recfolder = "/output/ARM_run_" + date + "/" + scen + "_" + window + add + "/EnKF_windspeed_" + EnKF + "_" + date2 + "/"


print("start")
# omit warnings (see above)
warnings.filterwarnings("ignore")

# create validation output folder
outdir = wind_recfolder + "/validation_" + time.strftime("%Y-%m-%d") + "/"

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
wind_obsfile = "/data/wind/cosmo/cosmo_20162020_daymean_windspeed_LV95.nc"
wspeed_obs = xr.open_dataset(wind_obsfile, engine ='netcdf4').chunk({"time": -1, "easting" : "auto", "northing" : "auto"})
wspeed_obs["time"] = pd.to_datetime(wspeed_obs.time.dt.strftime("%Y-%m-%d").values)
wspeed_obs = wspeed_obs["ws"]

# recon data
calc = False
if calc:
    if arm:
        print("read arm reconstruction data")
        wind_recfiles = glob.glob(f'{wind_recfolder}*wind_analogs_cosmo_*')
        wspeed_rec = xr.open_mfdataset(wind_recfiles).chunk({"time": -1})
        wspeed_rec = wspeed_rec['ws']

    else:
        print("read enkf reconstruction data")
        uwind_files = glob.glob(f'{wind_recfolder}*Uwind_EnKF*')
        vwind_files = glob.glob(f'{wind_recfolder}*Vwind_EnKF*')

        uwind_rec = xr.open_mfdataset(uwind_files).chunk({"time": -1})
        vwind_rec = xr.open_mfdataset(vwind_files).chunk({"time": -1})
        wspeed_rec = (uwind_rec.ws**2 + vwind_rec.ws**2)**(1/2)
        wspeed_rec.to_netcdf(wind_recfolder +  date2 +'_windspeed_EnKF_2016-01-01-2020-12-31.nc')
        
else:    
    if arm:
        print("read arm reconstruction data")
        wspeed_recfiles = glob.glob(f'{wind_recfolder}*windspeed_ARM_*')
        wspeed_rec = xr.open_mfdataset(wspeed_recfiles).chunk({"time": -1})
    else:
        print("read enkf reconstruction data")
        wspeed_recfiles = glob.glob(f'{wind_recfolder}*_windspeed_EnKF_2016-01-01-2020-12-31.nc')[0]
        wspeed_rec = xr.open_dataset(wspeed_recfiles).chunk({"time": -1})["ws"]


wspeed_rec = wspeed_rec.reindex(northing=wspeed_obs.northing)
wspeed_rec = wspeed_rec.reindex(time = wspeed_obs["time"])

#%% define seasons
seasons = {"JJA":np.arange(6,9),"DJF": [12,1,2],"MAM": np.arange(3,6),"SON":np.arange(9,12)}

#%% run validation
print(seasons)

for name,value in seasons.items():
    print(value)
    print(name)
    
    obs_seas = wspeed_obs.sel(time=np.isin(wspeed_obs.time.dt.month,value))
    rec_seas = wspeed_rec.sel(time=np.isin(wspeed_rec.time.dt.month,value))
    
    qt_dims = ("time")
    qt_values = np.arange(0.1,1.1,0.1)

    ds_qt = obs_seas.quantile(qt_values, dim=qt_dims)
    ## reorder to vecotr and remove NA values 
    rmse_out = []
    bias_out = []

    for qt in np.arange(0,len(qt_values)):
        print(qt)
        if qt == 0:
            selqs_obs = obs_seas.where(obs_seas < ds_qt[qt,:,:]).stack(points = ("easting","northing"))
            selqs_rec = rec_seas.where(obs_seas < ds_qt[qt,:,:]).stack(points = ("easting","northing"))    
        else:
            selqs_obs = obs_seas.where((obs_seas < ds_qt[qt,:,:]) & (obs_seas > ds_qt[qt - 1,:,:])).stack(points = ("easting","northing"))
            selqs_rec = rec_seas.where((obs_seas < ds_qt[qt,:,:]) & (obs_seas > ds_qt[qt - 1,:,:])).stack(points = ("easting","northing"))  
        print("rmse")
        rmse_out.append(rmse_calculation(selqs_obs, selqs_rec, "points").compute().to_numpy())   
        print("bias")
        bias_out.append(bias_calculation(selqs_obs, selqs_rec, "points").compute().to_numpy())
       
    bias_mat = np.array(bias_out)
    rmse_mat = np.array(rmse_out)
    
    np.savetxt(outdir + "bias_byquantile_" + name + ".csv",  bias_mat,  delimiter = ",")
    np.savetxt(outdir + "rmse_byquantile_" + name + ".csv",  rmse_mat,  delimiter = ",")
