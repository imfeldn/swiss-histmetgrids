#!/usr/bin/env python
# coding: utf-8
# -----------------------------------------------
# validation of wind speed
# -----------------------------------------------

# import packages 
import numpy as np
import pandas as pd
import xarray as xr

import sys
import glob
import time
import warnings 
import os

np.seterr(invalid='ignore')

## source functions
exec(open("A2_validation_functions.py").read())

#%% recalculate winddir/windspeed, or load it from directory
calc = False

#%% calculate for ARM only or for EnKF
arm = False

if arm:
    wind_recfolder = "/output/ARM_run_" + date + "/" + scen + "/"
    arm_str = "_ARM"
else:
    wind_recfolder = "/output/ARM_run_" + date + "/" + scen + "/EnKF_windspeed_" + EnKF + "_" + date2 + "/"
    arm_str = "_EnKF"
    
print(wind_recfolder)
print("start")
# omit warnings (see above)
warnings.filterwarnings("ignore")

# get start time
stm = time.time()

# check if file is valid
if not os.path.isdir(str(wind_recfolder)):
    sys.exit("file path invalid")

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

# orig data
print("read observation data")
wind_obsfile = "/data/wind/cosmo/cosmo_20162020_daymean_windspeed_LV95.nc"
wspeed_obs = xr.open_dataset(wind_obsfile, engine ='netcdf4').chunk({"time": -1, "easting" : "auto", "northing" : "auto"})
wspeed_obs["time"] = pd.to_datetime(wspeed_obs.time.dt.strftime("%Y-%m-%d").values)
wspeed_obs = wspeed_obs["ws"]

wind_obsfile = "data/wind/cosmo/cosmo_20162020_daymean_winddirection_LV95.nc"
wdir_obs = xr.open_dataset(wind_obsfile, engine ='netcdf4').chunk({"time": -1, "easting" : "auto", "northing" : "auto"})
wdir_obs["time"] = pd.to_datetime(wdir_obs.time.dt.strftime("%Y-%m-%d").values)
wdir_obs = wdir_obs["ws"]

#%% load recon data
if calc:
    if arm:
        uwind_files = glob.glob(f'{wind_recfolder}*Uwind_*{dist}*')
        vwind_files = glob.glob(f'{wind_recfolder}*Vwind_*{dist}*')
     
        uwind_rec = xr.open_mfdataset(uwind_files).chunk({"time": -1})
        vwind_rec = xr.open_mfdataset(vwind_files).chunk({"time": -1})
        wspeed_rec = (uwind_rec.ws**2 + vwind_rec.ws**2)**(1/2)
        wspeed_rec.to_netcdf(wind_recfolder +  date2 +'_windspeed_'+ dist + arm_str + '_2016-01-01-2020-12-31.nc')
     
        wdir_rec = np.degrees(np.arctan2(uwind_rec, vwind_rec)) % 360
        wdir_rec = wdir_rec["ws"]
        wdir_rec.to_netcdf(wind_recfolder +  date2 +'_wdir_' + dist + arm_str + '_2016-01-01-2020-12-31.nc')
    else:
        uwind_files = glob.glob(f'{wind_recfolder}*Uwind_*')
        vwind_files = glob.glob(f'{wind_recfolder}*Vwind_*')
     
        uwind_rec = xr.open_mfdataset(uwind_files).chunk({"time": -1})
        vwind_rec = xr.open_mfdataset(vwind_files).chunk({"time": -1})
        wspeed_rec = (uwind_rec.ws**2 + vwind_rec.ws**2)**(1/2)
        wspeed_rec.to_netcdf(wind_recfolder +  date2 +'_windspeed'+  arm_str + '_2016-01-01-2020-12-31.nc')
     
        wdir_rec = np.degrees(np.arctan2(uwind_rec, vwind_rec)) % 360
        wdir_rec = wdir_rec["ws"]
        wdir_rec.to_netcdf(wind_recfolder +  date2 +'_wdir'  + arm_str + '_2016-01-01-2020-12-31.nc')

else:
    if arm:
        print("read arm reconstruction data")
        wspeed_recfiles = glob.glob(f'{wind_recfolder}*windspeed_*{dist}*ARM')
        wspeed_rec = xr.open_mfdataset(wspeed_recfiles).chunk({"time": -1})
    else:
        print("read enkf reconstruction data")
        wspeed_recfiles = glob.glob(f'{wind_recfolder}*_windspeed*EnKF_2016-01-01-2020-12-31.nc')[0]
        wspeed_rec = xr.open_dataset(wspeed_recfiles).chunk({"time": -1})["ws"]

wspeed_rec = wspeed_rec.reindex(northing=wspeed_obs.northing)
wspeed_rec = wspeed_rec.reindex(time = wspeed_obs["time"])

#%% wind speed
seasons = {"JJA":np.arange(6,9),"DJF": [12,1,2],"SON":np.arange(9,12),"MAM": np.arange(3,6)}

#%%
for name,value in seasons.items():
    print(value)
    print(name)

    obs_seas = wspeed_obs.sel(time=np.isin(wspeed_obs.time.dt.month,value))
    rec_seas = wspeed_rec.sel(time=np.isin(wspeed_rec.time.dt.month,value))
    print("pcorr")
    pcorr_out = pearson_correlation(obs_seas, rec_seas, "time").compute()
    print("rmse")
    rmse_out = rmse_calculation(obs_seas, rec_seas, "time").compute()
    print("rmse")
    norm_rmse_out = norm_rmse_calculation(obs_seas, rec_seas, "time").compute()
    print("bias")
    bias_out = bias_calculation(obs_seas, rec_seas, "time").compute()
    print("msess")
    msess_out = msess_calculation(obs_seas, rec_seas, "time").compute()

    metrics = xr.concat([pcorr_out, rmse_out, norm_rmse_out,bias_out, msess_out], dim ="metrics", fill_value = np.nan)
    newcoords = ["corr_" + name, "rmse_"  + name,"norm_rmse_" + name,"bias_"  + name,"msess_"  + name]
    metrics = metrics.assign_coords(metrics=newcoords)

    print("write netcdf")
    metrics.to_netcdf(outdir + time.strftime("%Y-%m-%d") + '_validation_maps_wspeed_' + name +'.nc')
    del metrics, newcoords

print("eval done")


#%%  autocorrelation
for name,value in seasons.items():
    print(value)
    print(name)

    obs_seas = wspeed_obs.sel(time=np.isin(wspeed_obs.time.dt.month,value))
    rec_seas = wspeed_rec.sel(time=np.isin(wspeed_rec.time.dt.month,value))
    print("pcorr")
    acorr_obs = auto_correlation(obs_seas, "time").compute()
    acorr_rec = auto_correlation(rec_seas, "time").compute()
    acorr_diff = acorr_rec - acorr_obs
    
    allvals = xr.concat([acorr_obs, acorr_rec, acorr_diff], dim ="autocorr", fill_value = np.nan)
    newcoords = ["obs_" + name, "rec_"  + name, "diff_"  + name]
    allvals = allvals.assign_coords(autocorr = newcoords)
    
    print("write netcdf")
    allvals.to_netcdf(outdir + time.strftime("%Y-%m-%d") + '_autocorr_maps_wspeed_' + name +'.nc')
   
print("eval done")

#%%  wind direction

for name,value in seasons.items():
    print(value)
    print(name)

    obs_seas = wdir_obs.sel(time=np.isin(wdir_obs.time.dt.month,value)).load()
    rec_seas = wdir_rec.sel(time=np.isin(wdir_rec.time.dt.month,value)).load()
    print("mae")
    mae_out = mae_calculation_wd(obs_seas, rec_seas, "time")
    print("bias")
    bias_out = bias_calculation_wd(obs_seas, rec_seas, "time")

    metrics = xr.concat([mae_out,bias_out], dim ="metrics", fill_value = np.nan)
    newcoords = ["mae_"  + name, "bias_"  + name]
    metrics = metrics.assign_coords(metrics=newcoords)

    print("write netcdf")
    metrics.to_netcdf(outdir + time.strftime("%Y-%m-%d") + '_validation_maps_wdir_' + name +'.nc')
    del metrics, newcoords

print("eval done")

