# create yearly netCDF from best analogue dates

# import packages
import numpy as np
import pandas as pd
import xarray as xr
import time
import sys
import glob
import os

print(sys.prefix)
print(sys.argv)

# get start time
stm = time.time()

#%%
# check if input arguments exist
if len(sys.argv) < 2:
  sys.exit("indicate folder path to analogue dates")

year= int(sys.argv[2])
start = str(sys.argv[3])
dist = str(sys.argv[4])
scen = str(sys.argv[5])
date = str(sys.argv[6])

os.chdir("output/ARM_run_" + date)

#
# read folder
dates_ifolder = str(sys.argv[1])  # get first argument passed to script (folder of analogue dates)
print(dates_ifolder)

# create folder for fields
outdir = "output/ARM_run_" + date + "/" + scen   + "/ARM_SrelD_" + date + "/"
print(outdir)
if not os.path.isdir(outdir):
     os.makedirs(outdir)

#%%--------------------------------------------------------------------------------------------------------------------
# READ DATA
# --------------------------------------------------------------------------------------------------------------------
print("read readme")
print(os.getcwd())
print(glob.glob(f"{dates_ifolder}*readme*{dist}*{scen}.txt"))
print(glob.glob(f'{dates_ifolder}*analog.dates*{dist}*{scen}.txt'))

analog_readme = np.array(pd.read_csv(glob.glob(f"{dates_ifolder}*readme*{dist}*{scen}.txt")[0], sep=" "))  # read csv with "readme" in filename
endCal = pd.to_datetime(analog_readme[np.argwhere(analog_readme[:,0]=='endCal:').item(0),1])
analog_dates = pd.read_csv(glob.glob(f"{dates_ifolder}*analog.dates*{dist}*{scen}.txt")[0], sep=" ")  # read csv with "analog.dates" in filename

# get target dates and best analogue dates
dates = pd.to_datetime(analog_dates["date"].values)
analog_1_dates = pd.to_datetime(analog_dates["1"], unit='D', origin=pd.Timestamp("1970-01-01"))  # get vector of best analogue dates
analog_1_dates = analog_1_dates.dt.strftime("%Y-%m-%d")

#print(endCal)
#print(dates)

selyear = dates.year == year
print(year)
analog_1_dates.index = dates


#%%###########################################
# read raster data
###########################################
obs_ifolder = "data/" # this loads the relativate sunshine duration data from meteoswiss (https://www.meteoswiss.admin.ch/dam/jcr:981891db-30d1-47cc-a2e1-50c270bdaf22/ProdDoc_SrelD.pdf)
srel_files = glob.glob(f'{obs_ifolder}*SrelD*1971*')[0]

#
print("read data")
srel_obs = xr.open_dataset(srel_files,  engine = 'netcdf4')
srel_obs["time"] = srel_obs.time.dt.strftime("%Y-%m-%d")
srel_analogs = srel_obs.reindex(time=analog_1_dates[selyear])    
srel_analogs["time"] = pd.to_datetime(analog_dates["date"][selyear]).values

#%% write raster data to netCDF files
print("write data to netCDF files")
dist = np.argwhere(analog_readme[:,0]=='distance:').item(0)
scal = np.argwhere(analog_readme[:,0]=='startCal:').item(0)
scen = np.argwhere(analog_readme[:,0]=='scenario:').item(0)
#

print(scen)
srel_analogs.to_netcdf(outdir + analog_readme[0, 1] + "_srel_analogs_" + str(sys.argv[2]) + "-01-01-" + str(sys.argv[2])  + "-12-31_" + analog_readme[dist, 1] + "_" + analog_readme[scen, 1]  + ".nc", format='NETCDF4')
#
#%%
del srel_analogs,srel_obs
etm = time.time()
print("srel fields created in: " + str(round((etm-stm)//60)) + " minutes, " + str(round((etm-stm) % 60, 2)) + " seconds")
print("files are in: " + outdir)
#
#### END
