# create yearly netCDF from best analogue dates

# import packages
import numpy as np
import pandas as pd
import xarray as xr
import time
import sys
import glob
import os

#%% uncomment for line based run
# sys.argv = ["", "/path/to/file/", 2016,"2016-01-01","Gower","scen3_80","2024-10-31"]

print(sys.prefix)
print(sys.argv)

# get start time
stm = time.time()

#%%
# check if input arguments exist
if len(sys.argv) < 2:
    sys.exit("indicate folder path to analogue dates")

start = str(sys.argv[3])
dist = str(sys.argv[4])
scen = str(sys.argv[5])
date = str(sys.argv[6])

os.chdir("ARM_run_" + date)

# read folder
dates_ifolder = str(sys.argv[1])  # get first argument passed to script (folder of analogue dates)
print(dates_ifolder)

# create folder for fields
outdir = "ARM_run_" + date + "/" + scen + "/ARM_relhum_" + date + "/"
print(outdir)
if not os.path.isdir(outdir):
    try:
         os.makedirs(outdir)
    except OSError:
        print("")
    else:
        print("")


#%%--------------------------------------------------------------------------------------------------------------------
# READ DATA
# --------------------------------------------------------------------------------------------------------------------
print("read readme")

analog_readme = np.array(pd.read_csv(glob.glob(f"{dates_ifolder}*readme*{dist}*{scen}.txt")[0], sep=" "))  # read csv with "readme" in filename
endCal = pd.to_datetime(analog_readme[np.argwhere(analog_readme[:,0]=='endCal:').item(0),1])
analog_dates = pd.read_csv(glob.glob(f"{dates_ifolder}*analog.dates*{dist}*{scen}.txt")[0], sep=" ")  # read csv with "analog.dates" in filename

# get target dates and best analogue dates
dates = pd.to_datetime(analog_dates["date"].values)
analog_1_dates = pd.to_datetime(analog_dates["1"], unit='D', origin=pd.Timestamp("1970-01-01"))  # get vector of best analogue dates
analog_1_dates = analog_1_dates.dt.strftime("%Y-%m-%d")

    
year= int(sys.argv[2])
selyear = dates.year == year
print(year)

#%%###########################################
# read raster data
###########################################
obs_ifolder = "/data/humidity/cosmo/"
w_files = glob.glob(f'{obs_ifolder}*daily_RELHUM_2M*LV95.nc')

print("read data")
relhum_obs = xr.open_mfdataset(w_files,  engine = 'netcdf4', chunks = dict(time = 10) )
relhum_obs["time"] = relhum_obs.time.dt.strftime("%Y-%m-%d")
#%%
relhum_analogs = relhum_obs.reindex(time=analog_1_dates[selyear])    
relhum_analogs["time"] = pd.to_datetime(analog_dates["date"][selyear]).values


#%% write raster data to netCDF files
print("write data to netCDF files")
dist = np.argwhere(analog_readme[:,0]=='distance:').item(0)
scal = np.argwhere(analog_readme[:,0]=='startCal:').item(0)
scen = np.argwhere(analog_readme[:,0]=='scenario:').item(0)

print(scen)
filename = analog_readme[0, 1] + "_relhum_analogs_cosmo_" + str(sys.argv[2]) + "-01-01-" + str(sys.argv[2])  + "-12-31_" + analog_readme[dist, 1] + "_" + analog_readme[scen, 1]  + ".nc"
relhum_analogs.to_netcdf(outdir +  "help_" + str(sys.argv[2]) + ".nc", format='NETCDF4')
#%%%

print(outdir + "help_" + str(sys.argv[2]) + ".nc")
os.system('cdo -f nc4 -z zip copy ' + outdir + "help_" + str(sys.argv[2]) + ".nc" + ' ' + outdir + filename)
os.remove(outdir + "help_" + str(sys.argv[2]) + ".nc")


#%%
del relhum_analogs, relhum_obs
etm = time.time()
print("relhum fields created in: " + str(round((etm-stm)//60)) + " minutes, " + str(round((etm-stm) % 60, 2)) + " seconds")
print("files are in: " + outdir)
#
#### END
