
# import packages 
import numpy as np
import xarray as xr
from natsort import natsorted

import glob
import time
import os

np.seterr(invalid='ignore')

tmax = True
rhmin = True

folder = str(np.where(rhmin,"rhmin",""))
folder = folder + str(np.where(tmax,"_tmax",""))
outdir =  "output/fwi/recon/" + folder + "/"

# execute funtions
exec(open("F0_fire_indices_functions.py").read())

#%%  list all files
# thses data inputs need to be downloaded from PANGAEA and BORIS, see publications for references
if tmax:
    temp_folder = "boris/tmax/"
    temp_files = glob.glob(f'{temp_folder}*tmax*')
else:    
    temp_folder = "/recon/temp/"
    temp_files = glob.glob(f'{temp_folder}*temp*')

prec_folder = "/recon/precip/"
prec_files = glob.glob(f'{prec_folder}*precip*')


wind_folder = "/boris/wspeed/"
wind_files =  glob.glob(f'{wind_folder}*windspeed*')
    

if rhmin:
    rh_folder = "/boris/rh/daymin/"
    rh_files = glob.glob(f'{rh_folder}*daymin_relhum*')
else:
    rh_folder = "/boris/rh/daymean/"
    rh_files = glob.glob(f'{rh_folder}*relhum_analogs*')

      

#%% load a file with latlon values
latlondat = xr.open_dataset("CH_temp_TabsD_2020.nc")

#%% set years for calculation | in total 250 years
yrlist = np.arange(2020,2021)
yrlist = np.arange(2019,2021)
#yrlist = [1763]

start = yrlist[0]
lat2 = "N"
#%%
for yr in yrlist:
    
    print(yr)
    tic = time.perf_counter()
    if yr < 1961:
        tvar = "tmax"
        rvar = "precip"
    else:
        tvar = "tmax"
        rvar= "RhiresD"
    
    # get yearly file/s
    tf = list(filter(lambda k: str(yr) in k, temp_files))[0]
    pf = list(filter(lambda k: str(yr) in k, prec_files))[0]
    rhf = list(filter(lambda k: str(yr) in k, rh_files))[0]
    wf = list(filter(lambda k: str(yr) in k, wind_files))[0]

    tdat = xr.open_mfdataset(tf).load()
    pdat = xr.open_mfdataset(pf).load()
    rhdat = xr.open_mfdataset(rhf).load()
    wdat = xr.open_mfdataset(wf).load()
    wdat = wdat.reindex(northing=wdat.northing[::-1])

    out_dmc = xr.full_like(tdat, np.nan, dtype=np.float32)
    out_dc = xr.full_like(tdat, np.nan, dtype=np.float32)
    out_bui = xr.full_like(tdat, np.nan, dtype=np.float32)
    out_ffmc = xr.full_like(tdat, np.nan, dtype=np.float32)
    out_isi = xr.full_like(tdat, np.nan, dtype=np.float32)
    out_fwi = xr.full_like(tdat, np.nan, dtype=np.float32)

    out_dmc = out_dmc.rename({tvar : "dmc"})
    out_dc = out_dc.rename({tvar : "dc"})
    out_bui = out_bui.rename({tvar : "bui"})
    out_ffmc = out_ffmc.rename({tvar : "ffmc"})
    out_isi = out_isi.rename({tvar : "isi"})
    out_fwi = out_fwi.rename({tvar : "fwi"})
    
    tsteps = np.arange(0,out_dmc.dims["time"])
    sellatlon = pdat["N"].isin(tdat[lat2])
#
    for ts in tsteps:
        print(ts)
        if (ts < 2 and yr == 1763):
            dmc0 = xr.full_like(out_dmc["dmc"][1,:,:], 6.0)
            dc0 = xr.full_like(out_dc["dc"][1,:,:], 15.0)
            ffmc0 = xr.full_like(out_ffmc["ffmc"][1,:,:], 85.0)  
        elif (ts == 0 and yr == start):
            dmc0 = xr.full_like(out_dmc["dmc"][1,:,:], 6.0)
            dc0 = xr.full_like(out_dc["dc"][1,:,:], 15.0)
            ffmc0 = xr.full_like(out_ffmc["ffmc"][1,:,:], 85.0)   
        elif (ts == 0 and yr > start):
            # if ds:
            #     ds = xr.open_dataset()    
            lastday = len(ds["time"]) -1 # get index of last day of previous year
            dmc0 = ds["dmc"][lastday,:,:]
            dc0 = ds["dc"][lastday,:,:]
            ffmc0 = ds["ffmc"][lastday,:,:]
            print("took last years value")
            del ds
            
        else:
            dmc0 = out_dmc["dmc"][ts-1,:,:]
            dc0 = out_dc["dc"][ts-1,:,:]
            ffmc0 = out_ffmc["ffmc"][ts-1,:,:]
        
        out_dmc["dmc"][ts,:,:] = dmc(tdat[tvar][ts,:,:].values, 
                                    pdat[rvar][ts,sellatlon,:].values, 
                                    rhdat["rh"][ts,:,:].values, 
                                    tdat["time"][ts].dt.month.values,
                                    latlondat["lat"].values,  
                                    dmc0.values)
        
        out_dc["dc"][ts,:,:] = dc(tdat[tvar][ts,:,:].values, 
                                  pdat[rvar][ts,sellatlon,:].values, 
                                  tdat["time"][ts].dt.month.values,
                                  latlondat["lat"].values,  
                                  dc0.values)
        
        out_bui["bui"][ts,:,:] = bui(out_dmc["dmc"][ts,:,:].values, 
                                    out_dc["dc"][ts,:,:].values)
        # set bui to 0 where DMC is 0, see: https://confluence.ecmwf.int/download/attachments/239344103/Fire_In_CDS.pdf?version=1&modificationDate=1634645741421&api=v2    
        out_bui['bui'][ts,:,:] = out_bui['bui'][ts,:,:].where((out_dmc["dmc"][ts,:,:].isnull()) | (out_dmc["dmc"][ts,:,:]>0), 0)

        out_ffmc["ffmc"][ts,:,:] = ffmc(tdat[tvar][ts,:,:].values, 
                                  pdat[rvar][ts,sellatlon,:].values,
                                  wdat["ws"][ts,:,:].values, 
                                  rhdat["rh"][ts,:,:].values, 
                                  ffmc0.values)
    
        out_isi["isi"][ts,:,:] = isi(wdat["ws"][ts,:,:].values, 
                                  out_ffmc["ffmc"][ts,:,:].values)  
    
        out_fwi["fwi"][ts,:,:] = fwi(out_isi["isi"][ts,:,:].values, 
                                  out_bui["bui"][ts,:,:].values)   
                               
    toc = time.perf_counter()
    print(f"one year took {toc - tic:0.4f} seconds")
                       
    ds = xr.merge([out_dmc, out_dc, out_bui, out_ffmc, out_isi, out_fwi])
    
    print("write netcdf")
    if yr == 2020:
        ds.to_netcdf(outdir +  'help_fwi_' + str(yr) +'.nc')
    
        filename = outdir + time.strftime("%Y-%m-%d") + '_fwi_indices_' + str(yr) +'.nc'
        os.system('cdo -f nc4 -z zip  copy ' + outdir +  'help_fwi_' + str(yr) +'.nc' + ' ' + filename)
        os.system('rm ' + outdir +  'help_fwi_' + str(yr) +'.nc')

    
    del out_dmc, out_dc, out_bui, out_ffmc, out_isi, out_fwi
    del rhdat, tdat, wdat, pdat
                    


