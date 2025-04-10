
# import packages 
import numpy as np
import xarray as xr

import glob
import time
import os

np.seterr(invalid='ignore')
outdir =  "output/fwi/cosmo/"

# execute funtions
exec(open("F0_fire_indices_functions.py").read())

#%%  list all files
cosmo_folder = "/data/cosmo_fwi/"
cosmo_files = glob.glob(f'{cosmo_folder}*fwi_vars_*12CET.nc')
prec_files = glob.glob(f'{cosmo_folder}*cosmo_TOT_PREC*')

prec_dat = xr.open_mfdataset(prec_files)["TOT_PREC"]
# select every third however and calculate daily sum based on this
prec_3h = prec_dat[np.arange(0,36584,3),:,:]
precip_12cet = prec_3h.resample(time='1D', loffset='11H', closed='right', label='right').sum()

#%% set years for calculation 
yrlist = np.arange(2016,2021)
start = yrlist[0]

#%%
for yr in yrlist:
    
    print(yr)
    tic = time.perf_counter()
    tvar = "T_2M"
    rvar = "TOT_PREC"
    rhvar = "RELHUM_2M"
    
    # get yearly file/s
    yrfile = list(filter(lambda k: str(yr) in k, cosmo_files))
    prfile = list(filter(lambda k: str(yr) in k, prec_files))
    
    tdat = xr.open_mfdataset(yrfile, drop_variables = ["TOT_PREC","U_10M","V_10M","RELHUM_2M"]).load()
    rhdat = xr.open_mfdataset(yrfile, drop_variables = ["T_2M", "U_10M","V_10M","TOT_PREC"]).load()[rhvar]
    uvdat = xr.open_mfdataset(yrfile, drop_variables = ["TOT_PREC","RELHUM_2M"]).load()
    pdat = precip_12cet.where(precip_12cet.time.dt.year == yr, drop=True)
   
    ## calc windspeed
    wdat = (uvdat.V_10M**2 + uvdat.U_10M**2)**(1/2)
    del uvdat
    
    latlonvals = tdat["lat_1"]
   
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
    
    # make array and change to Celsius
    tdat = tdat[tvar] - 273.15
    tsteps = np.arange(0,out_dmc.dims["time"])

    for ts in tsteps:
        print(ts)
        if (ts == 0 and yr == start):
            dmc0 = xr.full_like(out_dmc["dmc"][1,:,:], 6.0)
            dc0 = xr.full_like(out_dc["dc"][1,:,:], 15.0)
            ffmc0 = xr.full_like(out_ffmc["ffmc"][1,:,:], 85.0)
        elif (ts == 0 and yr > start):
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
        
        out_dmc["dmc"][ts,:,:] = dmc(tdat[ts,:,:].values, 
                                    pdat[ts,:,:].values, 
                                    rhdat[ts,:,:].values, 
                                    tdat["time"][ts].dt.month.values,
                                    tdat["lat_1"].values,  
                                    dmc0.values)
        
        out_dc["dc"][ts,:,:] = dc(tdat[ts,:,:].values, 
                                  pdat[ts,:,:].values, 
                                  tdat["time"][ts].dt.month.values,
                                  tdat["lat_1"].values,  
                                  dc0.values)
        
        out_bui["bui"][ts,:,:] = bui(out_dmc["dmc"][ts,:,:].values, 
                                    out_dc["dc"][ts,:,:].values)
        # set bui to 0 where DMC is 0, see: https://confluence.ecmwf.int/download/attachments/239344103/Fire_In_CDS.pdf?version=1&modificationDate=1634645741421&api=v2    
        out_bui['bui'][ts,:,:] = out_bui['bui'][ts,:,:].where((out_dmc["dmc"][ts,:,:].isnull()) | (out_dmc["dmc"][ts,:,:]>0), 0)

        out_ffmc["ffmc"][ts,:,:] = ffmc(tdat[ts,:,:].values, 
                                  pdat[ts,:,:].values,
                                  wdat[ts,:,:].values, 
                                  rhdat[ts,:,:].values, 
                                  ffmc0.values)
    
        out_isi["isi"][ts,:,:] = isi(wdat[ts,:,:].values, 
                                  out_ffmc["ffmc"][ts,:,:].values)  
    
        out_fwi["fwi"][ts,:,:] = fwi(out_isi["isi"][ts,:,:].values, 
                                  out_bui["bui"][ts,:,:].values)   
                               
    toc = time.perf_counter()
    print(f"one year took {toc - tic:0.4f} seconds")
                       
    ds = xr.merge([out_dmc, out_dc, out_bui, out_ffmc, out_isi, out_fwi])
    
    print("write netcdf")
    ds.to_netcdf(outdir +  'help_fwi_' + str(yr) +'.nc')
    
    filename = outdir + time.strftime("%Y-%m-%d") + '_cosmo_fwi_indices_12CET' + str(yr) +'.nc'
    os.system('cdo -f nc4 -z zip  copy ' + outdir +  'help_fwi_' + str(yr) +'.nc' + ' ' + filename)
    os.system('rm ' + outdir +  'help_fwi_' + str(yr) +'.nc')

    
    del out_dmc, out_dc, out_bui, out_ffmc, out_isi, out_fwi
    del rhdat, tdat, wdat, pdat
                    


