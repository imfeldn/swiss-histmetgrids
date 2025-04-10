###  CALCULATE Covariance matrix based on weather types ###
rm(list=ls())

library(ncdf4)
library(lubridate)
library(terra)

# files
ufile <- "data/cosmo_20162020_daymean_U_10M_LV95.nc"
vfile <- "data/cosmo_20162020_daymean_V_10M_LV95.nc"
psfile <- "data/pressure_anom_station_cosmo.rds"
tafile <-  "data/temp_deseas_station_cosmo.rds"
gradfile <-  "data/pressure_station_cosmo.rds"
TOTfile <- "data/analogue_stationdata_sp_presshom.RData"

### read inventory
station_coords <- XLConnect::readWorksheetFromFile("data/station_inventory.xlsx",sheet=1, startRow=2)[,c("ID","E_lv95_ref","N_lv95_ref")]
station_coords$E_lv95_ref <- as.numeric(station_coords$E_lv95_ref)
station_coords$N_lv95_ref <- as.numeric(station_coords$N_lv95_ref)

### read station file
load(TOTfile)

### load wind 
vdat <- terra::rast(vfile)
v_time <- as.Date(terra::time(vdat))
udat <- terra::rast(ufile)
u_time <- as.Date(terra::time(udat))

### load temperature gridded data
temp_deseas_gc <- readRDS(file = tafile)
colnames(temp_deseas_gc) <- as.character(as.Date(as.numeric(colnames(temp_deseas_gc))))

### load pressure gridded data
pressure_anom_gc <- readRDS(file = psfile)
colnames(pressure_anom_gc) <- as.character(as.Date(as.numeric(colnames(pressure_anom_gc))))

### get spatial pressure gradients with respect to Lugano
pressure_gc <- readRDS(file = gradfile)
spatdiff_gc <- pressure_gc
spatdiff_gc[which(!grepl("LUG_p",rownames(spatdiff_gc))),]  <- spatdiff_gc[which(!grepl("LUG_p",rownames(spatdiff_gc))),] - pressure_gc["LUG_p",]
spatdiff_gc["LUG_p",]  <- pressure_gc["LUG_p",] - pressure_gc["ALT_p",] 
rownames(spatdiff_gc) <- gsub("_p","_sp",rownames(pressure_gc))
spatdiff_anom_gc <- t(scale(t(spatdiff_gc), scale = F))
saveRDS(spatdiff_anom_gc,file = "output/spatdiff_sp_anom_gc.rds",compress = "xz")

### get spatial pressure gradients with respect to Milano
spatdiff_mil_gc <- pressure_gc
spatdiff_mil_gc[which(!grepl("MIL_p",rownames(spatdiff_mil_gc))),]  <- spatdiff_mil_gc[which(!grepl("MIL_p",rownames(spatdiff_mil_gc))),] - pressure_gc["MIL_p",]
spatdiff_mil_gc["MIL_p",]  <- pressure_gc["MIL_p",] - pressure_gc["SMA_p",] 
rownames(spatdiff_mil_gc) <- gsub("_p","_spm",rownames(spatdiff_mil_gc))
spatdiff_mil_anom_gc <- t(scale(t(spatdiff_mil_gc), scale = F))
saveRDS(spatdiff_mil_anom_gc,file = "output/spatdiff_spm_anom_gc.rds",compress = "xz")

### create separate variables for ew,ns gradients
pgrads <- t(data.frame(NSG_sp = pressure_gc["ALT_p",] - pressure_gc["LUG_p",],
                       EWG_sp = pressure_gc["GVE_p",] - pressure_gc["SMA_p",],
                       ESG_sp = pressure_gc["GVE_p",] - pressure_gc["SAE_p",]))
pgrads_anom <- t(scale(t(pgrads),scale = FALSE))

##### calculate covariance matrices for these files ######
rast_vals <- rbind(temp_deseas_gc, pressure_anom_gc, spatdiff_anom_gc, spatdiff_mil_anom_gc, pgrads_anom)
allnames <- c(rownames(temp_deseas_gc),rownames(pressure_anom_gc),rownames(spatdiff_anom_gc),rownames(spatdiff_mil_anom_gc),rownames(pgrads_anom))

v_vals <- as.matrix(values(terra::as.points(vdat)))
colnames(v_vals) <-  as.character(v_time)
rast_vals_v <- rbind(v_vals,rast_vals)
saveRDS(rast_vals_v,file = "output/rastvals_anom_V_2016-2020.rds",compress = "xz")

u_vals <- as.matrix(values(terra::as.points(udat)))
colnames(u_vals) <-  as.character(u_time)
rast_vals_u <- rbind(u_vals,rast_vals)
saveRDS(rast_vals_u,file = "output/rastvals_anom_U_2016-2020.rds",compress = "xz")

# remove unnecessary stuff
rm(vdat,v_vals,u_vals,udat,u_dat,temp_deseas_gc,spatdiff_anom_gc,spatdiff_gc,spatdiff_mil_gc,spatdiff_mil_anom_gc,rast_vals,pressure_gc,pgrads,pgrads_anom,
   pressure_anom_gc)

############################################################
TOT <- TOThom
wtly <- TRUE
sel <- "U"
samp = TRUE

### calculate covarinace matrix for u and v per WT
if(wtly){
  for(sel in c("U","V")){
    print(sel)
    rast_vals <- readRDS(paste0("output/rastvals_anom_",sel,"_2016-2020.rds"))
    h.i <-  match(allnames,rownames(rast_vals))
    xb <- as.matrix(rast_vals)
    
    ### calculate covariance matrix
    print("calculate covariance matrix by weather type")
    wts <- unique(TOT$WT_type1)
    PH <- array(NA, dim=c(nrow(rast_vals),length(h.i),length(wts)))
    TOTcut <- TOT[TOT$date %in% u_time, ]
    
    for(ww in 1:length(wts)){
      
      wtind = which(TOTcut$WT_type1 == ww)
      if(samp){
        wtind = sample(wtind, size = 122)
      }
      
      xb_mean <- rowMeans(xb[,wtind])
      xb_1_wt <- xb[,wtind] - xb_mean ## calculate anomaly with respect to mean over 50closest analogs
      nens = ncol(xb_1_wt)
      print(paste0(ww, ":",nens))
      
      for (j in seq_along(h.i)){
        PH[,j,ww] <- (xb_1_wt %*% t(xb_1_wt[h.i[j],,drop=F]) / (nens - 1))
      }
    }
    
    dimnames(PH) <- list(rownames(rast_vals),allnames,1:length(wts))
    
    print("write PH")
    filename= paste0("PH_WTnew7_anom_",ifelse(samp,"samp_",""),stringr::str_to_upper(sel),"wind_all.rds")
    saveRDS(PH,file = filename,compress = "xz")
    
  }
}
