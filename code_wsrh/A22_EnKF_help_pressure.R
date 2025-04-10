### calculate bias between stations and grid  ###
rm(list=ls())

library(ncdf4)
library(lubridate)
library(terra)

scen = "scen_hist2_wind"
date = "2025-04-02"
psfile <- "data/cosmo_20162020_daymean_PS_LV95.nc"

### read inventory
source("code/A0_scenarios.R")
station_coords <- XLConnect::readWorksheetFromFile("data/station_inventory.xlsx",sheet=1, startRow=2)[,c("ID","E_lv95_ref","N_lv95_ref")]
station_coords$E_lv95_ref <- as.numeric(station_coords$E_lv95_ref)
station_coords$N_lv95_ref <- as.numeric(station_coords$N_lv95_ref)

#################################################################
### 1. mean monthly difference between cosmo pressure and measurements ###
#################################################################

print("calc diff")
(load(paste0("data/analogue_stationdata_offset_dtrd1979_",scen,"_",date,".RData")))

selstns = grepl("_p",colnames(TOT)) & !grepl("_prob",colnames(TOT))
station_data <- data.frame("date"=as.Date(TOT$date),TOT[,selstns])
stnnam <- colnames(station_data[-1])
coords_sub <- station_coords[match(stnnam,station_coords$ID),]

psdat <- terra::rast(psfile, "PS")
ps_time <- as.Date(time(psdat))

# get the equivalent grid cells
stat <- cellFromXY(psdat, cbind(unlist(coords_sub[,2]), unlist(coords_sub[,3])))
names(stat) <- coords_sub$ID
test <- t(terra::extract(x = psdat, y = stat)/100)
colnames(test) <- coords_sub$ID

## add station data for areas where there is no grid cell
nastns <- which(apply(test,2, function(x) all(is.na(x))))
tind <- which(TOT$date %in% ps_time)
for(ss in 1:length(nastns)){test[,nastns[ss]] <- TOT[tind,names(nastns)[ss]]}

## calculate monthly differences between 
tind <- station_data$date %in% ps_time
diffvalues <- station_data[tind,-1] - test
mondiff <- sapply(1:ncol(diffvalues), function(x) {aggregate(diffvalues[,x], by = list(lubridate::month(ps_time)), mean, na.rm = T)$x})
colnames(mondiff) <- colnames(test)

write.table(mondiff, file=paste0("data/diff_station_raster_all_p.txt"), row.names = F, col.names = T)

## calculate anomalies of the selected grid cells/filled station data
## write an RDS file with the extract pressure data from the grid cells
pressure_anom <- t(scale(test, scale = FALSE))
colnames(pressure_anom) <- ps_time
saveRDS(pressure_anom, file = "data/pressure_anom_station_cosmo.rds")  
saveRDS(t(test), file = "data/pressure_station_cosmo.rds")  
