
rm(list=ls())

library(ncdf4)
library(lubridate)
library(terra)

scen = "scen_hist2_wind"
date = "2025-04-02"
tempfile <- "data/cosmo_20162020_daymean_temp_LV95.nc"

station_coords <- XLConnect::readWorksheetFromFile("data/station_inventory.xlsx",sheet=1, startRow=2)[,c("ID","E_lv95_ref","N_lv95_ref")]
station_coords$E_lv95_ref <- as.numeric(station_coords$E_lv95_ref)
station_coords$N_lv95_ref <- as.numeric(station_coords$N_lv95_ref)

#################################################################
### 1. mean monthly difference between TabsD and Measurements ###
#################################################################

print("calc diff")

(load(paste0("data/analogue_stationdata_detrended_refperiod_",scen,"_",date,".RData")))
selstns = grepl("_ta",colnames(TOT)) & !grepl("_prob",colnames(TOT)) 
station_data <- data.frame("date"=as.Date(TOT$date),TOT[,selstns])
stnnam <- colnames(station_data[-1])
coords_sub <- station_coords[match(stnnam,station_coords$ID),]

tempdat <- terra::rast(tempfile, "T_2M")
temp_time <- as.Date(time(tempdat))

# get the equivalent grid cells
stat <- cellFromXY(tempdat, cbind(unlist(coords_sub[,2]), unlist(coords_sub[,3])))
names(stat) <- coords_sub$ID
test <- t(terra::extract(x = tempdat, y = stat) - 273.15)
colnames(test) <- coords_sub$ID

nastns <- which(apply(test,2, function(x) all(is.na(x))))
tind <- which(TOT$date %in% temp_time)
for(ss in 1:length(nastns)){test[,nastns[ss]] <- TOT[tind,names(nastns)[ss]]}

## calculate monthly differences between 
tind <- station_data$date %in% temp_time
diffvalues <- station_data[tind,-1] - test
mondiff <- sapply(1:ncol(diffvalues), function(x) {aggregate(diffvalues[,x], by = list(lubridate::month(temp_time)), mean, na.rm = T)$x})
colnames(mondiff) <- colnames(test)

write.table(mondiff, file=paste0("data/diff_station_raster_all_temp.txt"), row.names = F, col.names = T)

## deasonalize this temperature data
## write an RDS file with the extract temperature data from the grid cells
temp_deseas <- t(apply(test,2,rm_seas2, dates = temp_time))
colnames(temp_deseas) <- temp_time
saveRDS(temp_deseas, file = "data/temp_deseas_station_cosmo.rds")  

