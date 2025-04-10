### Calculate bias for tmin, tmax, and tabs; run manually for tmin/tmax
rm(list = ls())

library(ncdf4)
library(lubridate)
library(terra)

## Define file names
tasel <- "tmin"  # select "tmin" or "tmax"
varnam <- ifelse(tasel == "tmin", "TminD", "TmaxD")
tempfile <- "data/TabsD_ch01r.swiss.lv95_1971-2020.nc"
tmaxfile <- paste0("data/", tasel, "/", varnam, "_ch01r.swiss.lv95_detrended_1971-01-01-2020-12-31.nc")
era5trend <- "data/era5_zonmean_trend_1971-2020.RData"
offset <- TRUE

### Read inventory
station_coords <- XLConnect::readWorksheetFromFile("data/station_inventory.xlsx", sheet = 1, startRow = 2)[, c("ID", "E_lv95_ref", "N_lv95_ref", "lat_hist")]
station_coords$E_lv95_ref <- as.numeric(station_coords$E_lv95_ref)
station_coords$N_lv95_ref <- as.numeric(station_coords$N_lv95_ref)
rownames(station_coords) <- station_coords$ID

#################################################################
### 1. Mean monthly difference between TabsD and Measurements ###
#################################################################

print("calc diff")

# Load analogue station data
load(paste0("data/analogue_stationdata_offset_dtrd1971_", scen, "_", date, ".RData"))

# Select temperature columns
selstns <- grepl("_ta", colnames(TOT)) & !grepl("_prob", colnames(TOT)) 
station_data <- data.frame("date" = as.Date(TOT$date), TOT[, selstns])
stnnam <- colnames(station_data[-1])

# Subset coordinates
coords_sub <- station_coords[match(stnnam, station_coords$ID), ]
coords_lat <- as.numeric(coords_sub$lat_hist)

# Load TabsD raster
tempdat <- terra::rast(tempfile, "TabsD")
temp_time <- as.Date(time(tempdat))

# Extract raster cells for stations
stat <- cellFromXY(tempdat, cbind(unlist(coords_sub[, 2]), unlist(coords_sub[, 3])))
names(stat) <- coords_sub$ID
test <- t(terra::extract(x = tempdat, y = stat))
colnames(test) <- coords_sub$ID

# Add dates
test <- cbind(temp_time, test)

# Load and apply zonal mean trend
load(era5trend)
dnams <- dimnames(ref_trend)

# Remove trend based on latitude
nearlats <- sapply(coords_lat, function(x) { which.min(abs(x - as.numeric(dnams[[1]]))) })
test_dtr <- cbind(test[, 1], test[, -1] - t(ref_trend[nearlats, ]))

# Restore NA stations
nastns <- which(apply(test_dtr, 2, function(x) all(is.na(x))))
tind <- which(TOT$date %in% temp_time)
for (ss in 1:length(nastns)) {
  test_dtr[, nastns[ss]] <- TOT[tind, names(nastns)[ss]]
}

# Calculate monthly differences
tind <- station_data$date %in% temp_time
diffvalues <- station_data[tind, -1] - test_dtr[, -1]
mondiff <- sapply(1:ncol(diffvalues), function(x) {
  aggregate(diffvalues[, x], by = list(lubridate::month(temp_time)), mean, na.rm = TRUE)$x
})
colnames(mondiff) <- colnames(test)[-1]

write.table(mondiff, file = "data/diff_station_raster_tabsd.txt", row.names = FALSE, col.names = TRUE)

# Deseasonalize temperature data
temp_deseas <- t(apply(test_dtr, 2, rm_seas2, dates = temp_time))
colnames(temp_deseas) <- temp_time
rownames(test_dtr) <- temp_time

saveRDS(temp_deseas, file = "data/temp_deseas_station_tabs.rds")
saveRDS(t(test_dtr), file = "data/temp_station_tabs.rds")

terra::free_RAM()

########################################################################
### 2. Calculate differences with respect to tmin/tmax raster data ###
########################################################################

# Load raster and station data
tempdat <- terra::rast(tmaxfile, varnam)
temp_time <- as.Date(time(tempdat))

selstns <- grepl(paste0("_", tasel), colnames(TOT)) & !grepl("_prob", colnames(TOT))
station_data <- data.frame("date" = as.Date(TOT$date), TOT[, selstns])
stnnam <- colnames(station_data[-1])
stnnam <- gsub(tasel, ifelse(tasel == "tmax", "tx", "tn"), stnnam)

coords_sub <- station_coords[stnnam, ]

# Extract raster cells for stations
stat <- cellFromXY(tempdat, cbind(unlist(coords_sub[, 2]), unlist(coords_sub[, 3])))
names(stat) <- coords_sub$ID
test <- t(terra::extract(x = tempdat, y = stat))
colnames(test) <- coords_sub$ID

# Add dates
test <- cbind(temp_time, test)

# Calculate monthly differences
tind <- station_data$date %in% temp_time
diffvalues <- station_data[tind, -1] - test[, -1]
mondiff <- sapply(1:ncol(diffvalues), function(x) {
  aggregate(diffvalues[, x], by = list(lubridate::month(temp_time)), mean, na.rm = TRUE)$x
})
colnames(mondiff) <- colnames(test)[-1]

write.table(mondiff, file = paste0("data/diff_station_raster_", tasel, "g_", tasel, "s.txt"), row.names = FALSE, col.names = TRUE)
