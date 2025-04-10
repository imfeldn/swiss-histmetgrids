### Ensemble Kalman fitting (temperature)


#############################################################################################################
# 1. READ DATA
#    1.1 READ ANALOGUE DATES
#    1.2 READ RASTER DATA
#    1.3 READ STATION DATA
# 2. KALMAN FITTING
#    2.1 DEFINE VALUES
#    2.2 CALCULATE DISTANCE MATRIX
#    2.3 Kalman Filter for each day in analogue dates
# 3. WRITE RESULTS TO FILES
#############################################################################################################
library(lubridate)
library(ncdf4)
library(terra)

rm(list=ls())
RhpcBLASctl::blas_set_num_threads(1)

### the following lines take arguments (year) from console
ptm <- proc.time()

args <- commandArgs(TRUE)
if (length(args)>0) {
  in_year = as.numeric(args[1])
  scen = args[2]
  date = args[3]
  startCal = args[4]
  endCal = args[5]
  L = as.numeric(args[6])
  Z = as.numeric(args[7])
  climoff = args[8]
  dist = args[9]
  window = args[10]
  tasel = args[11]
  hist = args[12]
}

print(in_year)
print(scen)

source("code_txtnsd/A02_scenarios_txtn.R")
source("code_txtnsd/A03_helpfun.R")
detr = T

## adjust error list based on your errors
errors_list = list("hist1" = "txtn_all_temp_txtn_1763-01-01_1863-12-31",
                   "hist2" = "txtn_all_temp_txtn_1864-01-01_2020-12-31")

varnam <- ifelse(tasel == "tmax", "TmaxD","TminD")

###################################################################################################
### KALMAN FITTING
###################################################################################################
# define setting for Kalman fitting
kalmanset <- list(
  analog_dates =        paste0("output/ARM_run_",date,"/",scen,"_",window,"/",date,"_analog.dates_",startCal,"-",endCal,"_deseas",switch(detr+1,"_","_detrend_"),dist,"_",scen,"_",window,".txt"),
  station_data =        paste0("data/analogue_stationdata_",ifelse(climoff=="T","offset_",""),"dtrd1971_",scen,"_",date,".RData"),        
  input_raster_deseas =       paste0("data/",tasel,"/",varnam,"_ch01r.swiss.lv95_detrended_deseas_1971-01-01-2020-12-31.nc"), ## tmax tmin data from MeteoSwiss
  input_climatology =       paste0("data/",tasel,"/",varnam,"_ch01r.swiss.lv95_detrended_clim_1971-01-01-2020-12-31.nc"), ## climatology calculated from tmax tmin data from MeteoSwiss
  offset_file =         "data/runclim_offset_year_modera_till1971.txt",
  stationset = scen,
  Rref = hist, # define measurement error R either as value [°K] or calculated based on errors_list (line73)
  da_type = c("ta"),
  covmat = "members",
  scaling = "anom",
  nens = 50,
  restrictions = TRUE,
  localize = FALSE
)


if (!exists("in_year")) message("define year for which EnKF should be run")
stopifnot(exists("in_year"))

## create output folder for validation
rdir <- paste0("R",kalmanset$Rref)
loc <- ifelse(kalmanset$localize,paste0(kalmanset$cutoff,"km_",kalmanset$Z_cutoff,"m_",rdir,"_"),paste0(rdir,"_"))
subDir <- paste0("EnKF_",tasel,"_",paste0(kalmanset$da_type,collapse = "_"),"_rmseassep_stns_as_analogues_",ifelse(kalmanset$scaling!="",paste0(kalmanset$scaling,"_"),""),kalmanset$covmat,"_",loc,Sys.Date())
dir.create(file.path(dirname(kalmanset$analog_dates),subDir))
outfolder <- paste(file.path(dirname(kalmanset$analog_dates),subDir),"/",sep="")
rm(subDir)

#############################################################################################################
### 1. READ DATA
#############################################################################################################

#----------------------------------------------------------------------------------------
### 1.1 READ ANALOGUE DATES
#----------------------------------------------------------------------------------------
analog.dates <- read.table(kalmanset$analog_dates,as.is=T)

#----------------------------------------------------------------------------------------
### 1.2 READ RASTER DATA
#----------------------------------------------------------------------------------------
tmax_deseas <- terra::rast(kalmanset$input_raster_deseas)
## rename the layers of tmax_deseas
names(tmax_deseas) <- terra::time(tmax_deseas)
tmax_climatology <- terra::rast(kalmanset$input_climatology)
climatology_vals <- as.matrix(terra::values(terra::as.points(tmax_climatology)))

if(climoff){
  offset <- read.table(kalmanset$offset_file)
}

restricted = paste0(c("ALT","BAS","BER","CDF","CHM","GRH","LUZ","GSB","GVE","OTL","SAE","SMA","STG","GRC","BUS","VEV","SHA",
                          "DAV","ELM", "LUG", "SIA"),"_ta")

## restrict number of tmean assimilations
if(kalmanset$restrictions){
  stationset <- scenarios[[kalmanset$stationset]]
  ## if ta and not in restriction remove from stationset
  stationset <- stationset[stationset %in% restricted | grepl(tasel,stationset)]
  stnsel <- paste0(paste0("_",kalmanset$da_type), collapse = "|")
  stationset <- stationset[grepl(stnsel,stationset)]
} else {
  stationset <- scenarios[[kalmanset$stationset]]
  stnsel <- paste0(paste0("_",kalmanset$da_type), collapse = "|")
  stationset <- stationset[grepl(stnsel,stationset)]
}

### add station set again to the info-list, to display them in the readme file
kalmanset$stationset <- stationset

#----------------------------------------------------------------------------------------
### 1.3 READ STATION DATA
#----------------------------------------------------------------------------------------
### Read station measurements (only trend removed, nothing else)
(load(kalmanset$station_data)) # this is called TOT
station_data_all <- data.frame("date"=as.Date(TOT$date),TOT[,stationset])
station_coords <- XLConnect::readWorksheetFromFile(kalmanset$station_xy,sheet=1, colTypes="character",startRow=2)[,c("ID","lon_ref","lat_ref","E_lv95_ref","N_lv95_ref")]
station_coords$ID <- ifelse(grepl("tx",station_coords$ID), gsub("tx","tmax", station_coords$ID), gsub("tn","tmin", station_coords$ID))
station_coords <- station_coords[match(stationset,station_coords$ID),]

allmon = month(station_data_all$date)
### Remove Bias from tmax observations, for tmean not needed 
if(tasel %in% kalmanset$da_type ){
  diff <- read.table(paste0("data/diff_station_raster_",tasel,"g_",tasel,"s.txt"), head=T)
  colnames(diff) <- gsub("tx|tn",tasel, colnames(diff))
  id_stn <- which(colnames(station_data_all) %in% colnames(diff))
  
  for(j in colnames(station_data_all)[id_stn]){
    print(j)
    for(i in c(1:12)){
      station_data_all[allmon == i,j] <- station_data_all[allmon == i,j] - diff[i,j]
    }
  }
}
## remove seasonality from all stations for two periods separately
p1 <- station_data_all$date < "1971-01-01"
p2 <-station_data_all$date >= "1971-01-01"
all_da <- grepl(paste(kalmanset$da_type,collapse ="|"),colnames(station_data_all))
doys = lubridate::yday(station_data_all$date)
yrs =  lubridate::year(station_data_all$date)
station_data_all[p1,all_da] <- apply(X = station_data_all[p1,all_da],MARGIN = 2,FUN = rm_seas, doys=doys[p1],yrs= yrs[p1])
station_data_all[p2,all_da] <- apply(X = station_data_all[p2,all_da],MARGIN = 2,FUN = rm_seas, doys=doys[p2],yrs= yrs[p2])
rownames(station_data_all) <- station_data_all$date

### select assimilation year
station_data <- station_data_all[which(yrs%in%in_year),]
analog.dates <- analog.dates[which(substr(analog.dates$date,1,4)%in%in_year),]

### Load a temporary file
targetdays <- analog.dates$date
tmax_kalman <- tmax_deseas[[1:nrow(analog.dates)]]
names(tmax_kalman) <- targetdays
tmax_deseas_kalman <- tmax_kalman

#############################################################################################################
### 2. KALMAN FITTING
#############################################################################################################

#----------------------------------------------------------------------------------------
### 2.1 DEFINE VALUES
#----------------------------------------------------------------------------------------
### define measurement error 
print("create observation error covariance matrix")
if(is.numeric(kalmanset$Rref)){
  R <- matrix(rep(kalmanset$Rref,length(stationset)), ncol=1)
  rownames(R) <- stationset
} else{
  Rhelp <- list()
  for(var in kalmanset$da_type){
    print(var)
    error_path = paste0( "data/obserrors/obserror_",errors_list[[kalmanset$Rref]])
    Rhelp[[var]] <- data.frame(read.table(paste0(error_path,"_",var,"_",date,".csv"), sep = ","))
  }
  R <- matrix(unlist(sapply(Rhelp, function(x) x[,2])),ncol = 1)
  rownames(R) <- as.character(unlist(sapply(Rhelp, function(x) x[,1])))
}

#----------------------------------------------------------------------------------------
### 2.3 Kalman Filter for each day in analogue dates
#----------------------------------------------------------------------------------------
day <- 1
nens <- kalmanset$nens

for (day in 1:nrow(analog.dates)){
  
  print(analog.dates$date[day])
  
  ### subset raster for analog days
  analogdays_date <- as.Date(as.numeric(analog.dates[day,-1]),origin = "1970-01-01")[1:nens]
  
  if(all(is.na(analogdays_date))){
    values(tmax_kalman[[day]])[!is.na(values(tmax_kalman[[day]]))]<- NA
    values(tmax_deseas_kalman[[day]])[!is.na(values(tmax_deseas_kalman[[day]]))]<- NA
  } else {
    
    ## select analog values, remove NAs
    analogdays <- as.character(analogdays_date)
    tmax_analogs_deseas <- terra::subset(tmax_deseas, as.character(analogdays_date))
    
    doy <- paste0(varnam,"_",lubridate::yday(analog.dates$date[day]))
    
    ## change climatology depending on period
    clim <- climatology_vals[,doy]
    
    if(in_year < 1971){
      print("remove offset")
      clim <- climatology_vals[,doy] - offset[as.character(in_year),"x"]
    }
    
    ## create data matrix of analog days and add climatology
    rast_vals <- as.matrix(terra::values(terra::as.points(tmax_analogs_deseas)))

    ### subset station data and create observation vector
    station_data_sub <- station_data[station_data$date==analog.dates$date[day],c("date",stationset)]
    station_data_sub <- station_data_sub[,!is.na(station_data_sub)]
    obs <- station_data_sub[,2:ncol(station_data_sub)]
    
    ## add temp measurements using observations
    temp_obs <- t(station_data_all[as.character(analogdays_date),stationset[grepl("_ta",stationset)]])
    
    if(nrow(temp_obs)!=0){
      nastns <- apply(temp_obs, 1, function(x){any(is.na(x))})
      if(any(nastns)){
        temp_obs <- temp_obs[-which(nastns),]
        stationset <- stationset[-which(nastns)]
      } 
      rast_vals <- rbind(rast_vals,temp_obs)
    }
    
    ## select stations
    obs <- obs[colnames(obs) %in% rownames(temp_obs) | grepl(tasel, colnames(obs))]
    stnnam <- colnames(station_data_sub[2:ncol(station_data_sub)])
    stnnam_tasel <- stnnam[grepl(tasel, stnnam)]
    coords_sub <- station_coords[match(stnnam_tasel,station_coords$ID),]
    
    # get stations within tmax grid
    stat <- cellFromXY(tmax_analogs_deseas[[1]], cbind(as.numeric(unlist(coords_sub[,4])), as.numeric(unlist(coords_sub[,5]))))
    cell.nrs <- which(!is.na(terra::values(tmax_analogs_deseas[[1]])))
    
    ## create H vector with all stations or only inside
    if(nrow(temp_obs)!=0){
      outside <- match(colnames(obs),rownames(rast_vals))
      out_ind <- !is.na(outside)
      h.i <- c(outside[out_ind], match(stat,cell.nrs))
      names(h.i) <- c(colnames(obs)[out_ind], coords_sub$ID)
    } else{
      h.i <- match(stat,cell.nrs)
      names(h.i) <- coords_sub$ID
      outside <- c() ## must be empty
    }
    
    obsind <- !is.na(h.i)
    h.i <- h.i[obsind]
    ## order hi and obs so they are in aggreement
    obs <- obs[names(h.i)]
    H <- t(sapply(h.i, function(i) {out <- rep(0, nrow(rast_vals)) ; out[i] <- 1; out}))
    
    ### create observation covariance matrix
    Rcal <- R[names(h.i),]
    
    ### define background mean and variance
    xb <- as.matrix(rast_vals)
    xb_mean <- rowMeans(xb)
    xb_1 <- xb - xb_mean ## calculate anomaly with respect to mean over 50closest analogs
    
    ### calculate covariance matrix
    print("select covariance matrix")
    
    if (kalmanset$covmat == "members"){
      print("calculate covariance matrix based on 50 members")
      PH <- matrix(ncol=length(h.i),nrow = nrow(rast_vals))
      for (j in seq_along(h.i)){
        PH[,j] <- (xb_1 %*% t(xb_1[h.i[j],,drop=F]) / (nens - 1))
      }
      colnames(PH) <- names(obs)
      rownames(PH) <- rownames(rast_vals)
    } 
    # select needed gridcells
    PH <- PH[,colnames(obs)]
    
    ### do Kalman fitting
    print("apply Kalman fitting")
    HPHR <- (H %*% PH + diag(Rcal))
    K <- PH %*% solve(HPHR)
    eHPHR <- eigen(HPHR)
    sqrtHPHR <- eHPHR$vectors %*% diag(sqrt(eHPHR$values)) %*% t(eHPHR$vectors)
    Ktilde <- PH %*% t(solve(sqrtHPHR)) %*% solve(sqrtHPHR + diag(sqrt(Rcal)))
    
    ## update for mean
    xa <- as.vector(xb_mean + K %*% (as.numeric(obs[obsind]) - (H %*% xb_mean)))
    
    
    print("update anom")
    ## update for anomaly
    xb_1_small <- xb_1[h.i,]
    a_small <- Ktilde %*% H[,h.i]     
    b <- a_small %*% xb_1_small       
    xa1 <- xb_1 - b
    analysis <- xa + xa1
    
    print("analysis")
    ## remove the cells of the pressure information
    analysis_deseas <- analysis[1:(nrow(rast_vals) - nrow(temp_obs)),]
    analysis <- apply(analysis_deseas,2,function(x) x+clim) 
    
    ### write data to raster files (temperature & deseas temperature)
    terra::values(tmax_kalman[[day]])[!is.na(terra::values(tmax_kalman[[day]]))]<- as.numeric(analysis[,1])
    terra::values(tmax_deseas_kalman[[day]])[!is.na(terra::values(tmax_deseas_kalman[[day]]))]<- as.numeric(analysis_deseas[,1])

  }
  
}


#############################################################################################################
### 3. WRITE RESULTS TO FILES
#############################################################################################################
sdate <- as.character(analog.dates[1,1])
edate <- as.character(analog.dates[nrow(analog.dates),1])

if (kalmanset$localize==T) {
  props1 <- "_localized"
  props2 <- paste("_",kalmanset$cutoff,"km",sep="")
  props3 <- paste("_",kalmanset$Z_cutoff,"m",sep="")
} else {
  props1 <- c()
  props2 <- c()
  props3 <- c()
}

# write gridded files
print("write raster")

gridmap = c("grid_mapping_name=Oblique Mercator (LV95 - CH1903+)", 
            "longitude_of_projection_center=7.43958333", 
            "latitude_of_projection_center=46.9524056", 
            "false_easting=2600000.", 
            "false_northing=1200000.",
            "inverse_flattening=299.1528128",
            "semi_major_axis=6377397.155")

#rs = terra::rast(raster::stack(ws_kalman))
tt <- seq(as.Date(sdate),as.Date(edate),by="day")
terra::time(tmax_kalman) <- tt
terra::time(tmax_deseas_kalman) <- tt
longname <- ifelse(tasel == "tmax","Daily maximum temperature","Daily minimum temperature")

coords <- "+proj=somerc +lat_0=46.9524055555556 +lon_0=7.43958333333333 +k_0=1 +x_0=2600000 +y_0=1200000 +ellps=bessel +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs +type=crs"
terra::crs(tmax_kalman) <- coords
terra::crs(tmax_deseas_kalman) <- coords

terra::writeCDF(tmax_kalman, 
                filename = paste(outfolder,Sys.Date(),"_",tasel,"_EnKF_",sdate,"-",edate,props1,props2,props3,".nc",sep=""), 
                varname = tasel,
                longname = longname, 
                unit = "Celsius",
                missval = NA,
                prec= "float",
                gridmap = gridmap,
                overwrite = TRUE,
                compression = 9
)

terra::writeCDF(tmax_deseas_kalman, 
                filename = paste(outfolder,Sys.Date(),"_",tasel,"_deseas_EnKF_",sdate,"-",edate,props1,props2,props3,".nc",sep=""), 
                varname = tasel,
                longname = longname, 
                unit = "Celsius",
                missval = NA,
                prec= "float",
                gridmap = gridmap,
                overwrite = TRUE,
                compression = 9,
                
)

readme<-data.frame("description"=c("kalman filter","period","analog dates","station data","station coordinates","input_raster_deseas","input_raster_climatology","offset_file","stationset",
                                   "Robs","da_type","covmat","scaling","nens"),"parameters"=NA)

readme[1,2] <- as.character(as.Date(Sys.Date()))
readme[2,2] <- paste(sdate,edate,sep="-")
for(i in 3:(nrow(readme)-1)) {readme[i,2]<-as.character((paste(kalmanset[[i-2]], collapse = ", ")))}

write.table(readme,file=paste(outfolder,Sys.Date(),"_readme_EnKF_",sdate,"-",edate,props1,props2,props3,".txt",sep=""))

mtm <- proc.time()
print(paste("calculation time EnKF for ", in_year,":",(mtm[3]-ptm[3])%/%60,"minutes", round((mtm[3]-ptm[3])%%60),"seconds" ,sep=" "))
ptm <- mtm

# remove variables from environment
if(kalmanset$localize==T){rm(alt_mat,dist_mat,L,Z)}

# END

