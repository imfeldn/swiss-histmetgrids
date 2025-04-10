### Ensemble Kalman fitting 

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

RhpcBLASctl::blas_set_num_threads(1)

#library(raster)
library(lubridate)
library(ncdf4)

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
  UVsel = args[10]
  window = args[11]
}


print(UVsel)
print(in_year)
print(scen)

source("code_wsrh/A02_scenarios.R")
source("code_wsrh/A03_helpfun.R")
detr = T

## don't change this list, just add items if changing something
errors_list = list("hist1" = paste0(scen,"_1763-01-01_1863-12-31"),
                   "hist2" = paste0(scen,"_1864-01-01_2020-12-31"))

###################################################################################################
### KALMAN FITTING
###################################################################################################
# define setting for Kalman fitting
kalmanset <- list(
  analog_dates =        paste0("output/ARM_run_",date,"/",scen,"_",window,"/",date,"_analog.dates_",startCal,"-",endCal,"_deseas",switch(detr+1,"_","_detrend_"),dist,"_",scen,"_",window,".txt"),
  station_data =        paste0("data/analogue_stationdata_offset_dtrd1979_",scen,"_",date,".RData"),
  station_xy =          "data/station_inventory.xlsx", # define station coordinate file
  input_example =       "data/cosmo_20162020_daymean_U_10M_LV95.nc",
  stationset = scen,
  localize = F,
  cutoff = L, # define horizontal cutoff distance L [km] (if localize = T)
  Z_cutoff = Z, # define vertical cutoff distance z [m]  (if localize = T)
  Rref = "hist1", # define measurement error R either as value [°K] or calculated based on errors_list (line73)
  da_type = c("ta","p","sp","spm"),
  covmat = "blend",
  scaling = "anom",
  offset = climoff,
  nens = 20,
  restrictions = FALSE,
  climweight = 0.5 #weight of PHclim
)


if (!exists("in_year")) message("define year for which EnKF should be run")
stopifnot(exists("in_year"))

## create output folder for validation
rdir <- paste0("R",kalmanset$Rref)
loc <- ifelse(kalmanset$localize,paste0(kalmanset$cutoff,"km_",kalmanset$Z_cutoff,"m_",rdir,"_"),paste0(rdir,"_"))
subDir <- paste0("EnKF_windspeed_",paste0(kalmanset$da_type,collapse = "_"),"_",ifelse(kalmanset$scaling!="",paste0(kalmanset$scaling,"_"),""),
                 ifelse(kalmanset$offset == "T","offset_",""),
                 kalmanset$covmat,kalmanset$climweight,"_",loc,Sys.Date())
outfolder <- paste(file.path(dirname(kalmanset$analog_dates),subDir),"/",sep="")
dir.create(file.path(dirname(kalmanset$analog_dates),subDir))

rm(subDir)

#############################################################################################################
### 1. READ DATA
#############################################################################################################

#----------------------------------------------------------------------------------------
### 1.1 READ ANALOGUE DATES
#----------------------------------------------------------------------------------------
### load analogue dates
analog.dates <- read.table(kalmanset$analog_dates,as.is=T)
#----------------------------------------------------------------------------------------
### 1.2 READ RASTER DATA
#----------------------------------------------------------------------------------------
restricted_stns = FALSE

if(kalmanset$restrictions){
  stationset <- scenarios[[kalmanset$stationset]]
  stationset <- stationset %in% restricted_stns
  stnsel <- paste0(paste0("_",kalmanset$da_type), collapse = "|")
  stationset <- stationset[grepl(stnsel,stationset)]
  
} else {
  stationset <- scenarios[[kalmanset$stationset]]
  stnsel <- paste0(paste0("_",kalmanset$da_type), collapse = "|")
  stationset <- stationset[grepl(stnsel,stationset)]
  
}

### load rast value data and select relevant varibaels
input_all <- paste0("output/rastvals_",ifelse(kalmanset$scaling!="",paste0(kalmanset$scaling,"_"),""),UVsel,"_2016-2020.rds")
rast_val_all <- readRDS(input_all)
rast_val_all <- rast_val_all[rownames(rast_val_all) %in% c("",stationset),]

#----------------------------------------------------------------------------------------
### 1.3 READ STATION DATA
#----------------------------------------------------------------------------------------

### Read station measurements (only trend removed, nothing else)
load(kalmanset$station_data) # this is called TOT
## GSB as an error, set these to NA
station_data <- data.frame("date"=as.Date(TOT$date),TOT[,stationset])

### Temperature observations
if("ta" %in%  kalmanset$da_type){
  ## remove bias
  diff <- read.table(paste0("data/diff_station_raster_all_temp.txt"), head=T)
  id <- match(colnames(station_data),colnames(diff))
  id <- id[!is.na(id)]
  id_stn <- match(colnames(diff), colnames(station_data))
  id_stn <- id_stn[!is.na(id_stn)]
  
  for(j in 1:length(id_stn)){
    for(i in c(1:12)){
      station_data[month(station_data$date)==i,id_stn[j]] <- station_data[month(station_data$date)==i,id_stn[j]] - diff[i,id[j]]
    }
  }
  
  ## remove seasonality
  doys = lubridate::yday(station_data$date)
  yrs= lubridate::year(station_data$date)
  station_data[,id_stn] <- apply(X = station_data[,id_stn],MARGIN = 2,FUN = rm_seas, doys=doys, yrs=yrs)
  
}

### Pressure observations
if("p" %in%  kalmanset$da_type){
  ## remove bias
  diff <- read.table(paste0("data/diff_station_raster_all_p.txt"), head=T)
  id <- match(colnames(station_data),colnames(diff))
  id <- id[!is.na(id)]
  id_stn <- match(colnames(diff), colnames(station_data))
  id_stn <- id_stn[!is.na(id_stn)]
  
  for(j in 1:length(id_stn)){
    for(i in c(1:12)){
      station_data[month(station_data$date)==i,id_stn[j]] <- station_data[month(station_data$date)==i,id_stn[j]] - diff[i,id[j]]
    }
  }
  ## standardize
  station_data[,id_stn] <- sweep(station_data[,id_stn],2,colMeans(station_data[,id_stn], na.rm = T))
}

### spatial pressure observations
if(any(c("sp","spm") %in%  kalmanset$da_type)){
  id_stn <- grepl("sp",colnames(station_data))
  station_data[,id_stn] <- sweep(station_data[,id_stn],2,colMeans(station_data[,id_stn], na.rm = TRUE))
}

### select assimilation year
station_data <- station_data[which(year(station_data$date)%in%in_year),]
analog.dates <- analog.dates[which(substr(analog.dates$date,1,4)%in%in_year),]

### Load a temporary file
targetdays <- paste("X",gsub("-",".",analog.dates$date),sep="")
suppressWarnings(ws_kalman <- terra::rast(kalmanset$input_example)[[1:nrow(analog.dates)]])
names(ws_kalman) <- targetdays

#############################################################################################################
### 2. KALMAN FITTING
#############################################################################################################

#----------------------------------------------------------------------------------------
### 2.1 DEFINE VALUES
#----------------------------------------------------------------------------------------
### define horizontal (L) and vertical (Z) localisation parameters (sigma in km (L) and m (Z))
if (kalmanset$localize == T)  {
  L = kalmanset$cutoff
  Z = kalmanset$Z_cutoff
}

### define measurement error and cov matrix
print("create observation error covariance matrix")
if(is.numeric(kalmanset$Rref)){
  R <- matrix(rep(kalmanset$Rref,length(stationset)), ncol=1)
  rownames(R) <- stationset
} else{
  Rhelp <- list()
  for(var in kalmanset$da_type){
    error_path = paste0( "data/obserrors/obserror_",errors_list[[kalmanset$Rref]])
    Rhelp[[var]] <- data.frame(read.table(paste0(error_path,"_",var,"_",date,".csv"), sep = ","))
  }
  R <- matrix(unlist(sapply(Rhelp, function(x) x[,2])),ncol = 1)
  rownames(R) <- as.character(unlist(sapply(Rhelp, function(x) x[,1])))
}

print("load model error covariance matrix")
if(kalmanset$covmat == "blend"){
  PHfile <- paste0("output/PH_WTnew7_",ifelse(kalmanset$scaling!="",paste0(kalmanset$scaling,"_"),""),"samp_",UVsel,"wind_all.rds")
  PHdat <- readRDS(PHfile)
  PHdat <- PHdat[,stationset,]
} else if (kalmanset$covmat == "members") {
  PHdat <- NA
  PHfile <- "PH-calculated-from-members"
} else {
  PHfile <- paste0("output/PH_",kalmanset$covmat,"_",ifelse(kalmanset$scaling!="",paste0(kalmanset$scaling,"_"),""),UVsel,"wind_all.rds")
  PHdat <- readRDS(PHfile)
  PHdat <- PHdat[,stationset,]
}

#----------------------------------------------------------------------------------------
### 2.3 Kalman Filter for each day in analogue dates
#----------------------------------------------------------------------------------------
day <- 283
nens <- kalmanset$nens

for (day in 1:nrow(analog.dates)){
  
  print(analog.dates$date[day])
  
  ### subset raster for analog days
  analogdays_date <- as.Date(as.numeric(analog.dates[day,-1]),origin = "1970-01-01")[1:nens]
  
  if(all(is.na(analogdays_date))){
    values(ws_kalman[[day]])[!is.na(values(ws_kalman[[day]]))]<- NA
  } else {
    
    ## select analog values, remove NAs
    analogdays <- as.character(analogdays_date)
    analogdays <- analogdays[!is.na(analogdays)]
    rast_vals <- rast_val_all[,analogdays]
    
    ### subset station data and create observation vector
    station_data_sub <- station_data[station_data$date==analog.dates$date[day],]
    station_data_sub <- station_data_sub[,!is.na(station_data_sub)]
    obs <- station_data_sub[,2:ncol(station_data_sub)]
    
    ### only keep rast_vals where we have an observation
    rast_vals <- rast_vals[rownames(rast_vals) %in% c("",names(obs)),]
    
    ###### create H vector  ######
    h.i <-  match(colnames(obs),rownames(rast_vals))
    obsind <- !is.na(h.i)
    h.i <- h.i[obsind]
    names(h.i) <- colnames(obs)
    H <- t(sapply(h.i, function(i) {out <- rep(0, nrow(rast_vals)) ; out[i] <- 1; out}))
    
    ### create observation covariance matrix
    Rcal <- R[names(h.i),]
    
    ### define background mean and variance
    xb <- as.matrix(rast_vals)
    xb_mean <- rowMeans(xb)
    xb_1 <- xb - xb_mean ## calculate anomaly with respect to mean over 50closest analogs
    
    ### calculate covariance matrix
    print("select covariance matrix")
    
    if(kalmanset$covmat == "yday"){
      doy <- ifelse(day < 366, day, day-1)
      PH <- PHdat[,,doy]
    } else if (grepl("WT",kalmanset$covmat)){
      wt <- TOT$WT_type1[TOT$date == analog.dates$date[day]]
      PH <- PHdat[,,wt]
    } else if (kalmanset$covmat == "blend"){
      
      wt <- TOT$WT_type1[TOT$date == analog.dates$date[day]]
      PH1 <- PHdat[,,wt]
      PH1 <- PH1[rownames(PH1) %in% rownames(rast_vals),colnames(obs)]
      ## calc second PH
      print("calculate covariance matrix")
      PH2 <- matrix(ncol=length(h.i),nrow = nrow(rast_vals))
      for (j in seq_along(h.i)){
        PH2[,j] <- (xb_1 %*% t(xb_1[h.i[j],,drop=F]) / (nens - 1))
      }
      ## blend them together
      PH <- kalmanset$climweight*PH1 + (1 - kalmanset$climweight)*PH2
    } else if (kalmanset$covmat == "members"){
      print("calculate covariance matrix based on 20 members")
      PH <- matrix(ncol=length(h.i),nrow = nrow(rast_vals))
      for (j in seq_along(h.i)){
        PH[,j] <- (xb_1 %*% t(xb_1[h.i[j],,drop=F]) / (nens - 1))
      }
      colnames(PH) <- stationset
      rownames(PH) <- rownames(rast_vals)
    } else {
      PH <- PHdat
    }
    
    # select needed gridcells
    PH <- PH[rownames(PH) %in% rownames(rast_vals),colnames(obs)]
    
    ### do Kalman fitting
    print("apply Kalman fitting")
    #Rprof("test5.out")
    HPHR <- (H %*% PH + diag(Rcal))
    K <- PH %*% solve(HPHR)
    eHPHR <- eigen(HPHR)
    sqrtHPHR <- eHPHR$vectors %*% diag(sqrt(eHPHR$values)) %*% t(eHPHR$vectors)
    Ktilde <- PH %*% t(solve(sqrtHPHR)) %*% solve(sqrtHPHR + diag(sqrt(Rcal)))
    #Rprof(NULL)
    ## update for mean
    xa <- as.vector(xb_mean + K %*% (as.numeric(obs[obsind]) - (H %*% xb_mean)))
    
    ## update for anomaly
    xb_1_small <- xb_1[h.i,]
    a_small <- Ktilde %*% H[,h.i]     #[erstes Matrixprodukt]
    b <- a_small %*% xb_1_small       #[zweites Matrixprodukt]
    xa1 <- xb_1 - b
    analysis <- xa + xa1
    
    ## remove the cells of the pressure information
    analysis <- analysis[1:(nrow(rast_vals) - length(h.i)),] 
    
    ### write data to raster files (temperature & deseas temperature)
    terra::values(ws_kalman[[day]])[!is.na(terra::values(ws_kalman[[day]]))]<- as.numeric(analysis[,1])
    
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

terra::time(ws_kalman) <- as.Date(analog.dates$date) #seq(as.Date(sdate),as.Date(edate),by="day")

terra::writeCDF(ws_kalman, 
                filename = paste(outfolder,Sys.Date(),"_",UVsel,"wind_EnKF_",sdate,"-",edate,props1,props2,props3,".nc",sep=""), 
                varname = "ws",
                longname = "windspeed", 
                unit = "m/s",
                missval = NA,
                prec= "float",
                gridmap = gridmap,
                overwrite = TRUE,
                compression = 9
)


readme<-data.frame("description"=c("kalman filter","period","analog dates","station data","station coordinates","example_grid","stationset","localize","Z","L",
                                   "Robs","da_type","covmat","scaling","offset","nens","covmatfile","weights"),"parameters"=NA)

readme[1,2] <- as.character(as.Date(Sys.Date()))
readme[2,2] <- paste(sdate,edate,sep="-")
for(i in 3:(nrow(readme)-1)) {readme[i,2]<-as.character((paste(kalmanset[[i-2]], collapse = ", ")))}
readme[nrow(readme),2] <- PHfile

write.table(readme,file=paste(outfolder,Sys.Date(),"_readme_EnKF_",sdate,"-",edate,props1,props2,props3,".txt",sep=""))

mtm <- proc.time()
print(paste("calculation time EnKF for ", in_year,":",(mtm[3]-ptm[3])%/%60,"minutes", round((mtm[3]-ptm[3])%%60),"seconds" ,sep=" "))
ptm <- mtm

# remove variables from environment
if(kalmanset$localize==T){rm(alt_mat,dist_mat,L,Z)}

# END
terra::mem_info(terra::rast())
terra::free_RAM()
