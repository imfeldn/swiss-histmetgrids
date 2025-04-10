# 020 Preprocessing
# preprocessing of input data for analogue reconstructions
# 2022-07-29


### input needed:
# specifications from masterscript
# data frame named TOT with first col "date" (date format), second col "weather type"
# and station data with colnames paste("3-LETTER-ABBREVIATION_","ta","p","rr","rr0,sep="")

#############################################################################################################
#############################################################################################################
# 1. PRESELECTION (stations, variables)
# 2. TREND CORRECTION (temperature)
# 3. SEASONALITY CORRECTION (temperature)
# 4. CLIMATE OFFSET FOR REFERENCE STATIONS
# 5. SEASONALITY CORRECTION (temperature)
# 6. STANDARDISE DATA
# 7. CONVERT RR0 DATA TO LOGICAL IN CASE OF GOWER AS DISTANCE METRIC
# 8. CREATE DATA FRAMES FOR PREDICTOR AND PREDICTAND
#############################################################################################################
#############################################################################################################

### load data
(load("data/analogue_stationdata_tmax_tmin.RData"))

### read inventory
inventory <- readWorksheetFromFile("data/station_inventory.xlsx",sheet=1, colTypes="character",startRow=2)

#############################################################################################################
# 1. PRESELECTION (stations, variables, weather types)
#############################################################################################################
print("select stations and variables")

data <- TOT
selcols <- which(colnames(data) %in% scenarios[[specifications$scenario]])
cols <- c(grep("WT|date",colnames(data)),selcols)
data <- data[,cols]

#############################################################################################################
# 2. TREND CORRECTION (temperature)
#############################################################################################################
if ("ta"%in%specifications$variables ) {
  if (specifications$trend!=F) {
    print("removing temperature trend")
    
    t_dates <- data$date
    t_data <- data[,grepl("_ta|tmin|tmax",colnames(data))]
    
    ## coordinates of stations
    ids <- colnames(t_data)
    
    ## select ref period starting from 1979
    refind <- t_dates >= "1971-01-01"
    t_data <- t_data[refind,]
    
    ## load zonal trend
    if ("ta"%in%specifications$variables ){
      print("ta")
      (load("data/era5_zonmean_trend_1971-2020.RData"))
      ref_ta <- ref_trend; rm(ref_trend)
      dnams <- dimnames(ref_ta)
      
      ids_ta <- ids[grepl("_ta", ids)]
      coords_ta <- as.numeric(inventory[match(ids_ta,inventory$ID),c("lat_hist")])
      
      ## here, take the closest latitude for every station + select zonal mean trend
      nearlats <- sapply(coords_ta, function(x){which.min(abs(x-as.numeric(dnams[[1]])))})
      temp <- t_data[,grepl("_ta",colnames(t_data))] - t(ref_ta[nearlats,])
      data[refind,grepl("_ta",colnames(data))] <- temp  
      rm(ids_ta,coords_ta,ref_ta,temp,nearlats,dnams)
      
    }
    
    ## load zonal trend
    if ("tmin"%in%specifications$variables ){
      print("tmin")
      ref_tn <- read.table("data/era5_zonmean_tmin_trend_1971-2020.txt", sep = ",")
      dnams <- dimnames(ref_tn)
      
      ids_tn <- ids[grepl("_tmin", ids)]
      if(length(ids_tn) >0){
        ids_tn <- gsub("tmin","tn",ids_tn)
        coords_tn <- as.numeric(inventory[match(ids_tn,inventory$ID),c("lat_hist")])
        
        ## here, take the closest latitude for every station + select zonal mean trend
        nearlats <- sapply(coords_tn, function(x){which.min(abs(x-as.numeric(dnams[[1]])))})
        temp <- t_data[,grepl("_tmin",colnames(t_data))] - t(ref_tn[nearlats,])
        data[refind,grepl("_tmin",colnames(data))] <- temp      
        rm(ids_tn,coords_tn,ref_tn,temp,nearlats,dnams)
        
      }
    }

    ## load zonal trend
    if ("tmax"%in%specifications$variables ){
      print("tmax")
      ref_tx <- read.table("data/era5_zonmean_tmax_trend_1971-2020.txt", sep = ",")
      dnams <- dimnames(ref_tx)
      
      ids_tx <- ids[grepl("_tmax", ids)]
      if(length(ids_tx) >0){
        ids_tx <- gsub("tmax","tx",ids_tx)
        coords_tx <- as.numeric(inventory[match(ids_tx,inventory$ID),c("lat_hist")])
      
        ## here, take the closest latitude for every station + select zonal mean trend
        nearlats <- sapply(coords_tx, function(x){which.min(abs(x-as.numeric(dnams[[1]])))})
        temp <- t_data[,grepl("_tmax",colnames(t_data))] - t(ref_tx[nearlats,])
        data[refind,grepl("_tmax",colnames(data))] <- temp   
        #plot(1:length(temp$BER_tmax),temp$BER_tmax)
        #points(1:length(temp$BER_tmax),data$BER_tmax[refind], col = "blue")
        rm(ids_tx,coords_tx,ref_tx,temp,nearlats,dnams)
      }
    }
    
    TOT <- data
    save(TOT, file=paste0("data/analogue_stationdata_dtrd1971_refperiod_",scen,"_",date,".RData"))
  } else {
    print("save station data no detrending")
    TOT <- data
    save(TOT, file=paste0("data/analogue_stationdata_",scen,"_",Sys.Date(),".RData"))
  }
  
  #############################################################################################################
  # 4. CLIMATE OFFSET FOR REFERENCE STATIONS
  #############################################################################################################
  if(specifications$offset==T){
    
    print("add offset modera")
    
    t_dates <- data$date
    t_data <- data[,grepl("_ta|tmax|tmin",colnames(data))]
    
    offset <- read.table("code/runclim_offset_year_modera_till1971.txt", header=T)
    offset[is.na(offset)] <- 0
    offset_mat <- matrix(rep(offset[,1], ncol(t_data)),ncol=ncol(t_data), byrow = F)
    yind <- match(year(t_dates), as.numeric(rownames(offset)))
    temp_off <- t_data +  offset_mat[yind,]
    data[,grepl("_ta|tmax|tmin",colnames(data))] <- temp_off
    rm(offset,offset_mat, yind, temp_off, t_dates)
    
    save(TOT, file=paste0("data/analogue_stationdata_offset_dtrd1971_",scen,"_",date,".RData"))
    
    
  }
  
  
  #############################################################################################################
  # 5. SEASONALITY CORRECTION (temperature)
  #############################################################################################################
  if (specifications$seasonality==T) {
    
    print("filtering temperature seasonality")
    
    if(specifications$sep_ref_hist==T){
      
      ind <- list(ref = data$date >= specifications$startRef,
                  hist = data$date <= paste0(specifications$endHist,"-01-01"))
      
      for(per in 1:2){
        t_date <- data$date[ind[[per]]]
        t_data <- data[ind[[per]],grepl("ta",colnames(data))]
        # calculate seasonality-corrected temperature data
        doys <- lubridate::yday(t_date)
        yrs <- lubridate::year(t_date)
        deseas <- apply(t_data,2,rm_seas, doys = doys, yrs = yrs)
        data[ind[[per]],colnames(data)%in%colnames(deseas)] <- deseas
        rm(deseas,t_data,t_date)
      }
      
    } else {
      t_date <- data$date
      t_data <- data[,grepl("ta",colnames(data))]
      
      # calculate seasonality-corrected temperature data
      deseas <- apply(t_data,2,rm_seas, dates=t_date)
      data[,colnames(data)%in%colnames(deseas)] <- deseas
      
      rm(deseas,t_data,t_date)
    }
    
    TOT <- data
  }
}


#############################################################################################################
# 6. STANDARDISE DATA
#############################################################################################################
if(specifications$standardize==T) {
  print("standardizing data")
  selcols <- gsub(".*_","",colnames(data)) %in% c("p","ta","tmin","tmax")
  
  if(specifications$sep_ref_hist==T){
    
    ind <- list(ref = data$date >= specifications$startRef,
                hist = data$date <= paste0(specifications$endHist,"-01-01"))
    
    for(per in 1:2){
      data_stand <- data[ind[[per]],selcols]
      data_stand <- data.frame(scale(data_stand))
      data[ind[[per]],colnames(data)%in%colnames(data_stand)]<-data_stand
    }
    rm(data_stand, per,ind)
    
  } else {
    data_stand <- data[selcols]
    data_stand <- data.frame(scale(data_stand))
    data[,colnames(data)%in%colnames(data_stand)]<-data_stand
    rm(data_stand)
  }
  
  TOT <- data
}


#############################################################################################################
# 7. CONVERT RR0 DATA TO LOGICAL IN CASE OF GOWER AS DISTANCE METRIC
#############################################################################################################

if(specifications$distance=="Gower"){
  if(any(grepl("rr0",colnames(data)))){
    
    dfrr0 <- as.data.frame(data[,grep("rr0",colnames(data))])
    aold <- ifelse(dfrr0[,1]==1, T,F)
    if(dim(dfrr0)[2]>1){
      for(i in 2:ncol(dfrr0)){
        a <- ifelse(dfrr0[,i]==1, T,F)
        aold <- cbind(aold,a)
      }
      rm(a)
    }
    data[,grep("rr0",colnames(data))] <- aold
    rm(aold,dfrr0)
  }
}

#############################################################################################################
# 8. CREATE DATA FRAMES FOR PREDICTOR AND PREDICTAND
#############################################################################################################

# ADD EQUAL PROBABILITIES TO WTS WHERE THEY ARE MISSING
naprob <- which(is.na(data$WT_probability1))[2:length(which(is.na(data$WT_probability1)))]
if(!all(is.na(naprob))){
  data[naprob,grepl("probability", colnames(data))] <- 1/7
}

# define period of predictor & predictand
if(is.na(specifications$selDates[1])){
  predictordates <- seq(as.Date(specifications$startRef),as.Date(specifications$endRef),by="day")  
} else {
  predictordates <- specifications$selDates
}

predictor <- data[data$date%in%predictordates,!grepl("prob",colnames(data))]

# include the ref periods as well, if calonly is FALSE!
predictanddates <- switch(specifications$calonly+1,
                          unique(c(seq(as.Date(specifications$startCal[1]),as.Date(specifications$endCal[1]),by="day"), predictordates)),
                          c(seq(as.Date(specifications$startCal[1]),as.Date(specifications$endCal[1]),by="day")))
predictand <- data[data$date%in%predictanddates,]

print("preselection completed")

# remove variables
rm(predictordates,predictanddates,  selcols)


# END
