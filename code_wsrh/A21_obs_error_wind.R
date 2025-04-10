### Calculate error covariance for observations
# three files are created based on different periods.

rm(list=ls())

library(geosphere)
library(lubridate)

anomalies = TRUE # use anomalies or absolute values
minlength = 120 # minimum overlapping period in months
startCal = "1961-01-01"
endCal = "2020-12-31"
detr = T
scen = "scen_hist2_wind" # this does not call scenarios anymore, since it just uses all in TOT
date = "2025-04-02"


## Function to calculate anomalies
anom <- function(x, nrmlizd=FALSE) {
  xclim <- aggregate(x[,3], list(x[,1]), mean, na.rm=TRUE)
  stdclim <- aggregate(x[,3], list(x[,1]), sd, na.rm=TRUE)
  for (i in 1:12) {
    if (nrmlizd) {
      x[x[,1]==i,3] <- (x[x[,1]==i,3] - xclim$x[xclim[,1]==i]) / stdclim$x[stdclim[,1]==i]
    } else {
      x[x[,1]==i,3] <- x[x[,1]==i,3] - xclim$x[xclim[,1]==i]
    }
  }
  return(x)
}

station_coords <- XLConnect::readWorksheetFromFile("data/station_inventory.xlsx",sheet=1, startRow=2)[,c("ID","lon_ref","lat_ref")]
load(paste0("data/analogue_stationdata",switch(detr+1,"_","_dtrd1979_"),scen,"_",date,".RData"))
TOT <- TOT[,which(!grepl("WT_",colnames(TOT)))]

var_thresholds <- data.frame("p" = 20,"ta" = 10,"sp" = 35,"spm" = 70)

for(variable in c("p","ta","sp","spm")){
  
  var_th <- as.vector(var_thresholds[variable])
  selection <- paste0("_", variable)
  
  TOT_sel <- TOT[,c(1,grep(selection,colnames(TOT)))]
  if(variable == "sp"){ ## separate selection for sp and sp
    TOT_sel <- TOT_sel[!grepl("spm",colnames(TOT_sel))]
  }
  n <- ncol(TOT_sel)-1
  
  ## Read all data and metadata into a list
  Data_ref  <- list()
  for (i in 1:n) {
    Data_ref[[i]] <- list()
    id <- colnames(TOT_sel)[i+1]
    print(id)
    if (variable %in% c("sp","spm")){
      idsel <- station_coords$ID==gsub("_sp|_spm","_p",id)
    } else{
      idsel <- station_coords$ID==id
    }
    Data_ref[[i]]$coords <- as.numeric(station_coords[idsel,c("lon_ref","lat_ref")])
    yind <- which(year(TOT_sel$date)<=year(endCal) & year(TOT_sel$date)>=year(startCal))
    x <- TOT_sel[yind,c("date",id)]
    Data_ref[[i]]$dates <- ISOdate(year(x$date), month(x$date), 1)
    narows <- sum(!is.na(x[,2]))
    if (anomalies & narows>=minlength) {
      x <- cbind(month(x$date),year(x$date), x[,2]) 
      Data_ref[[i]]$ta <- anom(x,nrmlizd = F)[,3]
    } else Data_ref[[i]]$ta <- x$Value
  }
  
  names(Data_ref) <- colnames(TOT_sel)[2:ncol(TOT_sel)]
  
  ## Calculate distance and correlation for every pair of stations
  cormat_ref <- distmat <- array(dim=c(n,n))
  
  colnames(cormat_ref) <- colnames(TOT_sel)[2:ncol(TOT_sel)]
  rownames(cormat_ref) <- colnames(TOT_sel)[2:ncol(TOT_sel)]
  fitvec_ref <- array(dim=n)
  
  for (i in 1:n) {
    message(i)
    for (j in 1:n) {
      if (i != j) {
        distmat[i,j] <- (distGeo(Data_ref[[i]]$coords, Data_ref[[j]]$coords) / 1000)**2
        if (is.na(distmat[i,j])) distmat[i,j] <- 1e+12
        
        x_ref <- Data_ref[[i]]$ta[which(Data_ref[[i]]$dates %in% Data_ref[[j]]$dates)]
        y_ref <- Data_ref[[j]]$ta[which(Data_ref[[j]]$dates %in% Data_ref[[i]]$dates)]
        
        if (sum(!is.na(x_ref-y_ref)) >= minlength & distmat[i,j]<1000**2) {
          if (sd(x_ref,na.rm=TRUE) > 0 & sd(y_ref,na.rm=TRUE) > 0) {
            cormat_ref[i,j] <- var(x_ref-y_ref, use="pairwise.complete.obs")
          }
        }
      }
    }
    # Linear fit (for error estimation after Wartenburger et al.)
    #print(cormat_ref[i,])
    y_ref <- cormat_ref[i,]
    y_ref[which(y_ref>var_th[1])] <- NA
    print(y_ref)
    
    if (sum(!is.na(y_ref)) >= 5) {
      fitvec_ref[i] <- abs(lm(y_ref ~ distmat[i,])$coefficients[1] / 2)
    }
  }
  
  names(fitvec_ref) <- colnames(TOT_sel)[2:ncol(TOT_sel)]
  startCalchange <- startCal#"2016-01-01"
  write.table(fitvec_ref, paste0("data/obserrors/obserror_",scen,"_", startCalchange,"_",endCal,"_",variable,"_",date,".csv"), row.names=T, col.names =  F, sep = ",")
  
}
