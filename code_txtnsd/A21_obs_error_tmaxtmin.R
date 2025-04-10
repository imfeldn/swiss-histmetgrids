### Calculate error covariance for observations
# Two files are created based on different periods defined in startCal and endCal.

rm(list = ls())

library(geosphere)
library(lubridate)

# Settings
anomalies <- TRUE             # Use anomalies or absolute values
minlength <- 120              # Minimum overlapping period in months
startCal <- "1763-01-01"      # Calibration period start
endCal <- "1863-12-31"        # Calibration period end
detr <- TRUE
scen <- "all_hist"
date <- ""

## Function to calculate anomalies (monthly climatology based)
anom <- function(x, nrmlizd = FALSE) {
  xclim <- aggregate(x[, 3], list(x[, 1]), mean, na.rm = TRUE)
  stdclim <- aggregate(x[, 3], list(x[, 1]), sd, na.rm = TRUE)
  for (i in 1:12) {
    if (nrmlizd) {
      x[x[, 1] == i, 3] <- (x[x[, 1] == i, 3] - xclim$x[xclim[, 1] == i]) / stdclim$x[stdclim[, 1] == i]
    } else {
      x[x[, 1] == i, 3] <- x[x[, 1] == i, 3] - xclim$x[xclim[, 1] == i]
    }
  }
  return(x)
}

# Read station coordinates
station_coords <- XLConnect::readWorksheetFromFile("data/station_inventory.xlsx", sheet = 1, startRow = 2)[, c("ID", "lon_ref", "lat_ref")]
station_coords$ID <- ifelse(grepl("tx", station_coords$ID),
                            gsub("tx", "tmax", station_coords$ID),
                            gsub("tn", "tmin", station_coords$ID))

# Load station data
load(paste0("data/analogue_stationdata_offset_dtrd1971_", scen, "_txtn_", date, ".RData"))
TOT <- TOT[, !grepl("WT_", colnames(TOT))]

# Variables of interest
var_th <- list("ta" = 10, "tmax" = 16, "tmin" = 15)
vars <- c("ta")

# Main loop for each variable
for (variable in vars) {
  selection <- paste0("_", variable)
  TOT_sel <- TOT[, c(1, grep(selection, colnames(TOT)))]
  n <- ncol(TOT_sel) - 1
  
  # Read data and metadata into a list
  Data_ref <- list()
  
  for (i in 1:n) {
    Data_ref[[i]] <- list()
    id <- colnames(TOT_sel)[i + 1]
    print(id)
    
    idsel <- station_coords$ID == id
    Data_ref[[i]]$coords <- as.numeric(station_coords[idsel, c("lon_ref", "lat_ref")])
    
    yind <- which(year(TOT_sel$date) <= year(endCal) & year(TOT_sel$date) >= year(startCal))
    x <- TOT_sel[yind, c("date", id)]
    Data_ref[[i]]$dates <- ISOdate(year(x$date), month(x$date), 1)
    
    narows <- sum(!is.na(x[, 2]))
    if (anomalies & narows >= minlength) {
      x <- cbind(month(x$date), year(x$date), x[, 2])
      Data_ref[[i]]$ta <- anom(x, nrmlizd = FALSE)[, 3]
    } else {
      Data_ref[[i]]$ta <- x$Value
    }
  }
  
  names(Data_ref) <- colnames(TOT_sel)[2:ncol(TOT_sel)]
  
  # Calculate distance and variance of pairwise differences
  cormat_ref <- distmat <- array(dim = c(n, n))
  colnames(cormat_ref) <- rownames(cormat_ref) <- colnames(TOT_sel)[2:ncol(TOT_sel)]
  fitvec_ref <- array(dim = n)
  
  for (i in 1:n) {
    message(i)
    for (j in 1:n) {
      if (i != j) {
        distmat[i, j] <- (distGeo(Data_ref[[i]]$coords, Data_ref[[j]]$coords) / 1000)^2
        if (is.na(distmat[i, j])) distmat[i, j] <- 1e12
        
        x_ref <- Data_ref[[i]]$ta[which(Data_ref[[i]]$dates %in% Data_ref[[j]]$dates)]
        y_ref <- Data_ref[[j]]$ta[which(Data_ref[[j]]$dates %in% Data_ref[[i]]$dates)]
        
        if (sum(!is.na(x_ref - y_ref)) >= minlength & distmat[i, j] < 1000^2) {
          if (sd(x_ref, na.rm = TRUE) > 0 & sd(y_ref, na.rm = TRUE) > 0) {
            cormat_ref[i, j] <- var(x_ref - y_ref, use = "pairwise.complete.obs")
          }
        }
      }
    }
    
    # Linear fit (following Wartenburger et al.)
    y_ref <- cormat_ref[i, ]
    y_ref[y_ref > var_th[[variable]]] <- NA
    print(y_ref)
    
    if (sum(!is.na(y_ref)) >= 5) {
      fitvec_ref[i] <- abs(lm(y_ref ~ distmat[i, ])$coefficients[1] / 2)
    }
  }
  
  names(fitvec_ref) <- colnames(TOT_sel)[2:ncol(TOT_sel)]
  startCalchange <- startCal
  write.table(
    fitvec_ref,
    paste0("data/obserrors/obserror_txtn_", scen, "_", startCalchange, "_", endCal, "_", variable, "_", date, ".csv"),
    row.names = TRUE,
    col.names = FALSE,
    sep = ","
  )
}
