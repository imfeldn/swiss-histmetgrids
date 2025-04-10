#### check ou ekf difference during 1750 to 2003

rm(list=ls())

require(PCICt)
require(lubridate)
require(ncdf4)

main <- "/home/nimfeld/wear/swiss_grids/"
setwd(main)

source('/scratch3/nimfeld/wear/swiss_arm/code/data_prep/scripts/000_helpfun.R')
source("/scratch3/nimfeld/wear/summer1940s/scripts/R/get_ncdf.R")

ekfpath <- "/scratch3/nimfeld/wear/swiss_arm/data/mode-ra/swiss_cell_ModE-RA_ensmean_temp2_abs_1420-2009.nc"
outdir <- "data/tabs/"

refyear <- c(1978,2008)
refyear2 <- c(1971,2008)
histyear <- c(1763:2020)

aggmean <- function(x,time) {aggregate(x, list(month(time)), mean)$x}

print("read ekf")
# ekf_zonmean <- list()
# ekf_zonmean$time <- as.Date(as.character(as.Date(ncvar_get(nc_open(ekfpath),"time")*86400, origin="1600-01-01")))
# ekf_zonmean$dat <- ncvar_get(nc_open(ekfpath),"air_temperature")
# ekf_zonmean$lat <- ncvar_get(nc_open(ekfpath),"lat")

ekf_zonmean <- list()
thelp <- paste0(rep(1421:2009, each = 12),"-",pad2(1:12,2), "-15")
ekf_zonmean$time <- as.Date(thelp[ncvar_get(nc_open(ekfpath), "time") + 1])
ekf_zonmean$dat <- ncvar_get(nc_open(ekfpath),"temp2")
modera_yrs <- unique(year(ekf_zonmean$time))

ref_year_list <- c(refyear,refyear2)

for(rfyrs in ref_year_list){
  
  yind_ref <- year(ekf_zonmean$time) >= rfyrs[1]
  mon_ref <- aggmean(x=ekf_zonmean$dat[yind_ref], time=ekf_zonmean$time[yind_ref])
  yr_ref <- mean(ekf_zonmean$dat[yind_ref])
  
  # calculate yearly offset
  offset_yr <- array(NA, dim=c(length(histyear)))
  yrspan <- 30
  yr2 <- 0
  
  for(yr in histyear){
    if(yr < rfyrs[1]){
      yr2 <- yr2 + 1
      yind_hist <- year(ekf_zonmean$time) %in% c((yr-yrspan+1):(yr+yrspan))
      if (sum(yind_hist) == 720){
        yr_hist <- mean(ekf_zonmean$dat[yind_hist])
        offset_yr[yr2] <- yr_ref - yr_hist
      }  
    }
  }
  
  names(offset_yr) = histyear
  write.table(x= offset_yr,file=paste0("code/runclim_offset_year_modera_till",rfyrs[1],".txt"), row.names = T)
  
  ## linear interpolate year till 2016-2020 mean
  aa <- data.frame(yr = c(1978,2020), off = c(offset_yr["1978"],0))
  mod <- lm(off ~yr, data = aa)
  interp <- predict(mod,data.frame( yr = 1979:2020))
  offset_yr[which(names(offset_yr) %in% as.character(1979:2020))] <- interp
  
  plot(histyear, offset_yr)
  
  #brks <- seq(0,1.4,0.1)
  #redblues <-  colorRampPalette(c("navyblue","cornflowerblue","yellow","red","darkred"))
  #cols <- redblues(length(brks)-1)
  
  #png("code/varia/offset/offset_yearly_bylatitude_ekf400.png", width=300, height = 1300, res = 120)
  # png("../thesis/defence/plots/climate_offset.png", width=300, height = 800, res = 120)
  # par(mar=c(3,3,1,2))
  # fields::image.plot(histyear,1, as.matrix(offset_yr), breaks=brks, col=cols, xaxt = "n", xlab = "", ylab = "", horizontal = T)
  # axis(side = 1, at = 1:5, labels = round(ekf_zonmean$lat,1), las = 2)
  # #title("years offset from EKF400 by latitude")
  # dev.off()
  
  ## write offset to file as txt
  write.table(x= offset_yr,file=paste0("code/runclim_offset_year_modera_linintrp.txt"), row.names = T)
  
}

########################################
## calculate a monthly offset
########################################
# 
# offset <- array(NA, dim=c(length(mon_ref),length(histyear)))
# yrspan <- 35
# yr2 <- 0
# 
# for(yr in histyear){
#   yr2 <- yr2 +1
#   yind_hist <- year(ekf_zonmean$time) %in% c((yr-yrspan+1):(yr+yrspan))
#   mon_hist <- aggmean(x=ekf_zonmean$dat[yind_hist], time=ekf_zonmean$time[yind_hist])
#   offset[,yr2] <- mon_ref - mon_hist
# }
# 
# offset_latavg <- apply(offset,c(1,3), mean)
# par(mar=c(3,3,2,2))
# brks <- seq(-2.6,2.6,0.1)
# redblues <-  colorRampPalette(c("navyblue","cornflowerblue","yellow","red","darkred"))
# cols <- redblues(length(brks)-1)
# png("code/varia/offset/offset_monthly_ekf400.png", width=800)
# fields::image.plot(histyear,1:12, t(offset), col=cols, breaks=brks)
# title("monthly offset from EKF400")
# dev.off()
# 
# png("code/varia/offset/zonalmean_ekf400.png", width=1300, height = 300, res = 120)
# par(mar=c(3,3.5,2,2))
# #plot(1730:2003, caTools::runmean(aggregate(ekf_zonmean$dat[3,], list(year(ekf_zonmean$time)), mean)$x, k = 10), type = "l")
# plot(1730:2003, aggregate(ekf_zonmean$dat[3,], list(year(ekf_zonmean$time)), mean)$x, type = "l", 
#      xlab = "", ylab = "", lwd = 2, col = "navyblue")
# arrows(x0= 1961, x1 = 2003, y0 = 3.5, y1 = 3.5, lwd =2, code =3, length = 0.05)
# mtext(side = 3, "EKF400 zonal mean 47.6°N", line = 0.1)
# dev.off()
# 
# 
# ## plot individual month time series
# redblues <-  colorRampPalette(c("navyblue","cornflowerblue","yellow","red","darkred","yellow","navyblue","cornflowerblue"))
# mocols <-  redblues(length(1:12))
# 
# png("code/varia/offset/zonmean_timeseries_monthly_ekf400.png", width = 1000)
# moind <- which(month(ekf_zonmean$time) == 1)
# plot(ekf_zonmean$time[moind], ekf_zonmean$dat[1,moind], type="l", ylim=c(-12,25), col=mocols[1], ylab="degrees celsius",xlab="time")
# text(x= ekf_zonmean$time[moind][1], ekf_zonmean$dat[1,moind][1],labels = 1)
# 
# for(m in 2:12){
#   moind <- which(month(ekf_zonmean$time) == m)
#   points(ekf_zonmean$time[moind], ekf_zonmean$dat[1,moind], type="l", col=mocols[m])
#   text(x= ekf_zonmean$time[moind][1], ekf_zonmean$dat[1,moind][1],labels = m)
# }
# title("zonal mean timeseries for every month 1-12 individually EKF400")
# dev.off()
# 
# ################################################
# ## calculate a monthly offset based on 20crv3
# ################################################
# 
# dat20cr  <- read.table("data/20crv3/20crv3_monthly.csv", sep = ",", colClasses = "character")
# dat20cr$V1 <- as.Date(as.character(dat20cr$V1), format="%Y-%m-%d")
# dat20cr$V2 <- as.numeric(dat20cr$V2)
# 
# yind_ref <- year(dat20cr$V1) %in% 1961:2003
# momean_ref <- aggmean(dat20cr$V2[yind_ref], time=dat20cr$V1[yind_ref]) -273.15
# 
# histyear2 <- 1836:1960
# 
# offset <- array(NA, dim=c(12,length(histyear2)))
# yrspan <- 15
# yr2 <- 0
# 
# for(yr in histyear2){
#   yr2 <- yr2 +1
#   yind_hist <- year(dat20cr$V1) %in% c((yr-yrspan+1):(yr+yrspan))
#   momean_hist <- aggmean(dat20cr$V2[yind_hist], time=dat20cr$V1[yind_hist]) - 273.15
#   offset[,yr2] <- momean_ref - momean_hist
# }
# 
# aggregate(dat20cr$V2[yind_hist], list(month(dat20cr$V1[yind_hist])), mean) - 273.15
# 
# 
# par(mar=c(3,3,2,2))
# brks <- seq(-2,2,0.2)
# redblues <-  colorRampPalette(c("navyblue","cornflowerblue","yellow","red","darkred"))
# cols <- redblues(length(brks)-1)
# png("code/varia/offset/offset_monthly_ekf400.png", width=800)
# fields::image.plot(histyear2,1:12, t(offset), col=cols, breaks=brks)
# title("monthly offset from 20crv3")
# dev.off()
# 
# ## plot individual month time series
# redblues <-  colorRampPalette(c("navyblue","cornflowerblue","yellow","red","darkred","yellow","navyblue","cornflowerblue"))
# mocols <-  redblues(length(1:12))
# 
# png("code/varia/offset/zonmean_timeseries_monthly_ekf400.png", width = 1000)
# moind <- which(month(ekf_zonmean$time) == 1)
# plot(ekf_zonmean$time[moind], ekf_zonmean$dat[1,moind], type="l", ylim=c(-12,25), col=mocols[1], ylab="degrees celsius",xlab="time")
# text(x= ekf_zonmean$time[moind][1], ekf_zonmean$dat[1,moind][1],labels = 1)
# 
# for(m in 2:12){
#   moind <- which(month(ekf_zonmean$time) == m)
#   points(ekf_zonmean$time[moind], ekf_zonmean$dat[1,moind], type="l", col=mocols[m])
#   text(x= ekf_zonmean$time[moind][1], ekf_zonmean$dat[1,moind][1],labels = m)
# }
# title("zonal mean timeseries for every month 1-12 individually EKF400")
# dev.off()
