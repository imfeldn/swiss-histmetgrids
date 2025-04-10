#### calculate the running mean difference for modera between the past and the present ####
## this file calculate the running mean difference and writes it to a file that is later read in the pre-processing before calculating the analogues.
## ModERA can be downloaded here: https://www.wdc-climate.de/ui/project?acronym=ModE and plotted here http://climeapp-modera.unibe.ch:3838/
rm(list=ls())

require(PCICt)
require(lubridate)
require(ncdf4)

source('000_helpfun.R')

path <- "data/swiss_cell_ModE-RA_ensmean_temp2_abs_1420-2009.nc"
outdir <- "data"

refyear <- c(1978,2008)
refyear2 <- c(1971,2008)
histyear <- c(1763:2020)

aggmean <- function(x,time) {aggregate(x, list(month(time)), mean)$x}

print("read modERA")

modE_zonmean <- list()
thelp <- paste0(rep(1421:2009, each = 12),"-",pad2(1:12,2), "-15")
modE_zonmean$time <- as.Date(thelp[ncvar_get(nc_open(path), "time") + 1])
modE_zonmean$dat <- ncvar_get(nc_open(path),"temp2")
modera_yrs <- unique(year(modE_zonmean$time))

ref_year_list <- c(refyear,refyear2)

for(rfyrs in ref_year_list){
  
  yind_ref <- year(modE_zonmean$time) >= rfyrs[1]
  mon_ref <- aggmean(x=modE_zonmean$dat[yind_ref], time=modE_zonmean$time[yind_ref])
  yr_ref <- mean(modE_zonmean$dat[yind_ref])
  
  # calculate yearly offset
  offset_yr <- array(NA, dim=c(length(histyear)))
  yrspan <- 30
  yr2 <- 0
  
  for(yr in histyear){
    if(yr < rfyrs[1]){
      yr2 <- yr2 + 1
      yind_hist <- year(modE_zonmean$time) %in% c((yr-yrspan+1):(yr+yrspan))
      if (sum(yind_hist) == 720){
        yr_hist <- mean(modE_zonmean$dat[yind_hist])
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
  
  ## write offset to file as txt
  write.table(x= offset_yr,file=paste0("code/runclim_offset_year_modera_linintrp.txt"), row.names = T)
  
}
