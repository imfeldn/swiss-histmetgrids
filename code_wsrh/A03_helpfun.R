### help fun genereal

require(lubridate)

leapyear <- function(x){
  return(ifelse((x%%4==0 & x%%100!=0) | x%%400==0,TRUE,FALSE))
}

pad2 <- function(x, width=2){
  stringr::str_pad(x,width,pad="0")
}


cllps <- function(x,pat=""){
  if(is.matrix(x)){
    sapply(1:nrow(x), function(z) paste0(x[z,],collapse = pat))
  }else{
    paste0(x,collapse = pat)
  }
}

t_col <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  invisible(t.col)
}


pdate <- function(dat){as.Date(paste0(dat[,1],"-",pad2(dat[,2]),"-",pad2(dat[,3])))}

load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}


mean_Xpercna <- function(x, na=0.25){
  if(sum(is.na(x))/length(x) <= na){
    return(mean(x, na.rm=T))
  } else{return(NA)}}

sum_Xpercna <- function(x, na=0.25){
  if(sum(is.na(x))/length(x) <= na){
    return(sum(x, na.rm=T))
  } else{return(NA)}}



yday_without_29feb <- function(predperiod){
  leapyears <- leapyear(lubridate::year(predperiod))  
  jind <- which(leapyears & substr(predperiod,6,10)>"02-28")
  jdays <- lubridate::yday(predperiod)
  if (length(jind) > 0) jdays[jind] <- jdays[jind]-1 
  return(jdays)
}



quant_window_fun <- function(dat,dates,doy){
  window <-8
  doysel <- (doy-window):(doy+window)
  doysel[which(doysel<1)] <- doysel[which(doysel<1)]+365
  doysel[which(doysel>365)] <- doysel[which(doysel>365)]-365
  doyind <- which(yday_without_29feb(dates) %in% doysel)
  q <- quantile(dat[doyind],probs=c(0.1,0.25,0.5,0.75,0.9), na.rm = T)
  return(q)
}




rm_seas <- function(tempdata,dates) {
  
  if(!all(is.na(tempdata))){
    
    doys <- yday(dates)
    ndoy <- lubridate::leap_year(year(dates))+365
    ndoy <- ndoy[1:length(dates)]
    xdoys <-doys/ndoy

    pcos <- cos(2*pi*xdoys)
    psin <- sin(2*pi*xdoys)
    pcos2 <- cos(4*pi*xdoys)
    psin2 <- sin(4*pi*xdoys)
    
    temp.jg.lm <- lm(tempdata ~ pcos + psin + pcos2 + psin2)
    temp.jg.cor <- as.numeric(tempdata[!is.na(tempdata)] - temp.jg.lm$fitted.values)
    tempdata[!is.na(tempdata)] <- temp.jg.cor
  }
  
  return(tempdata)
}   



rm_seas <- function(tempdata,doys, yrs) {
  
  if(!all(is.na(tempdata))){
    
    ndoy <- lubridate::leap_year(yrs)+365
    ndoy <- ndoy[1:length(doys)]
    xdoys <-doys/ndoy
    
    pcos <- cos(2*pi*xdoys)
    psin <- sin(2*pi*xdoys)
    pcos2 <- cos(4*pi*xdoys)
    psin2 <- sin(4*pi*xdoys)
    
    temp.jg.lm <- lm(tempdata ~ pcos + psin + pcos2 + psin2)
    temp.jg.cor <- as.numeric(tempdata[!is.na(tempdata)] - temp.jg.lm$fitted.values)
    tempdata[!is.na(tempdata)] <- temp.jg.cor
  }
  
  return(tempdata)
}   

rm_seas_sep <- function(tempdata,pcos,psin,pcos2,psin2) {
  
  if(!all(is.na(tempdata))){
    
    # doys <- yday(dates)
    # ndoy <- leapyear(year(dates))+365
    # ndoy <- ndoy[1:length(dates)]
    # xdoys <-doys/ndoy
    # 
    # pcos <- cos(2*pi*xdoys)
    # psin <- sin(2*pi*xdoys)
    # pcos2 <- cos(4*pi*xdoys)
    # psin2 <- sin(4*pi*xdoys)
    
    temp.jg.lm <- lm(tempdata ~ pcos + psin + pcos2 + psin2)
    temp.jg.cor <- as.numeric(tempdata[as.integer(names(temp.jg.lm$fitted.values))] - temp.jg.lm$fitted.values)
    tempdata[!is.na(tempdata)] <- temp.jg.cor
  }
  
  return(tempdata)
}   



rm_seas2 <- function(tempdata,dates) {
 
  if(!all(is.na(tempdata))){
    
    doys <- lubridate::yday(dates)
    ndoy <- leapyear(lubridate::year(dates))+365
    ndoy <- ndoy[1:length(dates)]
    xdoys <-doys/ndoy
    
    pcos <- cos(2*pi*xdoys)
    psin <- sin(2*pi*xdoys)
    pcos2 <- cos(4*pi*xdoys)
    psin2 <- sin(4*pi*xdoys)
    
    temp.jg.lm <- lm(tempdata ~ pcos + psin + pcos2 + psin2)
    temp.jg.cor <- as.numeric(tempdata[names(temp.jg.lm$fitted.values)] - temp.jg.lm$fitted.values)
    tempdata[!is.na(tempdata)] <- temp.jg.cor
  }
  
  return(tempdata)
}   

# get date factor from time
get_date_factors_obj<-function(t1d){
  factors<-list()
  factors$years<-factor(format(t1d,"%Y"))
  factors$yearmons<-factor(format(t1d,"%Y-%m"))
  jdays<-format(t1d,"%j")
  jdays<-replace_jday_29feb(jdays)
  factors$jdays=factor(jdays)
  return(factors)
}


replace_jday_29feb <- function(jdays){
  indices <- which(jdays == 366)
  if (length(indices) > 0)
    jdays[rep(indices, each = 366) + -365:0] <- pad2(c(1:59,59,60:365), width=3)
  return(jdays)
}


trend_linreg <-
  function (dat.df,conf.probs=c(0.025,0.975)) {
    # predictor
    year <- as.array(dat.df[,1])
    # predictand (responses)
    resp <- as.array(dat.df[,2])
    # calculate linear regression
    trd.lm <- stats::lm(resp ~ year,na.action=na.omit)
    # prepare elements of trend object:
    # intermediate calculations
    coef <- stats::coefficients(trd.lm)   # model coefficients
    # standard error of slope
    stderr <- summary(trd.lm)$coefficients[2,2]
    t2 <- year[length(year)]      # last year
    t1 <- year[1]                 # first year
    t0 <- (t2+t1)/2               # middle year
    v0 <- coef[1]+t0*coef[2]      # value at middle year
    # fitted values:
    indx <- !is.na(resp)       # fitted omits NAs in predictor !
    fits <- matrix(c(year[indx],fitted(trd.lm)),ncol=2)
    
    # trend magnitude as absolute change over the full period
    trd.abs <- ((coef[1]+t2*coef[2])-(coef[1]+t1*coef[2]))
    # trend magnitude as fractional change of average
    if (v0 == 0) {
      trd.rel <- NA
    } else {
      trd.rel <- ((coef[1]+t2*coef[2])-(coef[1]+t1*coef[2]))/v0
    }
    # trend magnitude as ratio between end/beginning
    if ((coef[1]+t1*coef[2]) == 0) {
      trd.rat <- NA
    } else {
      trd.rat <- ((coef[1]+t2*coef[2])/(coef[1]+t1*coef[2]))
    }
    
    # confidence interval of trend magnitude
    df <- summary(trd.lm)$df[2]
    tcrit <- stats::qt(p=conf.probs,df=df)
    trd.conf.abs <- ((t2-t1)*(coef[2]+stderr*tcrit))
    names(trd.conf.abs) <- as.character(conf.probs)
    #trd.conf.rel <- ((v0+(t2-t0)*(coef[2]+stderr*tcrit))-
    #                 (v0+(t1-t0)*(coef[2]+stderr*tcrit)))/v0
    #names(trd.conf.rel) <- as.character(conf.probs)
    
    # p-value of trend coefficient using standard T-Test
    pval <- summary(trd.lm)$coefficients[2,4]
    
    # return trend-object
    res <- list(data=dat.df,
                mod=trd.lm,
                trend.abs=trd.abs,
                trend.rel=trd.rel,
                trend.rat=trd.rat,
                trend.conf.abs=trd.conf.abs,
                #trend.conf.rel=trd.conf.rel,
                coef=coef,
                pval=pval,
                fitted=fits,
                method="linreg")
    return(res)
  }