# Calculate Analogs
# analogue resampling method (ARM)

### input needed:
# from preprocessing: predictor matrix, predictand matrix
# from masterscript:  specifications
### options: 
# select distance measure, number of analogues, seasonal window, 
# time window around target for leave-one-out validation

#############################################################################################################
#############################################################################################################
# 1. DEFINE DISTANCE FUNCTION
# 2. ANALOGUE CALCULATION
#    2.1 DEFINE TARGET AND ANALOGUE POOL
#    2.2 DEFINE ANALOGUE FUNCTIONS
#    2.3 CALCULATE ANALOGUE DATES AND DISTANCES
# 3. WRITE TABLES
#############################################################################################################
#############################################################################################################

# get time
ptm <- proc.time()

#############################################################################################################
### 1. DEFINE DISTANCE FUNCTIONS
#############################################################################################################
print(paste("define distance function:",specifications$distance))

distance <- function(predictor, target, type){
  
  if(type=="Euclid") {
    # distance function: euclidean distance
    # NA's are removed, if all NA values -> NA
    # distance is weighted with inverse number of parameters (penalty for less parameters available)
    # similar to RMSE which weights by square root of inverse number of parameters (weaker penalty)
    if (all(is.na(target))){
      dist <- NA
      }
    else {
      dist <- sqrt(sum((target - predictor) ^ 2, na.rm = T))/length(which(!is.na(predictor)))}
  } 
  if(type=="RMSE") {
       if (all(is.na(target))){
         dist <- NA
         }
      else {dist <- rmse(target,predictor, na.rm = specifications$narm)}
    } 
      ### Mahalanobis distance: not good, should be separated for individual groups of variables or stations
  if(type=="Mahalanobis") {
        cov.mat <- cov(data[,-c(1,2)], use = "pairwise.complete.obs")
          if (all(is.na(target))){
            dist <- NA
            }
          else {
            dist <- sqrt(mahalanobis(target,predictor,cov=cov.mat))/length(which(!is.na(predictor)))}
      } 
  ### MAE
  if(type=="MAE") {
            if (all(is.na(target))){
              dist <- NA
              }
            else {
              dist <- mean(abs(target-predictor))/length(which(!is.na(predictor)))}
  
        } 
  ### Gower  
  if(type=="Gower") {dist <- NA}
  return(dist)
}

#############################################################################################################
### 2. ANALOGUE CALCULATION
#############################################################################################################

#------------------------------------------------------------------------------------------------------------
### 2.1 DEFINE TARGET AND ANALOGUE POOL
#------------------------------------------------------------------------------------------------------------
Analogpool <- predictor
rownames(Analogpool) <- Analogpool$date
Target <- predictand
rownames(Target) <- Target$date
Target <- split(Target,Target$date)

#------------------------------------------------------------------------------------------------------------
### 2.2 DEFINE ANALOGUE FUNCTIONS
#------------------------------------------------------------------------------------------------------------
print("define analogue function")

## outcome is a list of days with "nAnalog" best analogue dates and corresponding distances
Analogs = function(target,nAnalog,dayselection,window,rr0crit) {
  
  # define date to reconstruct and vector of weather types
  day <- target[["date"]]
  if(is.na(wt_probs)){wt_probs = 7}
  
  if(all(is.na(target[16:ncol(target)]))){
    analogue.ls <- switch((rr0crit!="")+1,
                          list("analog.dates"=rep(NA,nAnalog),"analog.dist"=rep(NA,nAnalog)),
                          list("analog.dates"=rep(NA,nAnalog),"analog.dist"=rep(NA,nAnalog),"analog.rscore"=rep(NA,nAnalog)))
  } else {
    
    w_type <- target$WT_type1
    
    # create matrix with measurements of target day
    target_data <- target[names(target)[!grepl("WT",names(target)) & !grepl("date",names(target))]]
    
    # make preselection to select only days +-XXX days from target date
    preselection <- sapply(seq(year(specifications$startRef), year(specifications$endRef)), 
                           function(i) seq(as.Date(yday(day)-dayselection, origin=paste(i,"01-01",sep="-")),
                                           as.Date(yday(day)+dayselection, origin=paste(i,"01-01",sep="-")),
                                           by="day"))
    preselection <- c(as.Date(preselection,origin="1970-01-01"))
    
    ## if predictand within predictor period: remove predictand dates +-window from predictor
    pres_window <- seq(as.Date(day)-window, as.Date(day)+window, by="day")
    preselection <- preselection[!preselection%in%pres_window]
    
    preselection <- preselection[preselection >= as.Date(specifications$startRef)]
    preselection <- preselection[preselection <= as.Date(specifications$endRef)]
    
    ## select only days with fullfilled rrt
    if(rr0crit!=""){
      rr0_target <-  as.matrix(target_data[grepl("rr0",names(target_data))])
      targlen <- length(rr0_target[!is.na(rr0_target)])
      rr0_analog <- Analogpool[,grepl("rr0",colnames(Analogpool))]
      rr0score <- round(apply(sweep(as.matrix(rr0_analog), MARGIN=2, STATS=rr0_target),1,function(x){sum(x^2)/targlen}),2)
      rr0sel <- rr0score <= rr0crit
      rscore <- rr0score[Analogpool$WT_type1%in%w_type & Analogpool$date%in%preselection]
    }
    
    # select possible analogues according to weather types and "preselection"
    fullselection <- switch((rr0crit!="")+1,
                            Analogpool$WT_type1%in%w_type & Analogpool$date%in%preselection,
                            Analogpool$WT_type1%in%w_type & Analogpool$date%in%preselection & rr0sel)
    analog_selection <- Analogpool[fullselection,]
    
    # calculate distances and get date/distance matrices
    if(specifications$distance=="Gower"){
      dist <- gower.dist(data.x=target_data, data.y=analog_selection[, -c(1:2)])[1,]
      names(dist) <- analog_selection$date
    } else {
      dist<-apply(analog_selection[,-c(1:2)],1,distance,target=as.numeric(target_data), type=specifications$distance)
    }
    ord <- order(dist)
    dist.min <- dist[ord][1:nAnalog]
    dist.index <- as.Date(names(dist.min))
    dist <- as.numeric(dist.min)
    # create list with analog dates and distances for target day
    # the last value is only passed for evaluating the amount of fewer analogs
    analogue.ls <- switch((rr0crit!="")+1,
                          list("analog.dates"=as.numeric(dist.index),"analog.dist"=dist),
                          list("analog.dates"=as.numeric(dist.index),"analog.dist"=dist,"analog.rscore"=as.numeric(rscore)))
    
    
  }
  return(analogue.ls)
}

#------------------------------------------------------------------------------------------------------------
### 2.3 CALCULATE ANALOGUE DATES AND DISTANCES
#------------------------------------------------------------------------------------------------------------
print(paste("calculate best",specifications$nAnalog,"analogues"))

### create list of days, containing  analogue dates and distances
#Rprof("test5.out")
analogue.list<-mcmapply(Analogs,target=Target,nAnalog=specifications$nAnalog,
                        dayselection=specifications$dayselection,window=specifications$window,rr0crit=specifications$rr0crit, SIMPLIFY = F,
                        mc.cores = round(0.15*detectCores()))
#Rprof(NULL)

### extract analogue dates from list as matrix
analog.dates<-data.frame(as.Date(names(analogue.list)),t(modify_depth(analogue.list,1,~ .[["analog.dates"]]) %>% as.data.frame))
colnames(analog.dates)<-c("date",1:specifications$nAnalog)
rownames(analog.dates)<-c(1:nrow(analog.dates))

### extract analogue distances from list as matrix
analog.dist<-data.frame(as.Date(names(analogue.list)),t(modify_depth(analogue.list,1,~ .[["analog.dist"]]) %>% as.data.frame))
colnames(analog.dist)<-c("date",1:specifications$nAnalog)
rownames(analog.dist)<-c(1:nrow(analog.dist))

### rscore of preselection - keep list format
if(specifications$rr0crit!=""){
  analog.rscore<- modify_depth(analogue.list,1,~ .[["analog.rscore"]])
}

#############################################################################################################
### 3. WRITE TABLES
#############################################################################################################
print("write files to folder")

## get input for names of files
if (specifications$seasonality==T) {props1 <- "deseas_"} else {props1 <- c()}
if (specifications$trend!=F) {props2 <- "detrend_"} else {props2 <- c()}
props3 <- specifications$distance
if(specifications$scenario!=""){
  props4 <- specifications$scenario
  props4a <- paste0("_",specifications$dayselection)
  props5 <- ""
  
} else {
  props4 <- paste0(specifications$variables,collapse="_")
  props5 <- gsub(",","_",specifications$stations)
}
props6 <- specifications$rr0crit

## create output folder for validation
## use separate folder, if I do the station validation
if(length(stnval)>0 & stnval!="all"){
  
  subDir <- paste("ARM_run",date,sep="_")
  dir.create(file.path("output",subDir), showWarnings = F)
  
  subDir2 <- paste(props4,specifications$dayselection,stnval,sep="_")
  dir.create(file.path("output",subDir,subDir2), showWarnings = F)
  
  outfolder <- paste(file.path("output",subDir,subDir2),"/",sep="")
  rm(subDir,subDir2)
  
} else {
  
  subDir <- paste("ARM_run",date,sep="_")
  dir.create(file.path("output",subDir), showWarnings = F)
  
  subDir2 <- paste(props4,specifications$dayselection,sep="_")
  dir.create(file.path("output",subDir,subDir2), showWarnings = F)
  
  outfolder <- paste(file.path("output",subDir,subDir2),"/",sep="")
  rm(subDir,subDir2)
  
}

## save tables (dates, distances, readme)
sdate <- as.character(analog.dates[1,1])
edate <- as.character(analog.dates[nrow(analog.dates),1])

write.table(analog.dates, file=paste(outfolder,date,"_analog.dates_",sdate,"-",edate,"_",props1,props2,props3,"_",props4,props4a,props5,props6,".txt",sep=""))
write.table(analog.dist,file=paste(outfolder,date,"_analog.dist_",sdate,"-",edate,"_",props1,props2,props3,"_",props4,props4a,props5,props6,".txt",sep=""))
if(specifications$rr0crit!=""){
  save(analog.rscore,file=paste(outfolder,date,"_analog.rscore_",sdate,"-",edate,"_",props1,props2,props3,"_",props4,props4a,props5,props6,".RData",sep=""))
}
readme<-data.frame("description"=c("analog calculation",paste(names(specifications),":",sep="")),"parameters"=NA)
readme[1,2]<-as.character(as.Date(Sys.Date()))
for(i in 2:nrow(readme)) {readme[i,2]<-as.character((paste(specifications[[i-1]], collapse = ", ")))}
readme[2,2]<- paste(unique(substr(colnames(predictand),1,3))[-c(1,2)],collapse = ", ")
write.table(readme,file=paste(outfolder,date,"_readme_analogs_",sdate,"-",edate,"_",props1,props2,props3,"_",props4,props4a,props5,props6,".txt",sep=""))

#############################################################################################################
## return summary
print(paste("analogue reconstruction from",sdate,"to",edate,sep=" "))
etm <- proc.time()
print(paste("calculation time:",(etm[3]-ptm[3])%/%60,"minutes", round((etm[3]-ptm[3])%%60),"seconds" ,sep=" "))
print(paste("files are in:",outfolder,sep=" "))

## remove variables
rm(Analogpool, Target, props1, props2, props3, props4, props5, props6, sdate, edate, ptm, etm, i, 
   outfolder, Analogs, distance, readme, analogue.list)
if(specifications$distance=="Mahalanobis"){rm(cov.mat)}

# END