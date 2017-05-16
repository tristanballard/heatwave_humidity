## This code just figures out what the threshold tmax value should be at each grid point.
## threshold.part2.R then takes that matrix of values and sees when the heat wave days are
## by seeing whden the values exceed the thresholds calculated here. Code here also includes 
## the same thing but creating a matrix of thresholds for humidity instead.


suppressMessages(library(fields))
suppressMessages(library(ncdf4))
suppressMessages(library(RColorBrewer))
suppressMessages(library(abind))

yr.start=1979
yr.end=2010
n.years=yr.end-yr.start+1


#############################################################
###########          Read in the Data         ###############
############################################################# 
#the data has already been subsetted for the particular region of interest
#station.data=readRDS('/scratch/users/tballard/shum/station.data/station.data.USOnly.rds') #23375x3x89
#dates=readRDS('/scratch/users/tballard/shum/station.data/station.data.USOnly.dates.rds') #23375x89
station.data=readRDS('/scratch/users/tballard/shum/station.data/station.data.SE.rds') #23375x3x89
dates=readRDS('/scratch/users/tballard/shum/station.data/station.data.SE.dates.rds') #23375x89



#### Subset for the years of interest ####
t.start=paste(yr.start,"-01-01",sep=""); t.start=as.Date(t.start) #e.g. "1979-01-01"
t.end=paste(yr.end, "-12-30", sep=""); t.end=as.Date(t.end) #"2010-12-31
yr.index=dates[,1]>=t.start & dates[,1]<=t.end
station.data2=station.data[yr.index,,]  #note there is some odd bug where the date for "1967-04-08" is NA 
station.data2=station.data2[-1,,] #for some reason, also an odd bug, a NA is at the start of the subsetted data, so remove it
## the bug is that the # of TRUE's in yr.index is 11324, but station.data2 has dim=11325. Very odd.

#### Figure out when leap years are #####
n.days=c(NA) #initialize
for (i in yr.start:yr.end){
  if (sum(dates==paste(i,"-02-29",sep=""),na.rm=T)>0) {
    n.days=c(n.days,366)
  }
  else {
    n.days=c(n.days,365)
  }
}
n.days=n.days[-1] #remove initializing NA; n.days=c(365,366,365,...), length=n.years



##### Compute the 90th percentile threshold for 7 days forward + back distribution #####
##### Note this is like a moving average where the first 7 days and last 7 days    #####
##### are cut. These two tails are calculated in the code directly after this and  #####
##### instead use years 1980-2014 for e.g. the first 7 days values.                #####

tmax=station.data2[,3,] #11687x89

time.start=8 #29
time.end= 365-7#1460-28
threshold.fxn=function(data,n.days,time.start,time.end){
  for (i in time.start:time.end){
  preceding.days=(i-7):(i-1)
  #preceding.days=preceding.days[c(TRUE,rep(FALSE,3))] #length=7
  following.days=(i+1):(i+7)
  #following.days=following.days[c(rep(FALSE,3),TRUE)] #length=7
  index=c(preceding.days,i,following.days) #length=15
  days.index=c()
  for (i in 1:(length(n.days)-1)){
    b=index+sum(n.days[1:i]) #sequence above, but for all the 1,2,3,4...36 years forward
    days.index=c(days.index,b)
  }
  days.index=c(index,days.index) #length = 540 = 15days*36years
  
  #Now subset the data by those days, and find the 90th percentile of that distribution
  x=data[days.index]
  x=as.numeric(x) #needed b/c tmax/station.data values are characters originally " "
  threshold=quantile(x, .9, names=F, na.rm=T)
  thresholds=c(thresholds,threshold) #vector length=1404
  #print(length(thresholds))
  }
  return(thresholds)
}
thresholds=c()
threshold=apply(tmax,c(2),threshold.fxn, n.days=n.days, time.start=time.start, time.end=time.end) #7x192x94
#saveRDS(thresh,"/scratch/users/tballard/shum/percentile.threshold/threshold")
#threshold=readRDS("/scratch/users/tballard/shum/percentile.threshold/threshold")


#### Compute for the first 7 days #####
time.start=366 #jan 1 1980
time.end=372 #jan 7 1980
threshold.fxn=function(data,n.days,time.start,time.end){
  for (i in time.start:time.end){
  preceding.days=(i-7):(i-1)
  #preceding.days=preceding.days[c(TRUE,rep(FALSE,3))] #length=7
  following.days=(i+1):(i+7)
  #following.days=following.days[c(rep(FALSE,3),TRUE)] #length=7
  index=c(preceding.days,i,following.days) #length=15
  days.index=c()
  for (i in 1:(length(n.days)-1)){
    b=index+sum(n.days[2:i]) #sequence above, but for all the 1,2,3,4...36 years forward
    #note this line above changes compared to the original one for all dates
    days.index=c(days.index,b)
  }
  days.index=c(index,days.index) #length = 540 = 15days*36years
  
  #Now subset the data by those hours, and find the 90th percentile of that distribution
  x=data[days.index]
  x=as.numeric(x) #needed b/c tmax/station.data values are characters originally " "
  threshold=quantile(x, .9, names=F, na.rm=T)
  thresholds=c(thresholds,threshold) #vector length=1404
  #print(length(thresholds))
  }
  return(thresholds)
}
thresholds=c()
threshold1=apply(tmax,c(2),threshold.fxn, n.days=n.days, time.start=time.start, time.end=time.end)
#saveRDS(thresh,"/scratch/users/tballard/shum/percentile.threshold/threshold.part1")
#threshold1=readRDS("/scratch/users/tballard/shum/percentile.threshold/threshold.part1")



##### Compute for the last 7 days of the year #####
time.start=365-6
time.end=365
threshold.fxn=function(data,n.days,time.start,time.end){
  for (i in time.start:time.end){
    preceding.days=(i-7):(i-1)
    #preceding.days=preceding.days[c(TRUE,rep(FALSE,3))] #length=7
    following.days=(i+1):(i+7)
    #following.days=following.days[c(rep(FALSE,3),TRUE)] #length=7
    index=c(preceding.days,i,following.days) #length=15
    days.index=c()
    for (i in 1:(length(n.days)-1)){
      b=index+sum(n.days[2:i]) #sequence above, but for all the 1,2,3,4...36 years forward
      #note this line above changes compared to the original one for all dates
      days.index=c(days.index,b)
    }
    days.index=c(index,days.index) #length = 540 = 15days*36years
    
    #Now subset the data by those hours, and find the 90th percentile of that distribution
    x=data[days.index]
    x=as.numeric(x)
    threshold=quantile(x, .9, names=F, na.rm=T)
    thresholds=c(thresholds,threshold) #vector length=1404
    #print(length(thresholds))
  }
  return(thresholds)
}
thresholds=c()
threshold2=apply(tmax,c(2),threshold.fxn, n.days=n.days, time.start=time.start, time.end=time.end)
#saveRDS(thresh2,"/scratch/users/tballard/shum/percentile.threshold/threshold.part2")
#threshold2=readRDS("/scratch/users/tballard/shum/percentile.threshold/threshold.part2")


##### Now combine the 3 threshold data frames into 1 universal one #####
# dim(threshold1)=7x96
# dim(threshold)=351x96
# dim(threshold2)=7x96

#### First rearrange the dimensions ####
threshold1=aperm(threshold1,c(2,1)); threshold=aperm(threshold,c(2,1)); threshold2=aperm(threshold2,c(2,1)); 
threshold4=abind(threshold1,threshold,threshold2) #dim=96x365
#saveRDS(threshold4, "/scratch/users/tballard/shum/station.data/threshold.final")
saveRDS(threshold4, "/scratch/users/tballard/shum/station.data/threshold.final.SE.rds")





