suppressMessages(library(fields))
suppressMessages(library(ncdf4))
suppressMessages(library(RColorBrewer))
suppressMessages(library(abind))
suppressMessages(library(base))

#This takes the raw station observation data and creates a station.data.USOnly.rds file
#of the cleaned values where it is dimension = ndays x 4 variables (date, temp, shum, tmax) x n.stations
#The stations are subsetted at the start by picking the 'zone' or region of U.S. you want to look at.
#The data is cleaned up a bit, with all variables here removed except temp and shum. This
#also converts the hourly data into daily averages.

zone=c(30,44,255,290) #SE only
#zone=c(19,50,210,290) #whole US

meta=read.csv("/scratch/PI/omramom/CDMP/station_metadata.csv")
  #convert lon from -180,180 to 0,360
  lon = ifelse(meta$lon >= 0 , meta$lon, meta$lon + 360) #since all are negative you could also just go ahead and add 360 to all but w/e
  lat = meta$lat
  index=lat>=zone[1] & lat<=zone[2] & lon>=zone[3] & lon<=zone[4]
  #saveRDS(lon[index],'/scratch/users/tballard/shum/station.data/lon.rds')
  #saveRDS(lat[index],'/scratch/users/tballard/shum/station.data/lat.rds')
  saveRDS(lon[index],'/scratch/users/tballard/shum/station.data/lon.SE.rds')
  saveRDS(lat[index],'/scratch/users/tballard/shum/station.data/lat.SE.rds')
  
  
  t.start=max(meta$start[index]) #find a starting year common to all the stations w/in the zone
  t.start=paste(t.start,"-01-01",sep=""); t.start=as.Date(t.start) #e.g. "1947-01-01"
  t.end=min(meta$end[index]) #end time (=2010)
  t.end=paste(t.end, "-12-30", sep=""); t.end=as.Date(t.end) #"2010-12-30
  n.hours=(t.end-t.start+1)*24 #number of hours (nrows) you'll end up subsetting 
  n.days=n.hours/24 #23012 days
  
fileNames = list.files(path="/scratch/PI/omramom/CDMP/data/",pattern="*.dat")
fileNames = paste("/scratch/PI/omramom/CDMP/data/", fileNames, sep="")
fileNames.zone = fileNames[index] #filenames for the particular grid box zone
#dat=read.delim("/scratch/PI/omramom/CDMP/data/KABE_H.dat", head=F)

NA.matrix=matrix(rep(NA,n.days*5),nrow=n.days) #used to initialize, 5 is for the 5 variables ultimately stored
station.data=array(NA.matrix,c(n.days,5,1)) #n.daysx4x1 array of NA's to initialize
  for (i in 1:length(fileNames.zone)){
    a=read.delim(fileNames.zone[i], head=F) #read in the .dat file
    a=a[,-(7:9)] #remove last 3 columns b/c not relevant
    colnames(a)=c("date","hour","T","Td","rhum","shum")
    a=a[,-(4:5)] #remove Tdew point and rhum (optional); change the '5' above though if remove this line
    a[,1]=as.Date(as.character(a[,1])) #convert 1st column to 'date' R format
    years.index=a[,1]>=t.start & a[,1]<=t.end
    a=a[years.index,] #subset only the time period common to all 
    a[a==-999]=NA #set missing values of -999 to NA instead
    days=rep(1:n.days, each=24) #vector of 1,1,1,...2,2,2,...23012,23012... (23012 days in dataset here)
    aa=aggregate(a, by=list(days), mean, na.rm=T) #average every 24 hourly values into daily values
    aa=aa[,-1] #aggregate() adds this index useless column vector as first column, so remove it
    tmax=aggregate(a, by=list(days), max, na.rm=T) #find max instead of the mean to get daily tmax
    tmax=tmax[,4]
    aa=cbind(aa,tmax)
    aa[aa==-Inf]=NA #the tmax=aggregate() outputs some -Inf values here and there, so set them to NA
    #rm(a); rm(tmax)
    station.data=abind(station.data,aa)
  }  
  station.data=station.data[,,-1] #remove the NA.matrix used to initialize
  station.data=station.data[,-2,] #remove the 'hour' variable that we no longer need
  #dimensions of station.data are 23012 days x 4 variables (date, T, shum, tmax) x 96 stations
  dates=station.data[,1,]
  station.data=station.data[,-1,] #remove the 'date' variable b/c it needs to stay a string
  station.data=array(as.numeric(station.data),dim=dim(station.data)) #convert values from strings to numeric
  #saveRDS(station.data,'/scratch/users/tballard/shum/station.data/station.data.USOnly.rds')
  #saveRDS(dates, '/scratch/users/tballard/shum/station.data/station.data.USOnly.dates.rds')
  saveRDS(station.data,'/scratch/users/tballard/shum/station.data/station.data.SE.rds')
  saveRDS(dates, '/scratch/users/tballard/shum/station.data/station.data.SE.dates.rds')

  
  
  
  
  
  
  
  
  