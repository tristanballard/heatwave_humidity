suppressMessages(library(fields))
suppressMessages(library(ncdf4))
suppressMessages(library(RColorBrewer))
suppressMessages(library(abind))
suppressMessages(library(base))

### So here I'm importing the threshold tmax values calculated earlier and then
### creating a array of 1's and 0's for if the observed tmax exceeded those
### thresholds. I then multiply that "hw.index" matrix by the observed tmax
### and shum matrices (leap days removed) to get the values on those heat wave days
### with values on non-heat wave days set to NA "tmax.hw" and "shum.hw". 
yr.start=1979
yr.end=2010
n.years=yr.end-yr.start+1

#Read in threshold values
#threshold.final=readRDS("/scratch/users/tballard/shum/station.data/threshold.final")
#station.data=readRDS('/scratch/users/tballard/shum/station.data/station.data.USOnly.rds')
#dates=readRDS('/scratch/users/tballard/shum/station.data/station.data.USOnly.dates.rds')
threshold.final=readRDS("/scratch/users/tballard/shum/station.data/threshold.final.SE.rds")
station.data=readRDS('/scratch/users/tballard/shum/station.data/station.data.SE.rds')
dates=readRDS('/scratch/users/tballard/shum/station.data/station.data.SE.dates.rds')


#### Subset for the years of interest ####
t.start=paste(yr.start,"-01-01",sep=""); t.start=as.Date(t.start) #e.g. "1979-01-01"
t.end=paste(yr.end, "-12-30", sep=""); t.end=as.Date(t.end) #"2010-12-30
yr.index=dates[,1]>=t.start & dates[,1]<=t.end
station.data=station.data[yr.index,,]  #note there is some odd bug where the date for "1967-04-08" is NA 
station.data=station.data[-1,,] #for some reason, also an odd bug, a NA is at the start of the subsetted data, so remove it
dates=dates[yr.index,]
dates=dates[-1,]
## the bug is that the # of TRUE's in yr.index is 11324, but station.data2 has dim=11325. Very odd.


### Edit original station data so as to not include Leap Days #####
### The threshold calculation intentionally skipped leap days
  year=c(yr.start:yr.end)
  leap=paste(year,"-02-29",sep="")
  leap.position=match(leap,dates[,1]) #row value of the leap days, vector length=32
  leap.position=leap.position[!is.na(leap.position)] #NA's if there wasn't a leap year that year, so remove the NAs from this vector
  station.data=station.data[-leap.position,,]
  
  
###Here's another issue: station.data is now dim=11679x4x96; but 32yr * 365 = 11680 bc removed 2010 dec 31
###So add in an NA value for Dec. 31 2010 to make computations easier later. Won't matter since heat waves don't happen in Dec. anyways
dims=dim(station.data)
NA.matrix=matrix(rep(NA,dims[2]*dims[3]),c(dims[2],dims[3])) #3x89
station.data=abind(station.data,NA.matrix, along=1)
dims=dim(station.data)
  
###### if tmax exceeds threshold, set = to TRUE #####
station.data2=array(station.data,c(365,dims[1]/365,dims[2],dims[3])) #365x32x3x89
tmax=station.data2[,,3,] #365x32x89
tmax=aperm(tmax,c(3,1,2))
dim.tmax=dim(tmax) #89x365x32
#tmax=array(as.numeric(tmax),dim.tmax)
hw.index=array(rep(NA,dims[1]*dims[3]),dim.tmax) #initialize; 89x365x32
for (i in 1:dim.tmax[3]){ 
  hw.index[,,i]=tmax[,,i]>=threshold.final 
}  
#Now convert array of T/F to 1/0's
hw.index=hw.index+0
sum(hw.index,na.rm=T)/length(hw.index)
##### Now satisfy the constraint that you need 3+ days of 1's in a row. (12 values for 6hr data) #####
#the if statement is saying if it currently exceeds the threshold and the next 2 days do as well,
#keep it as a 1, else turn it to 0. But, if you do that the 2nd day in a 3 day event will fail and 
#be turned to a 0. So you have to include an OR clause that if it's a 1 and the preceding and following
#days are 1's, then it's a heatwave, and likewise if it's a 1 and the previous 2 days are 1's then
#it's also a heat wave. The simple example below shows it works.

#   a=c(1,1,1,0,1,1,0,0,0,1,1,1,1,1,0,1)
#   for (i in 3:(length(a)-2)){
#     if (a[i]==1 & a[i+1]==1 & a[i+2]==1
#         | a[i]==1 & a[i-1]==1 & a[i+1]==1
#         | a[i]==1 & a[i-1]==1 & a[i-2]==1){
#       a[i]=1 #keep it as a heat wave
#     }
#     else {
#       a[i]=0 #set as a non-event
#     }
#   }
#   a

#But now we need to do this for loop for all the data, yikes.
#For now ignore heat waves that span December 27-Jan 2nd
  is.hw=function(a){
    for (i in 3:(length(a)-11)){
      if (a[i]==1 & a[i+1]==1 & a[i+2]==1 
          | a[i]==1 & a[i-1]==1 & a[i+1]==1 
          | a[i]==1 & a[i-1]==1 & a[i-2]==1 
          ){
        a[i]=1 #keep it as a heat wave
      }
      else {
        a[i]=0 #set as a non-event
      }
    }
    return(a)
  }
  #So that function above cracks if there's NA's, so convert the NA's in hw.index (~1% of the data) to 0's 
  #which is fine b/c if the day is an NA then we don't have it as a heat event anyways
  hw.index[is.na(hw.index)]=0
  hw.index.2.0=apply(hw.index,c(1,3),is.hw)
  #saveRDS(hw.index.2.0,"/scratch/users/tballard/shum/station.data/hw.index")
  saveRDS(hw.index.2.0,"/scratch/users/tballard/shum/station.data/hw.index.SE.rds")

#hw.index=readRDS("/scratch/users/tballard/shum/station.data/hw.index") #365x96x32
hw.index=hw.index.2.0
hw.index=aperm(hw.index.2.0,c(1,3,2)) #change dimensions to 365x32x89

### Now apply the hw.index array of 1's and 0's to the observed tmax and shum   
tmax=station.data2[,,3,]; 
dim.tmax=dim(tmax) #365x32x89

shum=station.data2[,,2,] #365x32x89

tmax.hw=tmax*hw.index
tmax.hw[tmax.hw==0]=NA
#saveRDS(tmax.hw,"/scratch/users/tballard/shum/station.data/tmax.hw")
saveRDS(tmax.hw,"/scratch/users/tballard/shum/station.data/tmax.hw.SE.rds")

shum.hw=shum*hw.index
shum.hw[shum.hw==0]=NA #set all days that arent heat waves to NA instead of 0
#saveRDS(shum.hw,"/scratch/users/tballard/shum/station.data/shum.hw")
saveRDS(shum.hw,"/scratch/users/tballard/shum/station.data/shum.hw.SE.rds")

### Now do the same thing for shum, but create a matrix where it's NA's during
### the heat waves and observed values for non-heat wave days
hw.index.not=hw.index
hw.index.not[hw.index.not==0]=2 #there is probably a more clever way to swap 0's and 1's but w/e
hw.index.not[hw.index.not==1]=0
hw.index.not[hw.index.not==2]=1

shum.hw.not=shum*hw.index.not
saveRDS(shum.hw.not,"/scratch/users/tballard/shum/station.data/shum.hw.not")


