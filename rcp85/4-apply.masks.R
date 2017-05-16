suppressMessages(library(fields))
suppressMessages(library(ncdf4))
suppressMessages(library(RColorBrewer))
suppressMessages(library(abind))
suppressMessages(library(base))
suppressMessages(library(Rfit))
suppressMessages(library(MASS))
suppressMessages(library(gplots))

#request 140GB, takes hr
#check on ACCESS1-3 b/c it's reporting 0 hw events for rcp85 and way too many for historical
#run this on the args=historical args, not the rcp85 args
##model.name.rcp85=NA if it was a run/model that existed for historical but not rcp85. If there;s model/runs available for rcp85 but not historical, this code won't run that
## For non-plotting, I run the historical and rcp85 sections separately (though in same script) b/c there are unequal numbers of runs and it's the only way to get those extra runs historical has
## For plotting you are comparing model/run to model/run so it's fine. Run it on the rcp85 args though.

args=(commandArgs(TRUE)) #read in command line prompt index (args will equal, say '6'), indicating which model to run this on
args=as.numeric(args)
print(args)

dir='/scratch/users/tballard/shum/percentile.threshold/daily.data/cmip5/historical/'
dir.rcp85='/scratch/users/tballard/shum/percentile.threshold/daily.data/cmip5/rcp85/'
dir.out='/scratch/users/tballard/shum/giorgi.regions/output/'
threshold.filenames=read.csv(paste(dir,'threshold.filenames.csv',sep=""), head=T)

model=as.character(threshold.filenames[args[2],1])
model.name=unlist(strsplit(model, split='_', fixed=TRUE))[3]
run=unlist(strsplit(model, split='_', fixed=TRUE))[5]
print(paste(model.name,'.',run,sep=''))

model.rcp85=as.character(threshold.filenames[args[2],7])
model.name.rcp85=unlist(strsplit(model.rcp85, split='_', fixed=TRUE))[3]
run.rcp85=unlist(strsplit(model.rcp85, split='_', fixed=TRUE))[5]
print(paste(model.name.rcp85,'.',run.rcp85,sep=''))

yr.start=1 #1=1970
yr.end=36 #36=2005
n.years=yr.end-yr.start+1
yr.start2=76 #36=2041; 56=2061; 76=2081; 60=2065
yr.end2=95 #55=2060; 75=2080; 95=2100
n.years2=yr.end2-yr.start2+1

file.start=1970 #used in file naming of output
file.start2=2006
yr.start.real=file.start+yr.start-1
yr.end.real=file.start+yr.end-1
yr.start.real2=file.start2+yr.start2-1
yr.end.real2=file.start2+yr.end2-1

#-------------------------------------------------------------------------------------
regions=readRDS("/scratch/users/tballard/shum/giorgi.regions/giorgi.regions.rds") #25 x 4
region.names=rownames(regions)

#-------------------------------------------------------------------------------------

tmax.hw=readRDS(paste(dir,model,".tmax.hw.1970.2005",sep="")) #360*180*365*36, celsius
shum.hw=readRDS(paste(dir,model,".shum.hw.1970.2005",sep=""))*1000

if (!is.na(model.name.rcp85)){
tmax.hw.rcp85=readRDS(paste(dir.rcp85,model.rcp85,".tmax.hw.2006.2100",sep="")) #360*180*365*95, celsius
shum.hw.rcp85=readRDS(paste(dir.rcp85,model.rcp85,".shum.hw.2006.2100",sep=""))*1000
}

lon = readRDS((paste(dir,'lon.cmip5.rds',sep=""))) #n.lon values (0 to 359)
lat = readRDS((paste(dir,'lat.cmip5.rds',sep=""))) #n.lat values (-89.5 to 89.5)
n.lon=length(lon); n.lat=length(lat)

## Define month indices ##
jan=1:31; feb=32:59; mar=60:90; apr=91:120;may=121:151; jun=152:181;
jul=182:212; aug=213:243; sep=244:273; oct=274:304; nov=305:334; dec=335:365;
month.index=list(jan,feb,mar,apr,may,jun,jul,aug,sep,oct,nov,dec)
month.full=c("January","February","March","April","May","June","July","August","September","October","November","December")

#-------------------------------------------------------------------------------------


for(i in args[1]:args[1]){ 
  print(i)
  month.table=list(c(jan,feb,dec),c(mar,apr,may),c(jun,jul,aug),c(sep,oct,nov))
  season=month.table[[i]]
  season.names=c("DJF","MAM","JJA","SON")
  season.name=season.names[i]
  
  #### First extract data and apply masks for the first period #####
  ## Extract the shum/tmax values for the month of interest. 
  daily.tmax=array(tmax.hw[,,season,yr.start:yr.end],dim=c(n.lon,n.lat,n.years*length(season)))
  daily.shum=array(shum.hw[,,season,yr.start:yr.end],dim=c(n.lon,n.lat,n.years*length(season)))
  #dimensions = n.lon x n.lat x time (time=n.years (18) * n.days in month)
  #Now extract the same but for the 2nd time period; a bit more memory intensive this way
  if (!is.na(model.name.rcp85)){
  daily.tmax.2=array(tmax.hw.rcp85[,,season,yr.start2:yr.end2],dim=c(n.lon,n.lat,n.years2*length(season)))
  daily.shum.2=array(shum.hw.rcp85[,,season,yr.start2:yr.end2],dim=c(n.lon,n.lat,n.years2*length(season)))
  }
  
  for (region in 1:dim(regions)[1]){
    region.name=region.names[region]
    
    ##### Read in ocean/land mask file #####
    ##### Now mask out the ocean #####  
    fileName="/scratch/users/tballard/shum/percentile.threshold/daily.data/cmip5/lsmask.nc"
    land = ncvar_get(nc_open(fileName), "mask") #n.lon x n.lat
    land=land[,rev(seq_len(ncol(land)))] #note that this mask has the lat values upside down compared w/ our files, so flip it bck
    land=land+1 #swap 1's for water NA for land (original) with the other way round
    land[land==2]=NA
    
    tmax.mask = sweep(daily.tmax, 1, land, "*", check.margin=F)
    shum.mask = sweep(daily.shum, 1, land, "*", check.margin=F)
    
    ##### Mask out as well all the globe not in the US and southern Canada #####
    zone=regions[region,]
    ##### Define the box of interest for the region #####
    region.lat=lat>zone[1] & lat<zone[2] #True if lat is w/in bounds
    ##
    if (zone[3]<zone[4]){ #need to do a special for loop for longitude
      region.lon=lon>zone[3] & lon<zone[4] #True if lon is w/in bounds
    } else { #'else' part only happens if lon bounds include the 360:0 transition. 
      region.lon=lon>zone[3] | lon<zone[4]  
    }
    
    region.lat[region.lat==FALSE]=NA #Set False values to NA
    region.lon[region.lon==FALSE]=NA
    region.lat=region.lat+0; region.lon=region.lon+0 #Set TRUE values to 1
    
    #Now that you have a vector for lat and vector for lon of NA or 1's,
    #Combine that into a matrix format that you can multiply things by
    region.lat.mat=matrix(rep(region.lat,n.lon), nrow=n.lon, byrow=T)
    region.lon.mat=matrix(rep(region.lon,n.lat), nrow=n.lon)
    region.mask=region.lat.mat*region.lon.mat #n.lonxn.lat
    
    ####Below is applying the mask to your entire array of values
    tmax.mask.p1 = sweep(tmax.mask, 1, region.mask, "*", check.margin=F)
    shum.mask.p1 = sweep(shum.mask, 1, region.mask, "*", check.margin=F)
    
    
    ################################################
    #######  Repeat for the 2nd Time period  #######
    ################################################
    if (!is.na(model.name.rcp85)){
    rm(tmax.mask); rm(shum.mask)
    tmax.mask = sweep(daily.tmax.2, 1, land, "*", check.margin=F)
    shum.mask = sweep(daily.shum.2, 1, land, "*", check.margin=F)
    
    tmax.mask.p2 = sweep(tmax.mask, 1, region.mask, "*", check.margin=F)
    shum.mask.p2 = sweep(shum.mask, 1, region.mask, "*", check.margin=F)
    }
    ##########################################################  
    ##### Now that we have the grid of masked values,      ###  
    ##### extract only the non-NA's and non-0's (0s happen ###
    ##### over grid cells if there wasn't a heat wave)     ###
    ##########################################################
    tmax.med=apply(tmax.mask.p1, c(1,2), median, na.rm=T) #360 x 180
    tmax.med[is.na(tmax.med)]=0
    tmax.val.1=sweep(tmax.mask.p1, c(1,2), tmax.med, FUN="-")
    
    tmax.val.1=as.vector(tmax.val.1); 
    tmax.val.1=tmax.val.1[!is.na(tmax.val.1)]
    #tmax.val.1=tmax.val.1[tmax.val.1!=0]
    
    shum.med=apply(shum.mask.p1, c(1,2), median, na.rm=T) #360 x 180
    shum.med[is.na(shum.med)]=0
    shum.val.1=sweep(shum.mask.p1, c(1,2), shum.med, FUN="-")

    shum.val.1=as.vector(shum.val.1); 
    shum.val.1=shum.val.1[!is.na(shum.val.1)]
    #shum.val.1=shum.val.1[shum.val.1!=0]
    
    if (!is.na(model.name.rcp85)){
   
    tmax.val.2=sweep(tmax.mask.p2, c(1,2), tmax.med, FUN="-") #subtract the pixel by pixel medians from the first time period
    tmax.val.2=as.vector(tmax.val.2);
    tmax.val.2=tmax.val.2[!is.na(tmax.val.2)]
    #tmax.val.2=tmax.val.2[tmax.val.2!=0]
    
    shum.val.2=sweep(shum.mask.p2, c(1,2), shum.med, FUN="-") #subtract the pixel by pixel medians from the first time period
    shum.val.2=as.vector(shum.val.2); 
    shum.val.2=shum.val.2[!is.na(shum.val.2)]
    #shum.val.2=shum.val.2[shum.val.2!=0]
    }
    
    if (!is.na(model.name.rcp85)){
    saveRDS(tmax.val.2, paste(dir.out, model.name,'.',run,'.tmax.val.2.',region.names[region],'.',season.name,'.',yr.start.real2,'.',yr.end.real2,'.rds',sep=''))
    saveRDS(shum.val.2, paste(dir.out, model.name,'.',run,'.shum.val.2.',region.names[region],'.',season.name,'.',yr.start.real2,'.',yr.end.real2,'.rds',sep=''))
    }
    saveRDS(tmax.val.1, paste(dir.out, model.name,'.',run,'.tmax.val.1.',region.names[region],'.',season.name,'.',yr.start.real,'.',yr.end.real,'.rds',sep=''))
    saveRDS(shum.val.1, paste(dir.out, model.name,'.',run,'.shum.val.1.',region.names[region],'.',season.name,'.',yr.start.real,'.',yr.end.real,'.rds',sep=''))

    
    
    ## Now compute the # of hw events in each period
    n.events.1=length(tmax.val.1) #length is the same for shum.val.1
    saveRDS(n.events.1, paste(dir.out, model.name,'.',run,'.n.events.1.',region.names[region],'.',season.name,'.',yr.start.real,'.',yr.end.real,'.rds',sep=''))
    
    if (!is.na(model.name.rcp85)){
    n.events.2=length(tmax.val.2)
    saveRDS(n.events.2, paste(dir.out, model.name.rcp85,'.',run.rcp85,'.n.events.2.',region.names[region],'.',season.name,'.',yr.start.real2,'.',yr.end.real2,'.rds',sep=''))
    }
  
  }
}    
