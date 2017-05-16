suppressMessages(library(fields))
suppressMessages(library(ncdf4))
suppressMessages(library(RColorBrewer))
suppressMessages(library(abind))
suppressMessages(library(base))
suppressMessages(library(Rfit))
suppressMessages(library(MASS))
suppressMessages(library(gplots))

#62GB, 1hr
### This computes the 2d table of frequency for tmax and shum for reanalysis 1979-2014. 
dir.out="/scratch/users/tballard/shum/giorgi.regions/output/"
dir.plot="/scratch/users/tballard/shum/giorgi.regions/reanalysis/plots/"

#-------------------------------------------------------------------------------------
regions=readRDS("/scratch/users/tballard/shum/giorgi.regions/giorgi.regions.rds") #25 x 4
region.names=rownames(regions)

n.bins=75
xlims=c(-25,25)
ylims=c(-15,15)

yr.start=1
yr.end=36

n.lon=192
n.lat=94

tmax.hw=readRDS("/scratch/users/tballard/shum/percentile.threshold/daily.data/tmax.hw")
shum.hw=readRDS("/scratch/users/tballard/shum/percentile.threshold/daily.data/shum.hw")
lat=readRDS("/scratch/users/tballard/plots/global.attributes/shum/lat") 
lon=readRDS("/scratch/users/tballard/plots/global.attributes/shum/lon")
## Define month indices ##
jan=1:31; feb=32:59; mar=60:90; apr=91:120;may=121:151; jun=152:181;
jul=182:212; aug=213:243; sep=244:273; oct=274:304; nov=305:334; dec=335:365;


for (region in 1:dim(regions)[1]){
  print(region)
  for(i in 1:4){ 
    month.table=list(c(jan,feb,dec),c(mar,apr,may),c(jun,jul,aug),c(sep,oct,nov))
    season=month.table[[i]]
    season.names=c("DJF","MAM","JJA","SON")
    season.name=season.names[i]
    region.name=region.names[region]
    
    ## Extract the shum/tmax values for the month of interest.
    daily.shum=shum.hw[,,season,yr.start:yr.end] 
    dim(daily.shum) = c(dim(daily.shum)[1], dim(daily.shum)[2], dim(daily.shum)[3]*dim(daily.shum)[4]) #collapses from a 4dim array to 3dim array
    daily.tmax=tmax.hw[,,season,yr.start:yr.end] 
    dim(daily.tmax) = c(dim(daily.tmax)[1], dim(daily.tmax)[2], dim(daily.tmax)[3]*dim(daily.tmax)[4]) #collapses from a 4dim array to 3dim array
    
    
    ##### Read in ocean/land mask file 1's are land #####
    fileName="/scratch/PI/omramom/reanalysis/ncep-doe-r2/4xdaily/ocean.mask/land.sfc.gauss.nc"
    land = ncvar_get(nc_open(fileName), "land") #192 x 94
    land[land==0]=NA #set ocean values to NA instead of 0
    mask=function(data,land){
      new.data=land*data
      return(new.data)
    }
    tmax.mask=apply(daily.tmax, c(3), mask, land=land)
    tmax.mask=array(tmax.mask, dim=dim(daily.tmax))
    shum.mask=apply(daily.shum, c(3), mask, land=land)
    shum.mask=array(shum.mask, dim=dim(daily.shum))
    
    
    ##### Mask out everything except the region #####
    zone=regions[region,]
    ##### Read in lat/lon values you'll use to make the mask #####
    fileName="/scratch/PI/omramom/reanalysis/ncep-doe-r2/daily/shum/shum.2m.gauss.1979.nc"
    lon. = ncvar_get(nc_open(fileName), "lon") #192 values (0 to 358.12)
    lat. = ncvar_get(nc_open(fileName), "lat") #94 values (88.542 to -88.542)
    
    region.lat=lat.>zone[1] & lat.<zone[2] #True if lat is w/in bounds
    ##
    if (zone[3]<zone[4]){ #need to do a special for loop for longitude
      region.lon=lon.>zone[3] & lon.<zone[4] #True if lon is w/in bounds
    } else { #'else' part only happens if lon bounds include the 360:0 transition. 
      region.lon=lon.>zone[3] | lon.<zone[4]  
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
    mask=function(data, region.mask){
      new.data=region.mask*data
      return(new.data)
    }
    tmax.mask.p1=apply(tmax.mask, c(3), mask, region.mask=region.mask)
    tmax.mask.p1=array(tmax.mask.p1, dim=dim(daily.tmax))
    shum.mask.p1=apply(shum.mask, c(3), mask, region.mask=region.mask)
    shum.mask.p1=array(shum.mask.p1, dim=dim(daily.shum))
    
    
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
    
    
    matrix.fraction=function(tmax, shum, n.bins, x.grid, y.grid){
          counts=matrix(rep(NA, (n.bins-1)^2), nrow=n.bins-1)
          dims=n.bins-1
          for (i in 1:dims){
            for (j in 1:dims){
            counts[j,i]=sum(tmax>x.grid[i] & tmax<x.grid[i+1] & shum<y.grid[n.bins-j+1] & shum>y.grid[n.bins-j])
            }
          }
          matrix.fractions=counts/length(tmax) #make them into fractions instead of counts
          return(matrix.fractions)
        }
        
        x.grid=seq(xlims[1], xlims[2], length=n.bins)
        y.grid=seq(ylims[1], ylims[2], length=n.bins)
        
        fractions.1=matrix.fraction(tmax.val.1, shum.val.1, n.bins=n.bins, x.grid=x.grid, y.grid=y.grid) #table of fractions (counts / total number) within grid
        fractions.1[fractions.1==0]=NA
        saveRDS(fractions.1, paste(dir.out,'fractions.1979.2014.',region.names[region],'.',season.name,'.rds',sep=''))

        # q1=sum(tmax.val.1>0 & shum.val.1>0)/length(tmax.val.1) #fraction of values where both tmax and shum are positive (quadrant 1)
        # q2=sum(tmax.val.1<0 & shum.val.1>0)/length(tmax.val.1) #b/c I dont have an >= or <=, these 4 exclude instances where tmax or shum are 0, so they wont sum to 1
        # q3=sum(tmax.val.1<0 & shum.val.1<0)/length(tmax.val.1)
        # q4=sum(tmax.val.1>0 & shum.val.1<0)/length(tmax.val.1)
        # quadrant.fractions=c(q1,q2,q3,q4)
        # saveRDS(quadrant.fractions, paste(dir.out,'quadrant.fractions.1979.2014.',region.names[region],'.',season.name,'.rds',sep=''))

  }    
}
