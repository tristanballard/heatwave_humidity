suppressMessages(library(fields))
suppressMessages(library(ncdf4))
suppressMessages(library(RColorBrewer))
suppressMessages(library(abind))
suppressMessages(library(base))
suppressMessages(library(Rfit))
suppressMessages(library(MASS))
suppressMessages(library(gplots))

#62GB, 18min
### This computes the 2d table of frequency for tmax and shum for reanalysis 1979-2014 over the 2 half periods within, as well as the difference.
dir.out="/scratch/users/tballard/shum/giorgi.regions/output/"
dir.plot="/scratch/users/tballard/shum/giorgi.regions/reanalysis/plots/"

#-------------------------------------------------------------------------------------
regions=readRDS("/scratch/users/tballard/shum/giorgi.regions/giorgi.regions.rds") #25 x 4
region.names=rownames(regions)

n.bins=75
xlims=c(-25,25)
ylims=c(-15,15)


yr.start1=1979
yr.end1=1996
yr.start2=1997
yr.end2=2014

n.lon=192
n.lat=94


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
    
    tmax.val.1=readRDS(paste(dir.out,'tmax.val.1.',region.names[region],'.',season.name,'.rds',sep=''))
    shum.val.1=readRDS(paste(dir.out,'shum.val.1.',region.names[region],'.',season.name,'.rds',sep=''))
    tmax.val.2=readRDS(paste(dir.out,'tmax.val.2.',region.names[region],'.',season.name,'.rds',sep=''))
    shum.val.2=readRDS(paste(dir.out,'shum.val.2.',region.names[region],'.',season.name,'.rds',sep=''))
     
    
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
        fractions.2=matrix.fraction(tmax.val.2, shum.val.2, n.bins=n.bins, x.grid=x.grid, y.grid=y.grid) #table of fractions (counts / total number) within grid
        fractions.diff=fractions.2-fractions.1
        fractions.diff[fractions.diff==0]=NA
        fractions.1[fractions.1==0]=NA
        fractions.2[fractions.2==0]=NA
   
        saveRDS(fractions.1, paste(dir.out,'fractions.',yr.start1,'.',yr.end1,'.',region.names[region],'.',season.name,'.rds',sep=''))
        saveRDS(fractions.2, paste(dir.out,'fractions.',yr.start2,'.',yr.end2,'.',region.names[region],'.',season.name,'.rds',sep=''))
        saveRDS(fractions.diff, paste(dir.out,'fractions.diff',yr.start1,'.',yr.end2,'.',region.names[region],'.',season.name,'.rds',sep=''))

        ### Get the fractions within each quadrant for the joint pdf plot
        # q1=sum(tmax.val.1>0 & shum.val.1>0)/length(tmax.val.1) #fraction of values where both tmax and shum are positive (quadrant 1)
        # q2=sum(tmax.val.1<0 & shum.val.1>0)/length(tmax.val.1) #b/c I dont have an >= or <=, these 4 exclude instances where tmax or shum are 0, so they wont sum to 1
        # q3=sum(tmax.val.1<0 & shum.val.1<0)/length(tmax.val.1)
        # q4=sum(tmax.val.1>0 & shum.val.1<0)/length(tmax.val.1)
        # quadrant.fractions=c(q1,q2,q3,q4)
        # saveRDS(quadrant.fractions.1, paste(dir.out,'quadrant.fractions.',yr.start1,'.',yr.end1,'.',region.names[region],'.',season.name,'.rds',sep='')); rm(quadrant.fractions)
        # 
        # q1=sum(tmax.val.2>0 & shum.val.2>0)/length(tmax.val.2) #fraction of values where both tmax and shum are positive (quadrant 1)
        # q2=sum(tmax.val.2<0 & shum.val.2>0)/length(tmax.val.2) #b/c I dont have an >= or <=, these 4 exclude instances where tmax or shum are 0, so they wont sum to 1
        # q3=sum(tmax.val.2<0 & shum.val.2<0)/length(tmax.val.2)
        # q4=sum(tmax.val.2>0 & shum.val.2<0)/length(tmax.val.2)
        # quadrant.fractions=c(q1,q2,q3,q4)
        # saveRDS(quadrant.fractions, paste(dir.out,'quadrant.fractions.',yr.start2,'.',yr.end2,'.',region.names[region],'.',season.name,'.rds',sep=''))
        # 

  }    
}
