start.time=Sys.time()
suppressMessages(library(fields))
suppressMessages(library(ncdf4))
suppressMessages(library(RColorBrewer))
suppressMessages(library(abind))
suppressMessages(library(base))
suppressMessages(library(Rfit))
suppressMessages(library(MASS))
suppressMessages(library(gplots))
suppressMessages(library(data.table))
print("Let's Begin")

### This computes the 2d table of frequency for tmax and shum for cmip historical 1970-87, 1988-2005, and the difference. The tmax and shum
### for each region/season combo have been created in 4-apply.masks.R
### Saves the output for every model/run combo in one array. This script is run on every season/region combo.

#request 62GB takes 4.5hr
args=(commandArgs(TRUE)) #read in command line prompt index 
args=as.numeric(args)
print(args) #arg = 1, 2, 3, or 4 for which season to run this on, and then either 1, 6 for first 6 regions or 7, 12 for 2nd set ...

dir='/scratch/users/tballard/shum/percentile.threshold/daily.data/cmip5/historical/'
dir.out='/scratch/users/tballard/shum/giorgi.regions/output/'
threshold.filenames=read.csv(paste(dir,'threshold.filenames.csv',sep=""), head=T)
#threshold.filenames[3,1]="NA" #Setting ACCESS1-3 to NA so that it won't read in its data
threshold.filenames[c(3,4,5,6,42:57,64:69,73:76,82:85,89),1]="NA" #Setting ACCESS1-3, and a few others to NA so that it won't read in its data; this throws a warning but it's OK

n.bins=75
xlims=c(-25,25)
ylims=c(-15,15)

yr.start.1=1970
yr.end.1=1987
yr.start.2=1988
yr.end.2=2005
#-------------------------------------------------------------------------------------
regions=readRDS("/scratch/users/tballard/shum/giorgi.regions/giorgi.regions.rds") #25 x 4
region.names=rownames(regions)
#-------------------------------------------------------------------------------------

lon = readRDS((paste(dir,'lon.cmip5.rds',sep=""))) #n.lon values (0 to 359)
lat = readRDS((paste(dir,'lat.cmip5.rds',sep=""))) #n.lat values (-89.5 to 89.5)
n.lon=length(lon); n.lat=length(lat)

## Define month indices ##
jan=1:31; feb=32:59; mar=60:90; apr=91:120;may=121:151; jun=152:181;
jul=182:212; aug=213:243; sep=244:273; oct=274:304; nov=305:334; dec=335:365;
month.index=list(jan,feb,mar,apr,may,jun,jul,aug,sep,oct,nov,dec)
month.full=c("January","February","March","April","May","June","July","August","September","October","November","December")

#-------------------------------------------------------------------------------------
##Get model names and run, e.g. 'ACCESS1-0' and run='r2i1p1' that'll correspond to data you're reading in
model=as.character(threshold.filenames[,1])
parse.name=function(a,idx){
  aa=unlist(strsplit(a, split='_', fixed=TRUE))[idx]
  if(is.na(aa)){
    aa=c()
  }
  return(aa)
}
model.names=unlist(lapply(model,parse.name,idx=3))
run=unlist(lapply(model,parse.name,idx=5))

#-------------------------------------------------------------------------------------
#### Run the actual meat of the code now:
for(i in args[1]:args[1]){
  print(paste("Season:",i))
  month.table=list(c(jan,feb,dec),c(mar,apr,may),c(jun,jul,aug),c(sep,oct,nov))
  season=month.table[[i]] 
  season.names=c("DJF","MAM","JJA","SON") 
  season.name=season.names[i] 
  
  for (region in args[2]:args[3]){
    print(paste("Region:",region))
    region.name=region.names[region]
    
    NA.matrix=matrix(rep(NA, (n.bins-1)^2), nrow=n.bins-1) #used to initialize
    fractions.diff.all=NA.matrix #initialize
    fractions.1.all=NA.matrix #initialize
    fractions.2.all=NA.matrix #initialize
    quadrant.fractions.1.all=c()
    quadrant.fractions.2.all=c()

    for (i in 1:length(model.names)){
      #print(paste(model.names[i], '.', run[i]))
      filename.1=paste(dir.out, model.names[i],'.',run[i],'.tmax.val.1.',region.names[region],'.',season.name,'.',yr.start.1,'.',yr.end.1,'.rds',sep='')       
      filename.2=paste(dir.out, model.names[i],'.',run[i],'.shum.val.1.',region.names[region],'.',season.name,'.',yr.start.1,'.',yr.end.1,'.rds',sep='')
      filename.3=paste(dir.out, model.names[i],'.',run[i],'.tmax.val.2.',region.names[region],'.',season.name,'.',yr.start.2,'.',yr.end.2,'.rds',sep='')        
      filename.4=paste(dir.out, model.names[i],'.',run[i],'.shum.val.2.',region.names[region],'.',season.name,'.',yr.start.2,'.',yr.end.2,'.rds',sep='')

      if(file.exists(filename.1) & file.exists(filename.2) & file.exists(filename.3) & file.exists(filename.4)){
        tmax.val.1=readRDS(filename.1)
        shum.val.1=readRDS(filename.2)
        tmax.val.2=readRDS(filename.3)
        shum.val.2=readRDS(filename.4)
        
        matrix.fraction=function(tmax, shum, n.bins, x.grid, y.grid){
          counts=matrix(rep(NA, (n.bins-1)^2), nrow=n.bins-1)
          dims=n.bins-1
          for (i in 1:dims){
            for (j in 1:dims){
            counts[j,i]=sum(tmax>x.grid[i] & tmax<x.grid[i+1] & shum<y.grid[n.bins-j+1] & shum>y.grid[n.bins-j])
            }
          }
          matrix.fractions=counts/length(tmax)
          return(matrix.fractions)
        }
        
        x.grid=seq(xlims[1], xlims[2], length=n.bins)
        y.grid=seq(ylims[1], ylims[2], length=n.bins)
        
        fractions.1=matrix.fraction(tmax.val.1, shum.val.1, n.bins=n.bins, x.grid=x.grid, y.grid=y.grid) #Period 1 table of fractions (counts / total number) within grid
        fractions.2=matrix.fraction(tmax.val.2, shum.val.2, n.bins=n.bins, x.grid=x.grid, y.grid=y.grid) #Period 2 table of fractions (counts / total number) within grid
        
        fractions.diff = fractions.2-fractions.1
        fractions.diff[fractions.diff==0]=NA
        fractions.1[fractions.1==0]=NA
        fractions.2[fractions.2==0]=NA
        
        fractions.1.all = abind(fractions.1.all, fractions.1)
        fractions.2.all = abind(fractions.2.all, fractions.2)
        fractions.diff.all = abind(fractions.diff.all, fractions.diff)      
      
        ### Get the fractions within each quadrant for the joint pdf plot
        # q1=sum(tmax.val.1>0 & shum.val.1>0)/length(tmax.val.1) #fraction of values where both tmax and shum are positive (quadrant 1)
        # q2=sum(tmax.val.1<0 & shum.val.1>0)/length(tmax.val.1) #b/c I dont have an >= or <=, these 4 exclude instances where tmax or shum are 0, so they wont sum to 1
        # q3=sum(tmax.val.1<0 & shum.val.1<0)/length(tmax.val.1)
        # q4=sum(tmax.val.1>0 & shum.val.1<0)/length(tmax.val.1)
        # quadrant.fractions=c(q1,q2,q3,q4)
        # quadrant.fractions.1.all=rbind(quadrant.fractions.1.all, quadrant.fractions)
        # 
        # q1=sum(tmax.val.2>0 & shum.val.2>0)/length(tmax.val.2) #fraction of values where both tmax and shum are positive (quadrant 1)
        # q2=sum(tmax.val.2<0 & shum.val.2>0)/length(tmax.val.2) #b/c I dont have an >= or <=, these 4 exclude instances where tmax or shum are 0, so they wont sum to 1
        # q3=sum(tmax.val.2<0 & shum.val.2<0)/length(tmax.val.2)
        # q4=sum(tmax.val.2>0 & shum.val.2<0)/length(tmax.val.2)
        # quadrant.fractions=c(q1,q2,q3,q4)
        # quadrant.fractions.2.all=rbind(quadrant.fractions.2.all, quadrant.fractions)

      }
      
      else {
        fractions.diff.all = abind(fractions.diff.all, NA.matrix)
        fractions.1.all = abind(fractions.1.all, NA.matrix)
        fractions.2.all = abind(fractions.2.all, NA.matrix)
        # quadrant.fractions.1.all=rbind(quadrant.fractions.1.all, c(NA,NA,NA,NA))
        # quadrant.fractions.2.all=rbind(quadrant.fractions.2.all, c(NA,NA,NA,NA))
      }
      
    }
      fractions.diff.all=array(fractions.diff.all, dim=c(n.bins-1, n.bins-1, length(model.names)+1))
      fractions.diff.all=fractions.diff.all[,,-1] #remove initializing NA matrix
      fractions.1.all=array(fractions.1.all, dim=c(n.bins-1, n.bins-1, length(model.names)+1))
      fractions.1.all=fractions.1.all[,,-1] #remove initializing NA matrix
      fractions.2.all=array(fractions.2.all, dim=c(n.bins-1, n.bins-1, length(model.names)+1))
      fractions.2.all=fractions.2.all[,,-1] #remove initializing NA matrix
      
      saveRDS(fractions.diff.all, paste(dir.out, 'fractions.diff.all.',region.names[region],'.',season.name,'.',yr.start.1,'.',yr.end.1,'.',yr.start.2,'.',yr.end.2,'.rds',sep=''))
      saveRDS(fractions.1.all, paste(dir.out, 'fractions.1.all.',region.names[region],'.',season.name,'.',yr.start.1,'.',yr.end.1,'.rds',sep=''))
      saveRDS(fractions.2.all, paste(dir.out, 'fractions.2.all.',region.names[region],'.',season.name,'.',yr.start.2,'.',yr.end.2,'.rds',sep=''))
      # saveRDS(quadrant.fractions.1.all, paste(dir.out, 'quadrant.fractions.1.all.',region.names[region],'.',season.name,'.',yr.start.1,'.',yr.end.1,'.rds',sep=''))
      # saveRDS(quadrant.fractions.2.all, paste(dir.out, 'quadrant.fractions.2.all.',region.names[region],'.',season.name,'.',yr.start.2,'.',yr.end.2,'.rds',sep=''))
      # 
  
  print(paste("Time elapsed: ", Sys.time()-start.time))
  }
}