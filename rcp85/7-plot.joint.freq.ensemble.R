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
#request 62, goes very fast
# args=(commandArgs(TRUE)) #read in command line prompt index 
# args=as.numeric(args)
# print(args)

dir='/scratch/users/tballard/shum/percentile.threshold/daily.data/cmip5/historical/'
dir.rcp85='/scratch/users/tballard/shum/percentile.threshold/daily.data/cmip5/rcp85/'
dir.out='/scratch/users/tballard/shum/giorgi.regions/output/'
dir.plot='/scratch/users/tballard/shum/giorgi.regions/cmip5/rcp85/plots/'
threshold.filenames=read.csv(paste(dir,'threshold.filenames.csv',sep=""), head=T)
#threshold.filenames[3,1]="NA" #Setting ACCESS1-3 to NA so that it won't read in its data
threshold.filenames[c(3,4,5,6,42:57,64:69,73:76,82:85,89),1]="NA" #Setting ACCESS1-3 to NA so that it won't read in its data

n.bins=75 #the quadrant code will fail if this is an even number
xlims=c(-25,25)
ylims=c(-15,15)

yr.start=1 #1=1970
yr.end=36 #36=2005
n.years=yr.end-yr.start+1
#yr.start2=args[1] #36=2041; 56=2061; 76=2081; 60=2065
#yr.end2=args[2] #55=2060; 75=2080; 95=2100
#n.years2=yr.end2-yr.start2+1

yr.start.obs=1979
yr.end.obs=2014
file.start=1970 #used in file naming of output
file.start2=2006
yr.start.real=file.start+yr.start-1
yr.end.real=file.start+yr.end-1
#yr.start.real2=file.start2+yr.start2-1
#yr.end.real2=file.start2+yr.end2-1
#print(paste(yr.start.real2,"-",yr.end.real2))
yr.start.rcp=c(2041,2061,2081) #redundant I know
yr.end.rcp=c(2060,2080,2100)
yr.start.hist=c(1970)
yr.end.hist=c(2005)
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
    x=as.factor(model.names)
    levels(x)=1:length(levels(x))
    n.unique=length(levels(x)) #number of unique models, regardless of # of runs
    x=as.numeric(x) #instead of access 1-0, access1-0, access 1-3, it's c(1,1,2...)

#-------------------------------------------------------------------------------------
#### Run the actual meat of the code now:
for (region in 1:dim(regions)[1]){
    print(paste("Region:",region))
    region.name=region.names[region]
    
    ##So sloppy for loop things happening. Here, go through and take ensemble medians for the 3 time periods
      ##of each season, combine into one array. Then as i goes from 1 to 4, combine all those arrays to get a n.bins-1 x n.bins-1 x 12.
      NA.matrix=matrix(rep(NA,(n.bins-1)^2),nrow=n.bins-1) #used to initialize
      mm.ens.season=NA.matrix 
      q.ens.season=matrix(rep(NA,8), nrow=4)
    for(i in 1:4){
      print(paste("Season:",i))
      month.table=list(c(jan,feb,dec),c(mar,apr,may),c(jun,jul,aug),c(sep,oct,nov))
      season=month.table[[i]] 
      season.names=c("DJF","MAM","JJA","SON") 
      season.name=season.names[i] 
      
          fractions.ensemble=function(fractions.diff.all, n.unique, x, n.bins){
              NA.matrix=matrix(rep(NA,(n.bins-1)^2),nrow=n.bins-1) #used to initialize
              fractions.ens.med=NA.matrix
            ###Now ensemble the tables for models with multiple runs
            for (i in 1:n.unique){
              seq.start=head(which(x==i), n=1) #e.g. for x=c(1,1,2,2,2,3), and i=2, seq.start=3, seq.end=5
              seq.end=tail(which(x==i), n=1)
              if (sum(x==i)>1){ #if there's more than one run for each unique model
              a=apply(fractions.diff.all[,,seq.start:seq.end], c(1,2), median, na.rm=T) #take the median of them
              fractions.ens.med=abind(fractions.ens.med, a)
              } else { #Otherwise you dont need to take the median (only one grid), so just add it on
                fractions.ens.med = abind(fractions.ens.med, fractions.diff.all[,,seq.start:seq.end])
              }
            }
            fractions.ens.med=array(fractions.ens.med, dim=c(n.bins-1, n.bins-1, n.unique+1)) #fix dimensions
            fractions.ens.med=fractions.ens.med[,,-1] #dim = n.bins=1  x n.bins-1 x n.unique
            ##Now that there's 1 table per model, take the final ensemble median:
            fractions.mm.ens.med1=apply(fractions.ens.med, c(1,2), median, na.rm=T) #n.bins-1 x n.bins-1
            fractions.mm.ens.med1[fractions.mm.ens.med1==0]=NA #not needed for the rcp files (already set as NA), but needed for observations/cmip historical
            return(fractions.mm.ens.med1)
          }
         
          ### chance to debug: one of the 15 unique models gives a sum of exactly 0 for q1, q2, q3, and q4 
          quadrant.ensemble=function(fractions.diff.all, n.unique, x, n.bins){
              NA.matrix=matrix(rep(NA,(n.bins-1)^2),nrow=n.bins-1) #used to initialize
              fractions.ens.med=NA.matrix
            ###Now ensemble the tables for models with multiple runs
            for (i in 1:n.unique){
              seq.start=head(which(x==i), n=1) #e.g. for x=c(1,1,2,2,2,3), and i=2, seq.start=3, seq.end=5
              seq.end=tail(which(x==i), n=1)
              if (sum(x==i)>1){ #if there's more than one run for each unique model
              a=apply(fractions.diff.all[,,seq.start:seq.end], c(1,2), median, na.rm=T) #take the median of them
              fractions.ens.med=abind(fractions.ens.med, a)
              } else { #Otherwise you dont need to take the median (only one grid), so just add it on
                fractions.ens.med = abind(fractions.ens.med, fractions.diff.all[,,seq.start:seq.end])
              }
            }
            fractions.ens.med=array(fractions.ens.med, dim=c(n.bins-1, n.bins-1, n.unique+1)) #fix dimensions
            fractions.ens.med=fractions.ens.med[,,-1] #dim = n.bins=1  x n.bins-1 x n.unique
            
            ##Now take the sum in each quadrant
            fractions.ens.med[fractions.ens.med==0]=NA #repeated elsewhere but helpful just in case
            breaks=(n.bins-1)/2
            q1=c(); q2=c(); q3=c(); q4=c(); 
            for (i in 1:n.unique){ #get the sum for each quadrant for each model
              q1 = c(q1, sum(fractions.ens.med[1:breaks,(breaks+1):(n.bins-1),i], na.rm=T))
              q2 = c(q2, sum(fractions.ens.med[1:breaks,1:breaks,i], na.rm=T))
              q3 = c(q3, sum(fractions.ens.med[(breaks+1):(n.bins-1),1:breaks,i], na.rm=T))
              q4 = c(q4, sum(fractions.ens.med[(breaks+1):(n.bins-1),(breaks+1):(n.bins-1),i], na.rm=T))
            }
            q1.m=median(q1, na.rm=1); q2.m=median(q2, na.rm=1);  q3.m=median(q3, na.rm=1);  q4.m=median(q4, na.rm=1);  
            if(q1.m>0){
              q1.a=mean(q1>0, na.rm=T)*100 #percent that agree in sign 
            } else{
              q1.a=mean(q1<0, na.rm=T)*100 #percent that agree in sign 
            }
            
            if(q2.m>0){
              q2.a=mean(q2>0, na.rm=T)*100 #percent that agree in sign 
            } else{
              q2.a=mean(q2<0, na.rm=T)*100 #percent that agree in sign 
            }
            
            if(q3.m>0){
              q3.a=mean(q3>0, na.rm=T)*100 #percent that agree in sign 
            } else{
              q3.a=mean(q3<0, na.rm=T)*100 #percent that agree in sign 
            }
            
            if(q4.m>0){
              q4.a=mean(q4>0, na.rm=T)*100 #percent that agree in sign 
            }  else{
              q4.a=mean(q4<0, na.rm=T)*100 #percent that agree in sign 
            }
            
            q.results=cbind(c(q1.m, q2.m, q3.m, q4.m), c(q1.a, q2.a, q3.a, q4.a)) #first column is the median value in the quadrant, 2nd column is the %models agreeing in sign
            return(q.results)
          }
          
           #RCP time period 1:
            #Read in joint frequency table from part 6, an array of tables for all the model/runs available
            fractions.diff.all=readRDS(paste(dir.out, 'fractions.diff.all.',region.name,'.',season.name,'.',yr.start.real,'.',yr.end.real,'.',yr.start.rcp[1],'.',yr.end.rcp[1],'.rds',sep='')) #dim = n.bins-1 x n.bins-1 x length(x)
            fractions.mm.ens.med1=fractions.ensemble(fractions.diff.all, n.unique, x, n.bins)
            q.ens1=quadrant.ensemble(fractions.diff.all, n.unique, x, n.bins)
            
           #RCP time period 2:
            #Read in joint frequency table from part 6, an array of tables for all the model/runs available
            fractions.diff.all=readRDS(paste(dir.out, 'fractions.diff.all.',region.name,'.',season.name,'.',yr.start.real,'.',yr.end.real,'.',yr.start.rcp[2],'.',yr.end.rcp[2],'.rds',sep='')) #dim = n.bins-1 x n.bins-1 x length(x)
            fractions.mm.ens.med2=fractions.ensemble(fractions.diff.all, n.unique, x, n.bins)
            q.ens2=quadrant.ensemble(fractions.diff.all, n.unique, x, n.bins)

           #RCP time period 3:
            #Read in joint frequency table from part 6, an array of tables for all the model/runs available
            fractions.diff.all=readRDS(paste(dir.out, 'fractions.diff.all.',region.name,'.',season.name,'.',yr.start.real,'.',yr.end.real,'.',yr.start.rcp[3],'.',yr.end.rcp[3],'.rds',sep='')) #dim = n.bins-1 x n.bins-1 x length(x)
            fractions.mm.ens.med3=fractions.ensemble(fractions.diff.all, n.unique, x, n.bins)
            q.ens3=quadrant.ensemble(fractions.diff.all, n.unique, x, n.bins)

            mm.ens.season=abind(mm.ens.season, fractions.mm.ens.med1, fractions.mm.ens.med2, fractions.mm.ens.med3)
            q.ens.season=abind(q.ens.season, q.ens1, q.ens2, q.ens3)
      
    } #closes the season for-loop
    mm.ens.season=array(mm.ens.season, dim=c(n.bins-1, n.bins-1, 13)) # 13 = 4 seasons * 3 periods + 1 initializing matrix
    mm.ens.season=mm.ens.season[,,-1] #remove initializing matrix; dim=n.bins-1 x n.bins-1 x 12
    q.ens.season=array(q.ens.season, dim=c(4,2,13))
    q.ens.season=q.ens.season[,,-1]
    
    
    
#-------------------------------------------------------------------------------------
#Now read in the CMIP historical data 1979-2014 and run the ensemble; only 1 time period so no need to do tedious for loop 
#CMIP Historical:
#Read in joint frequency table from part 6, an array of tables for all the model/runs available
    fractions.all=readRDS(paste(dir.out,'fractions.all.',region.name,'.',season.names[1],'.',yr.start.hist,'.',yr.end.hist,'.rds',sep=''))
    fractions.hist.djf=fractions.ensemble(fractions.all, n.unique, x, n.bins)
            
    fractions.all=readRDS(paste(dir.out,'fractions.all.',region.name,'.',season.names[2],'.',yr.start.hist,'.',yr.end.hist,'.rds',sep=''))
    fractions.hist.mam=fractions.ensemble(fractions.all, n.unique, x, n.bins)

    fractions.all=readRDS(paste(dir.out,'fractions.all.',region.name,'.',season.names[3],'.',yr.start.hist,'.',yr.end.hist,'.rds',sep=''))
    fractions.hist.jja=fractions.ensemble(fractions.all, n.unique, x, n.bins)

    fractions.all=readRDS(paste(dir.out,'fractions.all.',region.name,'.',season.names[4],'.',yr.start.hist,'.',yr.end.hist,'.rds',sep=''))
    fractions.hist.son=fractions.ensemble(fractions.all, n.unique, x, n.bins)

#-------------------------------------------------------------------------------------
#Read in the observations joint frequency table; no need to ensemble since only one 'run'; used later in plotting
    fractions.obs.djf=readRDS(paste(dir.out,'fractions.',yr.start.obs,'.',yr.end.obs,'.',region.name,'.',season.names[1],'.rds',sep=''))
    fractions.obs.mam=readRDS(paste(dir.out,'fractions.',yr.start.obs,'.',yr.end.obs,'.',region.name,'.',season.names[2],'.rds',sep=''))
    fractions.obs.jja=readRDS(paste(dir.out,'fractions.',yr.start.obs,'.',yr.end.obs,'.',region.name,'.',season.names[3],'.rds',sep=''))
    fractions.obs.son=readRDS(paste(dir.out,'fractions.',yr.start.obs,'.',yr.end.obs,'.',region.name,'.',season.names[4],'.rds',sep=''))
    
#-------------------------------------------------------------------------------------
##### Get the grids needed to plot the region mask overlayed on blank globe #####
      lat...=readRDS("/scratch/users/tballard/plots/global.attributes/shum/lat") #used to make plot showing where study region is
      lon...=readRDS("/scratch/users/tballard/plots/global.attributes/shum/lon")
        fileName="/scratch/PI/omramom/reanalysis/ncep-doe-r2/4xdaily/ocean.mask/land.sfc.gauss.nc"
        land... = ncvar_get(nc_open(fileName), "land") #192 x 94
        land...[land...==0]=NA #set ocean values to NA instead of 0
        fileName="/scratch/PI/omramom/reanalysis/ncep-doe-r2/daily/shum/shum.2m.gauss.1979.nc"
        lon. = ncvar_get(nc_open(fileName), "lon") #192 values (0 to 358.12)
        lat. = ncvar_get(nc_open(fileName), "lat") #94 values (88.542 to -88.542)
        
          zone=regions[region,] #lon,lat coordinates, a 4-vector
          region.lat=lat.>zone[1] & lat.<zone[2] #True if lat is w/in bounds
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
          region.lat.mat=matrix(rep(region.lat,192), nrow=192, byrow=T)
          region.lon.mat=matrix(rep(region.lon,94), nrow=192)
          region.mask=region.lat.mat*region.lon.mat #192x94


   
#############################################################
###########           Plotting Time           ###############
#############################################################
    
    #plot.name=paste(dir.plot,"plot.joint.pdf.ensemble.diff.",region.names[region],'.',season.name,'.',yr.start.real,'.',yr.end.real,'.',yr.start.real2,'.',yr.end.real2,".pdf",sep="")
    plot.name=paste(dir.plot,"plot.joint.pdf.ensemble.diff.rcp85.",region.name,".pdf",sep="")
    colorbar=rev(colorRampPalette(brewer.pal(11,"RdBu"))(15))
    colorbar[8]="white"
    colorbar=c(rep(colorbar[1],4),colorbar,rep(colorbar[length(colorbar)],4))
    
    color1=rev(colorRampPalette(brewer.pal(9,"Purples"))(15))
    color2=colorRampPalette(brewer.pal(9,"Oranges"))(15)
    #colorbar2=c(color1,rep("white",1),color2)
    #colorbar2=c(rep(colorbar2[1],1),colorbar2,rep(colorbar2[length(colorbar2)],1))
    colorbar2=colorRampPalette(brewer.pal(11,"BrBG"))(15); colorbar2[8]="white"
    colorbar2=c(colorbar2[1], colorbar2, colorbar2[length(colorbar2)])
    
      
    colorbar.giorgi=c("deepskyblue","red")
    
    colorbar.hist=c("white",colorRampPalette(brewer.pal(9,"YlOrRd"))(15))

    legend.lim.djf=max(abs(quantile(as.vector(mm.ens.season[,,1:3]), c(.01,.99), na.rm=T))) #find limits for plotting consistent for all 3 time periods
    legend.lim.mam=max(abs(quantile(as.vector(mm.ens.season[,,4:6]), c(.01,.99), na.rm=T)))
    legend.lim.jja=max(abs(quantile(as.vector(mm.ens.season[,,7:9]), c(.01,.99), na.rm=T)))
    legend.lim.son=max(abs(quantile(as.vector(mm.ens.season[,,10:12]), c(.01,.99), na.rm=T)))
    legend.lim.historical.djf=quantile(c(as.vector(fractions.obs.djf),as.vector(fractions.hist.djf)), c(.99), na.rm=T) #get upper limit, lower limit is 0
    legend.lim.historical.mam=quantile(c(as.vector(fractions.obs.mam),as.vector(fractions.hist.mam)), c(.99), na.rm=T)
    legend.lim.historical.jja=quantile(c(as.vector(fractions.obs.jja),as.vector(fractions.hist.jja)), c(.99), na.rm=T)
    legend.lim.historical.son=quantile(c(as.vector(fractions.obs.son),as.vector(fractions.hist.son)), c(.99), na.rm=T)
      
    x.grid=seq(xlims[1], xlims[2], length=n.bins)
    y.grid=seq(ylims[1], ylims[2], length=n.bins)
    

    plot.joint.pdf=function(fractions.mm.ens.med, colorbar, x.grid, y.grid, yr.start.real2, yr.end.real2, 
                            yr.start.real, yr.end.real, n.bins, legend.lims, q){
      #First set any values exceeding or below the legend.limits equal to the max/min so they don't get whited out
      fractions.mm.ens.med[fractions.mm.ens.med>legend.lims]=legend.lims; fractions.mm.ens.med[fractions.mm.ens.med<(-1*legend.lims)]=(-1*legend.lims)
      image.plot(x.grid, y.grid, t(fractions.mm.ens.med)[,(n.bins-1):1], col=colorbar, las=1,
              legend.width=1.2, horizontal=T, bty='n', zlim=c(-1*legend.lims, legend.lims),
              xlab=expression(Temperature~(degree~C)), ylab="Humidity (g/kg)",
              main=paste("PDF Change (",yr.start.real2,'-',yr.end.real2,") - (", yr.start.real, "-", yr.end.real, ")", sep=""))
            abline(v=0, lty=2, col='grey30', lwd=2.5)
            abline(h=0, lty=2, col='grey30', lwd=2.5)
            legend("topright", paste(round(q[1,1], 3), " (", round(q[1,2]), "%)", sep=""), bty='n', cex=1.2)
            legend("topleft", paste(round(q[2,1], 3), " (", round(q[2,2]), "%)", sep=""), bty='n', cex=1.2) 
            legend("bottomleft", paste(round(q[3,1], 3), " (", round(q[3,2]), "%)", sep=""), bty='n', cex=1.2) 
            legend("bottomright", paste(round(q[4,1], 3), " (", round(q[4,2]), "%)", sep=""), bty='n', cex=1.2)
           
    } 
    
    plot.joint.pdf.historical=function(fractions.mm.ens.med, colorbar, x.grid, y.grid, yr.start, yr.end, 
                              n.bins, legend.lims, cmip){
      #First set any values exceeding or below the legend.limits equal to the max/min so they don't get whited out
      fractions.mm.ens.med[fractions.mm.ens.med>legend.lims]=legend.lims; fractions.mm.ens.med[fractions.mm.ens.med<(-1*legend.lims)]=(-1*legend.lims)
      if(cmip==TRUE){
        main=paste("CMIP5 Historical PDF (",yr.start,"-",yr.end,")",sep="")
        }
      else{
         main=paste("Reanalysis Historical PDF (",yr.start,"-",yr.end,")",sep="") 
        }
      image.plot(x.grid, y.grid, t(fractions.mm.ens.med)[,(n.bins-1):1], col=colorbar, las=1,
              legend.width=1.2, horizontal=T, bty='n', zlim=c(0, legend.lims),
              xlab=expression(Temperature~(degree~C)), ylab="Humidity (g/kg)",
              main=main)
            abline(v=0, lty=2, col='grey30', lwd=2.5)
            abline(h=0, lty=2, col='grey30', lwd=2.5)
    }

    
    pdf(plot.name, width=36, height=28)
    par(mfrow=c(4,6), mar=c(12.1,4.1,8.1,2.1)) #first mar number drops the legend as increase it, 3rd one increases space betwen rows
    ## Row 1: DJF
    suppressWarnings(quilt.plot(lon..., lat..., as.vector(region.mask*land...),nx = length(lon.), ny = length(lat.), 
           main="Study Region",col=colorbar.giorgi, bty='n',axes=F, add.legend=F, legend.shrink=.75))
         world(add = TRUE, col = "black",lwd=1.2)
        plot.joint.pdf.historical(fractions.obs.djf, colorbar=colorbar.hist, x.grid, y.grid, yr.start.obs, yr.end.obs, n.bins=n.bins, legend.lims=legend.lim.historical.djf, cmip=F)
        plot.joint.pdf.historical(fractions.hist.djf, colorbar=colorbar.hist, x.grid, y.grid, yr.start.real, yr.end.real, n.bins=n.bins, legend.lims=legend.lim.historical.djf, cmip=T)
        plot.joint.pdf(fractions.mm.ens.med = mm.ens.season[,,1], colorbar=colorbar2, x.grid, y.grid, yr.start.rcp[1], yr.end.rcp[1], yr.start.real, yr.end.real, n.bins=n.bins, legend.lims=legend.lim.djf, q=q.ens.season[,,1])
        plot.joint.pdf(fractions.mm.ens.med = mm.ens.season[,,2], colorbar=colorbar2, x.grid, y.grid, yr.start.rcp[2], yr.end.rcp[2], yr.start.real, yr.end.real, n.bins=n.bins, legend.lims=legend.lim.djf, q=q.ens.season[,,2])
        plot.joint.pdf(fractions.mm.ens.med = mm.ens.season[,,3], colorbar=colorbar2, x.grid, y.grid, yr.start.rcp[3], yr.end.rcp[3], yr.start.real, yr.end.real, n.bins=n.bins, legend.lims=legend.lim.djf, q=q.ens.season[,,3])
    ## Row 1: MAM
        plot(1, axes=F, xlab="",ylab="",col="white")
            text(1,1,paste("MAM"),cex=3.5)
        plot.joint.pdf.historical(fractions.obs.djf, colorbar=colorbar.hist, x.grid, y.grid, yr.start.obs, yr.end.obs, n.bins=n.bins, legend.lims=legend.lim.historical.djf, cmip=F)
        plot.joint.pdf.historical(fractions.hist.djf, colorbar=colorbar.hist, x.grid, y.grid, yr.start.real, yr.end.real, n.bins=n.bins, legend.lims=legend.lim.historical.djf, cmip=T)
        plot.joint.pdf(fractions.mm.ens.med = mm.ens.season[,,4], colorbar=colorbar2, x.grid, y.grid, yr.start.rcp[1], yr.end.rcp[1], yr.start.real, yr.end.real, n.bins=n.bins, legend.lims=legend.lim.mam, q=q.ens.season[,,4])
        plot.joint.pdf(fractions.mm.ens.med = mm.ens.season[,,5], colorbar=colorbar2, x.grid, y.grid, yr.start.rcp[2], yr.end.rcp[2], yr.start.real, yr.end.real, n.bins=n.bins, legend.lims=legend.lim.mam, q=q.ens.season[,,5])
        plot.joint.pdf(fractions.mm.ens.med = mm.ens.season[,,6], colorbar=colorbar2, x.grid, y.grid, yr.start.rcp[3], yr.end.rcp[3], yr.start.real, yr.end.real, n.bins=n.bins, legend.lims=legend.lim.mam, q=q.ens.season[,,6])
    ## Row 1: JJA
        plot(1, axes=F, xlab="",ylab="",col="white")
            text(1,1,paste("JJA"),cex=3.5)
        plot.joint.pdf.historical(fractions.obs.jja, colorbar=colorbar.hist, x.grid, y.grid, yr.start.obs, yr.end.obs, n.bins=n.bins, legend.lims=legend.lim.historical.jja, cmip=F)
        plot.joint.pdf.historical(fractions.hist.jja, colorbar=colorbar.hist, x.grid, y.grid, yr.start.real, yr.end.real, n.bins=n.bins, legend.lims=legend.lim.historical.jja, cmip=T)
        plot.joint.pdf(fractions.mm.ens.med = mm.ens.season[,,7], colorbar=colorbar2, x.grid, y.grid, yr.start.rcp[1], yr.end.rcp[1], yr.start.real, yr.end.real, n.bins=n.bins, legend.lims=legend.lim.jja, q=q.ens.season[,,7])
        plot.joint.pdf(fractions.mm.ens.med = mm.ens.season[,,8], colorbar=colorbar2, x.grid, y.grid, yr.start.rcp[2], yr.end.rcp[2], yr.start.real, yr.end.real, n.bins=n.bins, legend.lims=legend.lim.jja, q=q.ens.season[,,8])
        plot.joint.pdf(fractions.mm.ens.med = mm.ens.season[,,9], colorbar=colorbar2, x.grid, y.grid, yr.start.rcp[3], yr.end.rcp[3], yr.start.real, yr.end.real, n.bins=n.bins, legend.lims=legend.lim.jja, q=q.ens.season[,,9])
    ## Row 1: SON
        plot(1, axes=F, xlab="",ylab="",col="white")
            text(1,1,paste("SON"),cex=3.5)
        plot.joint.pdf.historical(fractions.obs.son, colorbar=colorbar.hist, x.grid, y.grid, yr.start.obs, yr.end.obs, n.bins=n.bins, legend.lims=legend.lim.historical.son, cmip=F)
        plot.joint.pdf.historical(fractions.hist.son, colorbar=colorbar.hist, x.grid, y.grid, yr.start.real, yr.end.real, n.bins=n.bins, legend.lims=legend.lim.historical.son, cmip=T)
        plot.joint.pdf(fractions.mm.ens.med = mm.ens.season[,,10], colorbar=colorbar2, x.grid, y.grid, yr.start.rcp[1], yr.end.rcp[1], yr.start.real, yr.end.real, n.bins=n.bins, legend.lims=legend.lim.son, q=q.ens.season[,,10])
        plot.joint.pdf(fractions.mm.ens.med = mm.ens.season[,,11], colorbar=colorbar2, x.grid, y.grid, yr.start.rcp[2], yr.end.rcp[2], yr.start.real, yr.end.real, n.bins=n.bins, legend.lims=legend.lim.son, q=q.ens.season[,,11])
        plot.joint.pdf(fractions.mm.ens.med = mm.ens.season[,,12], colorbar=colorbar2, x.grid, y.grid, yr.start.rcp[3], yr.end.rcp[3], yr.start.real, yr.end.real, n.bins=n.bins, legend.lims=legend.lim.son, q=q.ens.season[,,12])

    dev.off()


} #closes region for loop
    
    
    
    