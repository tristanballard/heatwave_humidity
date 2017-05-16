suppressMessages(library(fields))
suppressMessages(library(ncdf4))
suppressMessages(library(RColorBrewer))
suppressMessages(library(abind))
suppressMessages(library(base))
suppressMessages(library(Rfit))
suppressMessages(library(MASS))
suppressMessages(library(gplots))

#62GB, 36min
dir.out="/scratch/users/tballard/shum/giorgi.regions/output/"
dir.plot="/scratch/users/tballard/shum/giorgi.regions/reanalysis/plots/"

#-------------------------------------------------------------------------------------
regions=readRDS("/scratch/users/tballard/shum/giorgi.regions/giorgi.regions.rds") #25 x 4
region.names=rownames(regions)

yr.start=1
yr.end=18
yr.start2=19
yr.end2=36
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
    
    
    
    
    ################################################
    #######  Repeat for the 2nd Time period  #######
    ################################################
    daily.shum=shum.hw[,,season,yr.start2:yr.end2] 
    dim(daily.shum) = c(dim(daily.shum)[1], dim(daily.shum)[2], dim(daily.shum)[3]*dim(daily.shum)[4]) #collapses from a 4dim array to 3dim array
    daily.tmax=tmax.hw[,,season,yr.start2:yr.end2] 
    dim(daily.tmax) = c(dim(daily.tmax)[1], dim(daily.tmax)[2], dim(daily.tmax)[3]*dim(daily.tmax)[4]) #collapses from a 4dim array to 3dim array
    
    ##### Read in ocean/land mask file 1's are land #####
    mask=function(data,land){
      new.data=land*data
      return(new.data)
    }
    tmax.mask=apply(daily.tmax, c(3), mask, land=land)
    tmax.mask=array(tmax.mask, dim=dim(daily.tmax))
    shum.mask=apply(daily.shum, c(3), mask, land=land)
    shum.mask=array(shum.mask, dim=dim(daily.shum))
    
    
    ##### Mask out as well all the globe not in the US and southern Canada #####
    ####Below is applying the mask to your entire array of values
    mask=function(data, region.mask){
      new.data=region.mask*data
      return(new.data)
    }
    tmax.mask.p2=apply(tmax.mask, c(3), mask, region.mask=region.mask)
    tmax.mask.p2=array(tmax.mask.p2, dim=dim(daily.tmax))
    shum.mask.p2=apply(shum.mask, c(3), mask, region.mask=region.mask)
    shum.mask.p2=array(shum.mask.p2, dim=dim(daily.shum))
    
    
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
    #tmax.val.1=tmax.val.1[tmax.val.1!=0] #this step is removing some values, I think b/c median is the actual observed value
    
    shum.med=apply(shum.mask.p1, c(1,2), median, na.rm=T) #360 x 180
    shum.med[is.na(shum.med)]=0
    shum.val.1=sweep(shum.mask.p1, c(1,2), shum.med, FUN="-")

    shum.val.1=as.vector(shum.val.1); 
    shum.val.1=shum.val.1[!is.na(shum.val.1)]
    #shum.val.1=shum.val.1[shum.val.1!=0]
    
    tmax.val.2=sweep(tmax.mask.p2, c(1,2), tmax.med, FUN="-") #subtract the pixel by pixel medians from the first time period
    tmax.val.2=as.vector(tmax.val.2);
    tmax.val.2=tmax.val.2[!is.na(tmax.val.2)]
    #tmax.val.2=tmax.val.2[tmax.val.2!=0]
    
    shum.val.2=sweep(shum.mask.p2, c(1,2), shum.med, FUN="-") #subtract the pixel by pixel medians from the first time period
    shum.val.2=as.vector(shum.val.2); 
    shum.val.2=shum.val.2[!is.na(shum.val.2)]
    #shum.val.2=shum.val.2[shum.val.2!=0]
    
  
     saveRDS(tmax.val.1, paste(dir.out,'tmax.val.1.',region.names[region],'.',season.name,'.rds',sep=''))
     saveRDS(shum.val.1, paste(dir.out,'shum.val.1.',region.names[region],'.',season.name,'.rds',sep=''))
     saveRDS(tmax.val.2, paste(dir.out,'tmax.val.2.',region.names[region],'.',season.name,'.rds',sep=''))
     saveRDS(tmax.val.2, paste(dir.out,'shum.val.2.',region.names[region],'.',season.name,'.rds',sep=''))
  
    
    
    ## Now compute the # of hw events in each period
    n.events.1=length(tmax.val.1) #length is the same for shum.val.1
    n.events.2=length(tmax.val.2)
    saveRDS(n.events.1, paste(dir.out,'n.events.1.',region.names[region],'.',season.name,'.rds',sep=''))
    saveRDS(n.events.2, paste(dir.out,'n.events.2.',region.names[region],'.',season.name,'.rds',sep=''))
    
    

    
    # ### Now compute the empirical bivariate density 
    # #f1 = kde2d(tmax.val.1, shum.val.1, n = 100, lims = c(-25, 25, -15, 15)) #bounds in plot; n controls resolution
    # #f2 = kde2d(tmax.val.2, shum.val.2, n = 100, lims = c(-25, 25, -15, 15)) #bounds in plot; n controls resolution
    # 
    # f1 = kde2d(tmax.val.1, shum.val.1, n = 100, lims = c(-20, 20, -12, 12)) #bounds in plot; n controls resolution
    # f2 = kde2d(tmax.val.2, shum.val.2, n = 100, lims = c(-20, 20, -12, 12)) #bounds in plot; n controls resolution
    # 
    # saveRDS(f1, paste(dir.out, 'anomaly.f1.',region.names[region],'.',season.name,'.rds',sep=''))
    # saveRDS(f2, paste(dir.out, 'anomaly.f2.',region.names[region],'.',season.name,'.rds',sep=''))
    # 
    # ### Now subtract two bivariate densities and plot
    # f3=f2$z-f1$z #the 'z' component of the object is the probabilities
    # 
    # 
    # ### Now compute the empirical bivariate density but w/ histogram-like binning instead
    # df=data.frame(tmax.val.1,shum.val.1)    
    # nbins <- 50
    # x.bin <- seq(floor(min(tmax.val.1,tmax.val.2)), ceiling(max(tmax.val.1,tmax.val.2)), length=nbins)
    # y.bin <- seq(floor(min(shum.val.1,shum.val.2)), ceiling(max(shum.val.1,shum.val.2)), length=nbins)
    # 
    # freq <-  as.data.frame(table(findInterval(df[,1], x.bin),findInterval(df[,2], y.bin)))
    # freq[,1] <- as.numeric(freq[,1])
    # freq[,2] <- as.numeric(freq[,2])
    # 
    # freq2D.1 <- diag(nbins)*0
    # freq2D.1[cbind(freq[,1], freq[,2])] <- freq[,3]
    # 
    # df=data.frame(tmax.val.2,shum.val.2)    
    # freq <-  as.data.frame(table(findInterval(df[,1], x.bin),findInterval(df[,2], y.bin)))
    # freq[,1] <- as.numeric(freq[,1])
    # freq[,2] <- as.numeric(freq[,2])
    # 
    # freq2D.2 <- diag(nbins)*0
    # freq2D.2[cbind(freq[,1], freq[,2])] <- freq[,3]
    # 
    # ##Now subtract the two matrices of observed counts/frequencies
    # freq2D.3=freq2D.2-freq2D.1
    # 
    # #############################################################
    # ###########           Plotting Time           ###############
    # #############################################################
    # plot.name=paste(dir.plot,'plot.joint.pdf.',region.names[region],'.',season.name,'.png',sep="")
    # 
    # plot.width=3300*.9 #1.3 #units are pixels
    # plot.height=700*1.3*1.6
    # plot.res=200*1.1
    # colorbar=rev(colorRampPalette(brewer.pal(11,"RdBu"))(15))
    # colorbar[8]="white"
    # colorbar=c(rep(colorbar[1],4),colorbar,rep(colorbar[length(colorbar)],4))
    # # colorbar.region=c("yellowgreen","red")
    # colorbar.region=c("darkorange2","red")
    # # colorbar.region=c("deepskyblue","red")  
    # colorbar.2=c("white",colorRampPalette(brewer.pal(9,"YlOrRd"))(15))
    # colorbar.2=c(rep(colorbar.2[1],1),colorbar.2,rep(colorbar.2[length(colorbar.2)],3))
    # #legend.1=list(c(-.0055,.0055),c(-.018,.018),c(-.024,.024))
    # #legend.lim=unlist(legend.1[region])
    # legend.lim=c(-.01,.01) 
    # #legend.2=list(c(0,.015),c(0,.035),c(0,.05))
    # #legend.lim.2=unlist(legend.2[region])
    # legend.lim.2=c(0,.025)
    # #legend.3=list(c(0,130),c(0,30),c(0,20))
    # #legend.lim.3=unlist(legend.3[region])
    # legend.lim.3=c(0,40)
    # #legend.4=list(c(-65,65),c(-20,20),c(-10,10))
    # #legend.lim.4=unlist(legend.4[region])
    # legend.lim.4=c(-20,20)
    # #xlims=c(-25,25)
    # #ylims=c(-15,15)
    # xlims=c(-20,20)
    # ylims=c(-12,12)
    # 
    # png(plot.name, units="px", width=plot.width, height=plot.height, res=plot.res)
    # par(mfrow=c(2,4), mar=c(6.1,4.1,4.1,2.1)) #originally 5.1,4.2,4.1,2.1; adds a bit more vertical space
    # quilt.plot(lon, lat, as.vector(region.mask*land),nx = n.lon, ny = n.lat, 
    #            main=paste("Study Region",region.name),col=colorbar.region, bty='n',axes=F, add.legend=F, legend.shrink=.75)
    # world(add = TRUE, col = "black",lwd=1.2)
    # image.plot(f1$x,f1$y,f1$z, col=colorbar.2, 
    #            legend.width=1, horizontal=T, xlab="Temperature",
    #            ylab="Humidity", las=1, bty='n', zlim=legend.lim.2, xlim=xlims, ylim=ylims,
    #            main=paste(yr.start+1978,"-",yr.end+1978,"PDF"), legend.shrink=.75)
    # abline(v=0,lty=3,lwd=.5); abline(h=0,lty=3,lwd=.5)
    # text(6,2,paste("n=",n.events.1,sep=""),cex=1.5)
    # image.plot(f1$x,f1$y,f2$z, col=colorbar.2, 
    #            legend.width=1, horizontal=T, xlab="Temperature",
    #            ylab="Humidity", las=1, bty='n', zlim=legend.lim.2, xlim=xlims, ylim=ylims,
    #            main=paste(yr.start2+1978,"-",yr.end2+1978,"PDF"), legend.shrink=.75)
    # abline(v=0,lty=3,lwd=.5); abline(h=0,lty=3,lwd=.5)
    # text(6,2,paste("n=",n.events.2,sep=""),cex=1.5)
    # image.plot(f1$x,f1$y,f3, col=colorbar, 
    #            legend.width=1, horizontal=T, xlab="Temperature",
    #            ylab="Humidity", las=1, bty='n', zlim=legend.lim, xlim=xlims, ylim=ylims,
    #            main=paste("PDF Change"), legend.shrink=.75)
    # abline(v=0,lty=3,lwd=.5); abline(h=0,lty=3,lwd=.5)
    # 
    # ##Row 2 plots
    # plot(1, axes=F, xlab="",ylab="",col="white")
    # text(1,1,season.name,cex=2.5)
    # image.plot(x.bin, y.bin, freq2D.1, col=colorbar.2, horizontal=T,
    #            ylab="Humidity", las=1, bty='n',xlab="Temperature", xlim=xlims, ylim=ylims,
    #            main=paste("Empirical PDF"), zlim=legend.lim.3, legend.shrink=.75)
    # abline(v=0,lty=3,lwd=.5); abline(h=0,lty=3,lwd=.5)
    # image.plot(x.bin, y.bin, freq2D.2, col=colorbar.2, horizontal=T,
    #            ylab="Humidity", las=1, bty='n',xlab="Temperature", xlim=xlims, ylim=ylims,
    #            main=paste("Empirical PDF"), zlim=legend.lim.3, legend.shrink=.75)
    # abline(v=0,lty=3,lwd=.5); abline(h=0,lty=3,lwd=.5)
    # image.plot(x.bin, y.bin, freq2D.3, col=colorbar, horizontal=T,
    #            ylab="Humidity", las=1, bty='n',xlab="Temperature", xlim=xlims, ylim=ylims,
    #            main=paste("Empirical PDF Change"), zlim=legend.lim.4, legend.shrink=.75)
    # abline(v=0,lty=3,lwd=.5); abline(h=0,lty=3,lwd=.5)
    # dev.off()
    
    
    # 
    # #############################################################
    # ###########           Plotting Time           ###############
    # #############################################################
    # plot.name=paste(dir.plot,'plot.joint.pdf.',region.names[region],'.',season.name,'.png',sep="")
    # 
    # plot.width=3300*1.4 #1.3 #units are pixels
    # plot.height=700*1.3*1.6*2.3
    # plot.res=200*1.1
    # colorbar=rev(colorRampPalette(brewer.pal(11,"RdBu"))(15))
    # colorbar[8]="white"
    # colorbar=c(rep(colorbar[1],4),colorbar,rep(colorbar[length(colorbar)],4))
    # # colorbar.region=c("yellowgreen","red")
    # colorbar.region=c("darkorange2","red")
    # # colorbar.region=c("deepskyblue","red")  
    # colorbar.2=c("white",colorRampPalette(brewer.pal(9,"YlOrRd"))(15))
    # colorbar.2=c(rep(colorbar.2[1],1),colorbar.2,rep(colorbar.2[length(colorbar.2)],3))
    # #legend.1=list(c(-.0055,.0055),c(-.018,.018),c(-.024,.024))
    # #legend.lim=unlist(legend.1[region])
    # legend.lim=c(-.012,.012) 
    # #legend.2=list(c(0,.015),c(0,.035),c(0,.05))
    # #legend.lim.2=unlist(legend.2[region])
    # legend.lim.2=c(0,.03)
    # #legend.3=list(c(0,130),c(0,30),c(0,20))
    # #legend.lim.3=unlist(legend.3[region])
    # legend.lim.3=c(0,60)
    # #legend.4=list(c(-65,65),c(-20,20),c(-10,10))
    # #legend.lim.4=unlist(legend.4[region])
    # legend.lim.4=c(-32,32)
    # #xlims=c(-25,25)
    # #ylims=c(-15,15)
    # xlims=c(-20,20)
    # ylims=c(-12,12)
    # 
    # png(plot.name, units="px", width=plot.width, height=plot.height, res=plot.res)
    # par(mfrow=c(3,4), mar=c(6.1,4.1,4.1,2.1)) #originally 5.1,4.2,4.1,2.1; adds a bit more vertical space
    # quilt.plot(lon, lat, as.vector(region.mask*land),nx = n.lon, ny = n.lat, 
    #            main=paste("Study Region",region.name),col=colorbar.region, bty='n',axes=F, add.legend=F, legend.shrink=.75)
    # world(add = TRUE, col = "black",lwd=1.2)
    # image.plot(f1$x,f1$y,f1$z, col=colorbar.2, 
    #            legend.width=1, horizontal=T, xlab="Temperature",
    #            ylab="Humidity", las=1, bty='n', zlim=legend.lim.2, xlim=xlims, ylim=ylims,
    #            main=paste(yr.start+1978,"-",yr.end+1978,"PDF"), legend.shrink=.75)
    # abline(v=0,lty=3,lwd=.5); abline(h=0,lty=3,lwd=.5)
    # text(6,2,paste("n=",n.events.1,sep=""),cex=1.5)
    # image.plot(f1$x,f1$y,f2$z, col=colorbar.2, 
    #            legend.width=1, horizontal=T, xlab="Temperature",
    #            ylab="Humidity", las=1, bty='n', zlim=legend.lim.2, xlim=xlims, ylim=ylims,
    #            main=paste(yr.start2+1978,"-",yr.end2+1978,"PDF"), legend.shrink=.75)
    # abline(v=0,lty=3,lwd=.5); abline(h=0,lty=3,lwd=.5)
    # text(6,2,paste("n=",n.events.2,sep=""),cex=1.5)
    # image.plot(f1$x,f1$y,f3, col=colorbar, 
    #            legend.width=1, horizontal=T, xlab="Temperature",
    #            ylab="Humidity", las=1, bty='n', zlim=legend.lim, xlim=xlims, ylim=ylims,
    #            main=paste("PDF Change"), legend.shrink=.75)
    # abline(v=0,lty=3,lwd=.5); abline(h=0,lty=3,lwd=.5)
    # 
    # ##Row 2 plots
    # plot(1, axes=F, xlab="",ylab="",col="white")
    # text(1,1,season.name,cex=2.5)
    # image.plot(x.bin, y.bin, freq2D.1, col=colorbar.2, horizontal=T,
    #            ylab="Humidity", las=1, bty='n',xlab="Temperature", xlim=xlims, ylim=ylims,
    #            main=paste("Empirical PDF"), zlim=legend.lim.3, legend.shrink=.75)
    # abline(v=0,lty=3,lwd=.5); abline(h=0,lty=3,lwd=.5)
    # image.plot(x.bin, y.bin, freq2D.2, col=colorbar.2, horizontal=T,
    #            ylab="Humidity", las=1, bty='n',xlab="Temperature", xlim=xlims, ylim=ylims,
    #            main=paste("Empirical PDF"), zlim=legend.lim.3, legend.shrink=.75)
    # abline(v=0,lty=3,lwd=.5); abline(h=0,lty=3,lwd=.5)
    # image.plot(x.bin, y.bin, freq2D.3, col=colorbar, horizontal=T,
    #            ylab="Humidity", las=1, bty='n',xlab="Temperature", xlim=xlims, ylim=ylims,
    #            main=paste("Empirical PDF Change"), zlim=legend.lim.4, legend.shrink=.75)
    # abline(v=0,lty=3,lwd=.5); abline(h=0,lty=3,lwd=.5)
    # 
    # 
    # ##Row 3 plots
    # plot(1, axes=F, xlab="",ylab="",col="white")
    # image.plot(x.bin, y.bin, (freq2D.1/n.events.1)*100, col=colorbar.2, horizontal=T,
    #            ylab="Humidity", las=1, bty='n',xlab="Temperature", xlim=xlims, ylim=ylims,
    #            main=paste("Empirical PDF"), zlim=legend.lim.3/75, legend.shrink=.75)
    # abline(v=0,lty=3,lwd=.5); abline(h=0,lty=3,lwd=.5)
    # image.plot(x.bin, y.bin, (freq2D.2/n.events.2)*100, col=colorbar.2, horizontal=T,
    #            ylab="Humidity", las=1, bty='n',xlab="Temperature", xlim=xlims, ylim=ylims,
    #            main=paste("Empirical PDF"), zlim=legend.lim.3/75, legend.shrink=.75)
    # abline(v=0,lty=3,lwd=.5); abline(h=0,lty=3,lwd=.5)
    # freq3=(freq2D.2/n.events.2)*100-(freq2D.1/n.events.1)*100
    # image.plot(x.bin, y.bin, freq3, col=colorbar, horizontal=T,
    #            ylab="Humidity", las=1, bty='n',xlab="Temperature", xlim=xlims, ylim=ylims,
    #            main=paste("Empirical PDF Change"), zlim=legend.lim.4/75, legend.shrink=.75)
    # abline(v=0,lty=3,lwd=.5); abline(h=0,lty=3,lwd=.5)
    # dev.off()
    # 
  }    
}
