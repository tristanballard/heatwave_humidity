suppressMessages(library(fields))
suppressMessages(library(ncdf4))
suppressMessages(library(RColorBrewer))
suppressMessages(library(abind))
suppressMessages(library(base))
suppressMessages(library(Rfit))
suppressMessages(library(MASS))
suppressMessages(library(gplots))


s=Sys.time()
for (m in 5:8){
#m=7 #month from 1:12
  
#region=2 #1=S.US , 2=N.US+Canada, 3=Europe
yr.start=1
yr.end=16
yr.start2=19
yr.end2=32



tmax.hw=readRDS("/scratch/users/tballard/shum/station.data/tmax.hw.SE.rds") #365x32x89
shum.hw=readRDS("/scratch/users/tballard/shum/station.data/shum.hw.SE.rds") #365x32x89
dim.shum=dim(shum.hw) #365x32x89

### Extract only the values where we have recordings for both shum and tmax ###
### b/c ~3% of the shum measurements are NA on the heat wave days, and a    ###
### very few days have shum measurements but no tmax                        ###
tmax.hw.adj=tmax.hw
tmax.hw.adj[is.na(shum.hw)]=NA
shum.hw.adj=shum.hw
shum.hw.adj[is.na(tmax.hw)]=NA

#Check that the number of NA's are now equal:
if (sum(is.na(tmax.hw.adj))==sum(is.na(shum.hw.adj))){
  print("Good Job")
  } else {
    print("Ya done goofed!")
  }


## Define month indices ##
  jan=1:31; feb=32:59; mar=60:90; apr=91:120;may=121:151; jun=152:181;
  jul=182:212; aug=213:243; sep=244:273; oct=274:304; nov=305:334; dec=335:365;
month.index=list(jan,feb,mar,apr,may,jun,jul,aug,sep,oct,nov,dec)
## Extract the shum/tmax values for the month of interest. If I were better at matrix operations
## I might be able to avoid this, but running regression on 192x94x365x36 is tricky when
## you want 3 dimensions instead

  tmax.p1=tmax.hw.adj[month.index[[m]],yr.start:yr.end,]
  shum.p1=shum.hw.adj[month.index[[m]],yr.start:yr.end,]

################################################
#######  Repeat for the 2nd Time period  #######
################################################
  tmax.p2=tmax.hw.adj[month.index[[m]],yr.start2:yr.end2,]
  shum.p2=shum.hw.adj[month.index[[m]],yr.start2:yr.end2,]
  

##########################################################  
##### Now that we have our monthly values              ###  
##### extract only the non-NA's and non-0's (0s happen ###
##### over grid cells if there wasn't a heat wave)     ###
##########################################################
##Left off here, for some reason tmax.val.1 and shum.val.1 aren't equal length....no zero values in tmax.p1 or shum.p1, but diff # of NA's
tmax.val.1=as.vector(tmax.p1); 
tmax.val.1=tmax.val.1[!is.na(tmax.val.1)]
tmax.val.1=tmax.val.1[tmax.val.1!=0]

shum.val.1=as.vector(shum.p1); 
shum.val.1=shum.val.1[!is.na(shum.val.1)]
shum.val.1=shum.val.1[shum.val.1!=0]

tmax.val.2=as.vector(tmax.p2); 
tmax.val.2=tmax.val.2[!is.na(tmax.val.2)]
tmax.val.2=tmax.val.2[tmax.val.2!=0]

shum.val.2=as.vector(shum.p2); 
shum.val.2=shum.val.2[!is.na(shum.val.2)]
shum.val.2=shum.val.2[shum.val.2!=0]

## Now compute the # of hw events in each period
n.events.1=length(tmax.val.1) #length is the same for shum.val.1
n.events.2=length(tmax.val.2)

### Now compute the empirical bivariate density 
f1 = kde2d(tmax.val.1, shum.val.1, n = 50, lims = c(0, 55, 0, 45)) #bounds in plot; n controls resolution
f2 = kde2d(tmax.val.2, shum.val.2, n = 50, lims = c(0, 55, 0, 45)) #bounds in plot; n controls resolution

### Now subtract two bivariate densities and plot
f3=f2$z-f1$z #the 'z' component of the object is the probabilities


### Now compute the empirical bivariate density but w/ histogram-like binning instead
df=data.frame(tmax.val.1,shum.val.1)    
nbins <- 30 #50
x.bin <- seq(floor(min(tmax.val.1,tmax.val.2)), ceiling(max(tmax.val.1,tmax.val.2)), length=nbins)
y.bin <- seq(floor(min(shum.val.1,shum.val.2)), ceiling(max(shum.val.1,shum.val.2)), length=nbins)

freq <-  as.data.frame(table(findInterval(df[,1], x.bin),findInterval(df[,2], y.bin)))
freq[,1] <- as.numeric(freq[,1])
freq[,2] <- as.numeric(freq[,2])

freq2D.1 <- diag(nbins)*0
freq2D.1[cbind(freq[,1], freq[,2])] <- freq[,3]

df=data.frame(tmax.val.2,shum.val.2)    
freq <-  as.data.frame(table(findInterval(df[,1], x.bin),findInterval(df[,2], y.bin)))
freq[,1] <- as.numeric(freq[,1])
freq[,2] <- as.numeric(freq[,2])

freq2D.2 <- diag(nbins)*0
freq2D.2[cbind(freq[,1], freq[,2])] <- freq[,3]

##Now subtract the two matrices of observed counts/frequencies
freq2D.3=freq2D.2-freq2D.1

#############################################################
###########           Plotting Time           ###############
#############################################################
 

#In order to make the plot showing the region we're zoned in on, read in some long stuff from the original
#plot.joint.pdf.R w/ reanalysis data
  zone=c(30,44,255,290) #S. US
  ##### Read in lat/lon values you'll use to make the mask #####
  fileName="/scratch/PI/omramom/reanalysis/ncep-doe-r2/daily/shum/shum.2m.gauss.1979.nc"
  lon. = ncvar_get(nc_open(fileName), "lon") #192 values (0 to 358.12)
  lat. = ncvar_get(nc_open(fileName), "lat") #94 values (88.542 to -88.542)
    region.lat=lat.>zone[1] & lat.<zone[2] #True if lat is w/in bounds
    region.lon=lon.>zone[3] & lon.<zone[4] #True if lon is w/in bounds
    region.lat[region.lat==FALSE]=NA #Set False values to NA
    region.lon[region.lon==FALSE]=NA
    region.lat=region.lat+0; region.lon=region.lon+0 #Set TRUE values to 1
    #Now that you have a vector for lat and vector for lon of NA or 1's,
    #Combine that into a matrix format that you can multiply things by
    region.lat.mat=matrix(rep(region.lat,192), nrow=192, byrow=T)
    region.lon.mat=matrix(rep(region.lon,94), nrow=192)
    region.mask=region.lat.mat*region.lon.mat #192x94
  
    fileName="/scratch/PI/omramom/reanalysis/ncep-doe-r2/4xdaily/ocean.mask/land.sfc.gauss.nc"
    land = ncvar_get(nc_open(fileName), "land") #192 x 94
    land[land==0]=NA #set ocean values to NA instead of 0
    lat=readRDS("/scratch/users/tballard/plots/global.attributes/shum/lat") 
    lon=readRDS("/scratch/users/tballard/plots/global.attributes/shum/lon")


### Now we can plot

month.full=c("January","February","March","April","May","June","July","August","September","October","November","December")
dir="/scratch/users/tballard/plots/shum/station.data/" #where to save plots
    plot.name=paste(dir,"plot.joint.pdf.station.data.US.",month.full[m],".South.png",sep="")
    plot.width=3300*.9 #1.3 #units are pixels
    plot.height=700*1.3*1.6
    plot.res=200*1.1
    colorbar=rev(colorRampPalette(brewer.pal(11,"RdBu"))(15))
    colorbar[8]="white"
    colorbar=c(rep(colorbar[1],4),colorbar,rep(colorbar[length(colorbar)],4))
   # colorbar.giorgi=c("yellowgreen","red")
    colorbar.giorgi=c("darkorange2","red")
   # colorbar.giorgi=c("deepskyblue","red")  
    colorbar.2=c("white",colorRampPalette(brewer.pal(9,"YlOrRd"))(15))
    colorbar.2=c(rep(colorbar.2[1],1),colorbar.2,rep(colorbar.2[length(colorbar.2)],3))
    legend.lim=c(-.007,.007)
    legend.lim.2=c(0,.022)
    legend.lim.3=c(0,75)
    legend.lim.4=c(-45,45)
    xlims=c(0,45)
    ylims=c(0,25)
    
png(plot.name, units="px", width=plot.width, height=plot.height, res=plot.res)
  par(mfrow=c(2,4), mar=c(6.1,4.1,4.1,2.1)) #originally 5.1,4.2,4.1,2.1; adds a bit more vertical space
    quilt.plot(lon, lat, as.vector(region.mask*land),nx = 192, ny = 94, 
           main="Study Region",col=colorbar.giorgi, bty='n',axes=F, add.legend=F, legend.shrink=.75)
        world(add = TRUE, col = "black",lwd=1.2)
    image.plot(f1$x,f1$y,f1$z, col=colorbar.2, 
          legend.width=1, horizontal=T, xlab="Temperature",
          ylab="Humidity", las=1, bty='n', zlim=legend.lim.2, xlim=xlims, ylim=ylims,
          main=paste(month.full[m],yr.start+1978,"-",yr.end+1978,"PDF"), legend.shrink=.75)
          text(6,2,paste("n=",n.events.1,sep=""),cex=1.5)
    image.plot(f1$x,f1$y,f2$z, col=colorbar.2, 
          legend.width=1, horizontal=T, xlab="Temperature",
          ylab="Humidity", las=1, bty='n', zlim=legend.lim.2, xlim=xlims, ylim=ylims,
          main=paste(month.full[m],yr.start2+1978,"-",yr.end2+1978,"PDF"), legend.shrink=.75)
          text(6,2,paste("n=",n.events.2,sep=""),cex=1.5)
    image.plot(f1$x,f1$y,f3, col=colorbar, 
          legend.width=1, horizontal=T, xlab="Temperature",
          ylab="Humidity", las=1, bty='n', zlim=legend.lim, xlim=xlims, ylim=ylims,
          main=paste(month.full[m],"PDF Change"), legend.shrink=.75)
    
    ##Row 2 plots
    plot(1, axes=F, xlab="",ylab="",col="white")
    text(1,1,month.full[m],cex=4.5)
    image.plot(x.bin, y.bin, freq2D.1, col=colorbar.2, horizontal=T,
               ylab="Humidity", las=1, bty='n',xlab="Temperature", xlim=xlims, ylim=ylims,
               main=paste(month.full[m],"Empirical PDF"), zlim=legend.lim.3, legend.shrink=.75)
    image.plot(x.bin, y.bin, freq2D.2, col=colorbar.2, horizontal=T,
               ylab="Humidity", las=1, bty='n',xlab="Temperature", xlim=xlims, ylim=ylims,
               main=paste(month.full[m],"Empirical PDF"), zlim=legend.lim.3, legend.shrink=.75)
    image.plot(x.bin, y.bin, freq2D.3, col=colorbar, horizontal=T,
               ylab="Humidity", las=1, bty='n',xlab="Temperature", xlim=xlims, ylim=ylims,
               main=paste(month.full[m],"Empirical PDF Change"), zlim=legend.lim.4, legend.shrink=.75)
dev.off()





#### Marginal histograms + scatterplot 
    plot.name=paste(dir,"plot.joint.pdf.station.data.US.",month.full[m],".South.marginal.hist1.png",sep="")
    plot.width=3300*.9 #1.3 #units are pixels
    plot.height=700*1.3*1.6
    plot.res=200*1.1
    
png(plot.name, units="px", width=plot.width, height=plot.height, res=plot.res)
    scatterhist = function(x, y, xlab="", ylab=""){
      zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
      #layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
      layout(zones, c(3,1),c(1,3))
      xhist = hist(x, plot=FALSE)
      yhist = hist(y, plot=FALSE)
      top = max(c(xhist$counts, yhist$counts))
      par(mar=c(3,3,1,1))
      plot(x,y)
      #par(mar=c(0,3,1,1))
      par(mar=c(0,2,1,0))
      barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
      #par(mar=c(3,0,1,1))
      par(mar=c(c,0,.5,1))
      barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
      par(oma=c(3,3,0,0))
    }
    scatterhist(tmax.val.1,shum.val.1)
dev.off()

    plot.name=paste(dir,"plot.joint.pdf.station.data.US.",month.full[m],".South.marginal.hist2.png",sep="")
    png(plot.name, units="px", width=plot.width, height=plot.height, res=plot.res)
    scatterhist(tmax.val.2,shum.val.2)
    
dev.off()


 
}
t=Sys.time()-s


