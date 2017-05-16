suppressMessages(library(fields))
suppressMessages(library(ncdf4))
suppressMessages(library(RColorBrewer))
suppressMessages(library(abind))
suppressMessages(library(base))

yr.start=1979
yr.end=2010

#tmax.hw=readRDS("/scratch/users/tballard/shum/percentile.threshold/daily.data/tmax.hw") #192x94x365x36
shum.hw=readRDS("/scratch/users/tballard/shum/station.data/shum.hw")
dim.shum=dim(shum.hw) #365x32x89
#shum.hw.not=readRDS("/scratch/users/tballard/shum/percentile.threshold/daily.data/shum.hw.not")

## Define month indices ##
  jan=1:31; feb=32:59; mar=60:90; apr=91:120;may=121:151; jun=152:181;
  jul=182:212; aug=213:243; sep=244:273; oct=274:304; nov=305:334; dec=335:365;
  month.full=c("January","February","March","April","May","June","July","August","September","October","November","December")
  

##### Compute them OLS trends! #####
  #Function below computes OLS regression of shum vs year and outputs the slope value + SE
  #Then apply this function to every pixel using 'apply' command
  #Note the lm function automatically skips over NA's
  fit.lm=function(dataset, month){
    years=rep(c(yr.start:yr.end),each=length(month)) #1979, 1979, ... 2014, 2014 
    #shum=as.vector(dataset)
    a=tryCatch(summary(lm(dataset~years))$coefficient[2,1:2], error=function(e) c(NA,NA)) #slope for 'year' and SE
    return(a)
  }
 
  m=fit.lm(shum.hw[jan,,45],month=jan)
  
  #Array of lon,lat,month,results; results is 2D of the slope and its SE from running the regression
  #'aperm' rearranges order of arrays the aperm below switches the apply output from dim=2,192,94 to dim=192,94,2
  lm.trends=array(rep(NA,dim.shum[3]*12*2),c(dim.shum[3],12,2)) #initialize
  lm.trends[,1,]=aperm(apply(array(shum.hw[jan,,],c(length(jan)*dim.shum[2],dim.shum[3])), c(2), fit.lm, month=jan), c(2,1))   
  lm.trends[,2,]=aperm(apply(array(shum.hw[feb,,],c(length(feb)*dim.shum[2],dim.shum[3])), c(2), fit.lm, month=feb), c(2,1))   
  lm.trends[,3,]=aperm(apply(array(shum.hw[mar,,],c(length(mar)*dim.shum[2],dim.shum[3])), c(2), fit.lm, month=mar), c(2,1))   
  lm.trends[,4,]=aperm(apply(array(shum.hw[apr,,],c(length(apr)*dim.shum[2],dim.shum[3])), c(2), fit.lm, month=apr), c(2,1))   
  lm.trends[,5,]=aperm(apply(array(shum.hw[may,,],c(length(may)*dim.shum[2],dim.shum[3])), c(2), fit.lm, month=may), c(2,1))   
  lm.trends[,6,]=aperm(apply(array(shum.hw[jun,,],c(length(jun)*dim.shum[2],dim.shum[3])), c(2), fit.lm, month=jun), c(2,1))   
  lm.trends[,7,]=aperm(apply(array(shum.hw[jul,,],c(length(jul)*dim.shum[2],dim.shum[3])), c(2), fit.lm, month=jul), c(2,1))   
  lm.trends[,8,]=aperm(apply(array(shum.hw[aug,,],c(length(aug)*dim.shum[2],dim.shum[3])), c(2), fit.lm, month=aug), c(2,1))   
  lm.trends[,9,]=aperm(apply(array(shum.hw[sep,,],c(length(sep)*dim.shum[2],dim.shum[3])), c(2), fit.lm, month=sep), c(2,1))   
  lm.trends[,10,]=aperm(apply(array(shum.hw[oct,,],c(length(oct)*dim.shum[2],dim.shum[3])), c(2), fit.lm, month=oct), c(2,1))   
  lm.trends[,11,]=aperm(apply(array(shum.hw[nov,,],c(length(nov)*dim.shum[2],dim.shum[3])), c(2), fit.lm, month=nov), c(2,1))   
  lm.trends[,12,]=aperm(apply(array(shum.hw[dec,,],c(length(dec)*dim.shum[2],dim.shum[3])), c(2), fit.lm, month=dec), c(2,1))   
  
  
  saveRDS(lm.trends,"/scratch/users/tballard/shum/station.data/shum.hw.trend.rds")
  
  
  ##### Set non-significant trend values to NA ##### 
  w=lm.trends[,,1]+1.96*lm.trends[,,2] #96x12x2
  x=lm.trends[,,1]-1.96*lm.trends[,,2]
  
  lm.trends.adj=lm.trends[,,1] #
  lm.trends.adj[(w>0 & x<0)]=NA #set values where CI contains 0 to NA
  mean(is.na(lm.trends.adj)) #percent NA's (nonsignificant trends) ~70% (vs. 66ish for reanalysis)
 
  
#############################################################
###########           Plotting Time           ###############
#############################################################
  ##### Import lat and lon (made in clean.data.R) #####
  lat=readRDS("/scratch/users/tballard/shum/station.data/lat.rds")
  lon=readRDS("/scratch/users/tballard/shum/station.data/lon.rds")
  lon=lon-360
  
  dir="/scratch/users/tballard/plots/shum/station.data/" #where to save plots
  plot.name=paste(dir,"plot.shum.hw.monthly.trend.png",sep="")
  plot.width=2250*.9 #units are pixels
  plot.height=1500*.9*12
  plot.res=200*.9
  colorbar=rev(colorRampPalette(brewer.pal(11,"RdBu"))(20))
  colorbar=c(colorbar[1:10],"#FFFFFF",colorbar[11:20])
  legend.lim=c(-.05,.05)
  
  ##### Make the plot #####
  plot.slope=function(data,lon,lat,colorbar,legend.lim,month.full,title=""){
    quilt.plot(lon,lat,data,col=colorbar,zlim=legend.lim,nx=40,ny=40,
               main=title)
    US(add = TRUE, col = "black",lwd=1.2)
    par(font=2); legend("bottomright", month.full, bty='n', cex=2.5)
  }
  
  png(plot.name, units="px", width=plot.width, height=plot.height, res=plot.res)
    par(mfrow=c(12,1))
    plot.slope(lm.trends[,1,1],lon,lat,colorbar,legend.lim,month.full[1],title="OLS Trend in Humidity During Heat Events")
    plot.slope(lm.trends[,2,1],lon,lat,colorbar,legend.lim,month.full[2])
    plot.slope(lm.trends[,3,1],lon,lat,colorbar,legend.lim,month.full[3])
    plot.slope(lm.trends[,4,1],lon,lat,colorbar,legend.lim,month.full[4])
    plot.slope(lm.trends[,5,1],lon,lat,colorbar,legend.lim,month.full[5])
    plot.slope(lm.trends[,6,1],lon,lat,colorbar,legend.lim,month.full[6])
    plot.slope(lm.trends[,7,1],lon,lat,colorbar,legend.lim,month.full[7])
    plot.slope(lm.trends[,8,1],lon,lat,colorbar,legend.lim,month.full[8])
    plot.slope(lm.trends[,9,1],lon,lat,colorbar,legend.lim,month.full[9])
    plot.slope(lm.trends[,10,1],lon,lat,colorbar,legend.lim,month.full[10])
    plot.slope(lm.trends[,11,1],lon,lat,colorbar,legend.lim,month.full[11])
    plot.slope(lm.trends[,12,1],lon,lat,colorbar,legend.lim,month.full[12])
  
  dev.off()
  
  
  