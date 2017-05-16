suppressMessages(library(fields))
suppressMessages(library(ncdf4))
suppressMessages(library(RColorBrewer))
suppressMessages(library(abind))
suppressMessages()

yr.start=1979
yr.end=2010
n.years=yr.end-yr.start+1


#############################################################
###########          Read in the Data         ###############
############################################################# 
#the data has already been subsetted for the particular region of interest
station.data=readRDS('/scratch/users/tballard/shum/station.data/station.data.USOnly.rds') #23375x3x89
dates=readRDS('/scratch/users/tballard/shum/station.data/station.data.USOnly.dates.rds') #23375x89
dates=as.Date(dates[,1]) #take out the dates for a single station (same across all stations anyways)

shum=station.data[,2,]

month.full=c("January","February","March","April","May","June","July","August","September","October","November","December")
shum.median=c(rep(NA,dim(shum)[2]))
for (i in 1:12){
  month.index=months(dates)==month.full[i]
  shum.val=shum[month.index,]
  shum.med=apply(shum.val,2,median,na.rm=T)
  shum.median=cbind(shum.median,shum.med)
  #assign(paste('shum.median.',month.full[i],sep=""),shum.median)
  
}
shum.median=shum.median[,-1] #141stationsx12months of the median shum
saveRDS(shum.median, '/scratch/users/tballard/shum/station.data/shum.median.rds')








#############################################################
#######    Plotting time (On perconal computer)       #######
############################################################# 


setwd("~/Documents/shum/station.data")
lat=readRDS('lat.rds')
lon=readRDS('lon.rds')
shum.median=readRDS('shum.median.rds')

library(plotly)
Sys.setenv("plotly_username" = "333tristan")
Sys.setenv("plotly_api_key" = "h3bhvlgnr9")

# marker styling
m <- list(
  colorbar = list(title = "Humidity"),
  size = 8, opacity = 0.8, symbol = 'square',cmin=0,cmax=18#,showscale=F
)

# geo styling
g <- list(
  scope = 'usa',
  projection = list(type = 'albers usa'),
  showland = TRUE,
  landcolor = toRGB("gray95"),
  subunitcolor = toRGB("gray85"),
  countrycolor = toRGB("gray85"),
  countrywidth = 0.5,
  subunitwidth = 0.5
  
)
colorbar=(colorRampPalette(brewer.pal(11,"RdYlBu"))(20))
colorbar=c(colorbar,rep(colorbar[20],5))
#colorbar=rainbow(50)[1:45]
p=plot_ly(lat = lat, lon = lon, color = shum.median[,1],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Median Humidity 1979-2010<br>January', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.median.jan.png")
 
  p=plot_ly(lat = lat, lon = lon, color = shum.median[,2],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Median Humidity 1979-2010<br>February', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.median.feb.png")
  
  p=plot_ly(lat = lat, lon = lon, color = shum.median[,3],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Median Humidity 1979-2010<br>March', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.median.mar.png")
 
  
  p=plot_ly(lat = lat, lon = lon, color = shum.median[,4],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Median Humidity 1979-2010<br>April', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.median.apr.png")
 
  p=plot_ly(lat = lat, lon = lon, color = shum.median[,5],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Median Humidity 1979-2010<br>May', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.median.may.png")
 
  p=plot_ly(lat = lat, lon = lon, color = shum.median[,6],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Median Humidity 1979-2010<br>June', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.median.jun.png")
 
  p=plot_ly(lat = lat, lon = lon, color = shum.median[,7],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Median Humidity 1979-2010<br>July', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.median.jul.png")
 
  p=plot_ly(lat = lat, lon = lon, color = shum.median[,8],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Median Humidity 1979-2010<br>August', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.median.aug.png")
 
  p=plot_ly(lat = lat, lon = lon, color = shum.median[,9],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Median Humidity 1979-2010<br>September', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.median.sep.png")
 
  p=plot_ly(lat = lat, lon = lon, color = shum.median[,10],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Median Humidity 1979-2010<br>October', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.median.oct.png")
 
  p=plot_ly(lat = lat, lon = lon, color = shum.median[,11],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Median Humidity 1979-2010<br>November', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.median.nov.png")
 
  p=plot_ly(lat = lat, lon = lon, color = shum.median[,12],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Median Humidity 1979-2010<br>December', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.median.dec.png")
 
  

