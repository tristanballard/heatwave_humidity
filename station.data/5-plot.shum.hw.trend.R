setwd("~/Documents/shum/station.data")
lat=readRDS('lat.SE.rds')
lon=readRDS('lon.SE.rds')
lm.trends=readRDS('shum.hw.trend.rds') #make sure this is for the whole US not just the SE; 141x12x2

## Set non-significant trends to 0 (plotly doesn't plot NA's well), optional to plot this data or raw lm.trends data, both done below
  w=lm.trends[,,1]+1.96*lm.trends[,,2] #141x12
  x=lm.trends[,,1]-1.96*lm.trends[,,2]
  
  lm.trends.adj=lm.trends[,,1] #
  lm.trends.adj[(w>0 & x<0)]=0 #141x12




library(plotly)
Sys.setenv("plotly_username" = "333tristan")
Sys.setenv("plotly_api_key" = "h3bhvlgnr9")

# marker styling
m <- list(
  colorbar = list(title = "Trend in Humidity During Heat Events"),
  size = 8, opacity = 0.8, symbol = 'square',cmin=-.25,cmax=.25#,showscale=F
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
colorbar=rev(colorRampPalette(brewer.pal(11,"RdBu"))(20))
  colorbar=c(colorbar[1:10],"#FFFFFF",colorbar[11:20])
  colorbar=c(rep(colorbar[1],4),colorbar,rep(colorbar[21],4))

p=plot_ly(lat = lat, lon = lon, color = lm.trends[,1,1],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>January', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.jan.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends[,2,1],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>Febuary', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.feb.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends[,3,1],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>March', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.mar.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends[,4,1],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>April', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.apr.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends[,5,1],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>May', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.may.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends[,6,1],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>Jun', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.jun.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends[,7,1],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>July', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.jul.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends[,8,1],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>August', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.aug.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends[,9,1],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>September', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.sep.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends[,10,1],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>October', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.oct.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends[,11,1],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>November', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.nov.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends[,12,1],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>December', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.dec.png")

  
  
  
  
  
  
  
####################################################
#####  Now plot with the non-sig values as NA  #####
####################################################

Sys.setenv("plotly_username" = "333tristan")
Sys.setenv("plotly_api_key" = "h3bhvlgnr9")

# marker styling
m <- list(
  colorbar = list(title = "Trend in Humidity During Heat Events"),
  size = 8, opacity = 0.8, symbol = 'square',cmin=-.25,cmax=.25#,showscale=F
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
colorbar=rev(colorRampPalette(brewer.pal(11,"RdBu"))(20))
  colorbar=c(colorbar[1:10],"#FFFFFF",colorbar[11:20])
  colorbar=c(rep(colorbar[1],4),colorbar,rep(colorbar[21],4))

p=plot_ly(lat = lat, lon = lon, color = lm.trends.adj[,1],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>January', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.jan.sig.only.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends.adj[,2],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>Febuary', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.feb.sig.only.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends.adj[,3],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>March', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.mar.sig.only.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends.adj[,4],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>April', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.apr.sig.only.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends.adj[,5],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>May', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.may.sig.only.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends.adj[,6],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>Jun', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.jun.sig.only.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends.adj[,7],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>July', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.jul.sig.only.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends.adj[,8],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>August', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.aug.sig.only.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends.adj[,9],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>September', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.sep.sig.only.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends.adj[,10],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>October', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.oct.sig.only.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends.adj[,11],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>November', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.nov.sig.only.png")
  

p=plot_ly(lat = lat, lon = lon, color = lm.trends.adj[,12],
        type = 'scattergeo', locationmode = 'USA-states', mode = 'markers',
        marker = m, colors=colorbar) %>%
  layout(title = 'Humidity Trend During Heat Events<br>December', geo = g)
  plotly_IMAGE(p, format = "png", out_file = "plot.shum.hw.trend.dec.sig.only.png")

  
  
