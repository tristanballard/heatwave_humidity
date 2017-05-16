suppressMessages(library(RNetCDF))

fileNames=read.csv("/scratch/users/tballard/filenames.csv", head=F)
for (i in 28:30){
  nc = open.nc(paste('/scratch/PI/omramom/CMIP5/historical/daily/1deg/',as.character(fileNames[i,1]),sep=''))
  lw = var.get.nc(nc,'huss')
  saveRDS(lw,paste('/scratch/PI/omramom/CMIP5/historical/daily/1deg/',as.character(fileNames[i,1]),'.rds',sep=''))
  rm(lw)
}


#check 17
#done