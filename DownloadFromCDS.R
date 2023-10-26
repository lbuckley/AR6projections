library(ncdf4)
library(raster)
library(terra)

#Climate factory
#https://esgf-node.ipsl.upmc.fr/projects/cmip6-adjust/
#https://esgf-node.ipsl.upmc.fr/search/cmip6-adjust/
  
#raster strategy
library(raster)
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/")
nc.daily <- nc_open("tasmaxAdjust_day_GFDL-ESM4_ssp585_r1i1p1f1_gr010_TCDF-CDFT23-ERA5Land-1981-2010_19510101-19651231.nc")

lon <- ncvar_get(nc.daily,"lon");# extracts longitude
lat <- ncvar_get(nc.daily,"lat");# extracts latitude
time <- ncvar_get(nc.daily,"time");# extracts time

lon.ind= which(lon>251 & lon<258) 
lat.ind= which(lat>37 & lat<41)

nc.sub <- ncvar_get(nc.daily, "tasmaxAdjust",
                     start = c(lon.ind[1], lat.ind[1], 1),
                     count = c(length(lon.ind),length(lat.ind),length(time)))

nc.rast<- raster(nc.sub[,,1])
plot(nc.rast)

