library(ncdf4)
library(raster)
library(terra)
library(fields)

#Climate factory
#citation https://doi.org/10.1016/j.dib.2022.108669
#https://esgf-node.ipsl.upmc.fr/projects/cmip6-adjust/
#https://esgf-node.ipsl.upmc.fr/search/cmip6-adjust/

#CMIP5 has radiation
#https://help.theclimatedatafactory.com/en/article/weather-variables-ipcc-ar5-1qasjno/
#https://esgf-node.ipsl.upmc.fr/projects/c3s-cmip5-adjust/
#https://esgf-node.ipsl.upmc.fr/search/cmip5/
  
#raster strategy
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/")

#test on CMIP 5
#nc.daily <- nc_open("tasmaxAdjust_day_GFDL-ESM4_ssp585_r1i1p1f1_gr010_TCDF-CDFT23-ERA5Land-1981-2010_19510101-19651231.nc")

# CMIP 6
nc.daily <- nc_open("tasAdjust_day_IPSL-CM6A-LR_ssp245_r1i1p1f1_gr010_TCDF-CDFT23-ERA5Land-1981-2010_20060101-20151231.nc")

lon <- ncvar_get(nc.daily,"lon");# extracts longitude
lat <- ncvar_get(nc.daily,"lat");# extracts latitude
time <- ncvar_get(nc.daily,"time");# extracts time

lon.ind= which(lon>251 & lon<258) 
lat.ind= which(lat>37 & lat<41)

nc.sub <- ncvar_get(nc.daily, "tasAdjust", #"tasmaxAdjust"
                     start = c(lon.ind[1], lat.ind[1], 1),
                     count = c(length(lon.ind),length(lat.ind),500)) #abbreviate time

nc.rast<- raster(nc.sub[,,1], xmn=min(lon[lon.ind]), xmx=max(lon[lon.ind]), ymn=min(lat[lat.ind]), ymx=max(lat[lat.ind]))

plot(nc.rast)

#time series
plot(nc.sub[1,1,]  )

# #variables
# hursAdjust: relative humidity %
# prAdjust: precipitation_flux kg m-2 s-1
# sfcWindAdjust
# tasAdjust
# tasmaxAdjust
# tasminAdjust

#----------------
#MACA statistically downscaled
#https://toolkit.climate.gov/tool/maca-cmip5-statistically-downscaled-climate-projections
#https://www.climatologylab.org/maca.html
#https://www.climatologylab.org/maca.html

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/MACA/")

nc.daily <- nc_open("macav2livneh_tasmax_IPSL-CM5A-MR_r1i1p1_rcp85_2006_2025_CONUS_daily.nc")

lon <- ncvar_get(nc.daily,"lon");# extracts longitude
lat <- ncvar_get(nc.daily,"lat");# extracts latitude
time <- ncvar_get(nc.daily,"time");# extracts time

lon.ind= which(lon>251 & lon<258) 
lat.ind= which(lat>37 & lat<41)

nc.sub <- ncvar_get(nc.daily, "air_temperature",
                    start = c(lon.ind[1], lat.ind[1], 1),
                    count = c(length(lon.ind),length(lat.ind),500)) #abbreviate time

nc.rast<- raster(nc.sub[,,1], xmn=min(lon[lon.ind]), xmx=max(lon[lon.ind]), ymn=min(lat[lat.ind]), ymx=max(lat[lat.ind]))

plot(nc.rast)

#time series
plot(nc.sub[1,1,]  )

#----------------
#CA dynamic downscaling
#https://loca.ucsd.edu/loca-version-2-for-california-ca-may-2023/


