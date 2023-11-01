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
#WUS-D3
#https://registry.opendata.aws/wrf-cmip6/
# browse https://wrf-cmip6-noversioning.s3.amazonaws.com/index.html
# products: https://dept.atmos.ucla.edu/alexhall/downscaling-cmip6
# CNRM-ESM2-1 SSP3-7.0 r1i1p1f2 from 2015-2100*

#https://wrf-cmip6-noversioning.s3.amazonaws.com/index.html#downscaled_products/gcm/mpi-esm1-2-hr_r3i1p1f1_historical_bc/6hourly/1980/d03/
#wrfout_d01_1980-09-01_00:00:00

#Tier3
#wspd10max.daily.mpi-esm1-2-lr.r71ip1f1.ssp370.d02.2092.nc.
# /<data_name>/postprocessed/<domain>

#read data from AWS
#https://blog.djnavarro.net/posts/2022-03-17_using-aws-s3-in-r/

library(aws.s3)
library(tibble)

bucket_exists(
  bucket = "s3://wrf-cmip6-noversioning/", 
  region = "us-west-2"
)

wus_files <- get_bucket_df(
  bucket = "s3://wrf-cmip6-noversioning/downscaled_products/gcm/mpi-esm1-2-hr_r3i1p1f1_historical_bc/postprocess/d02/", 
  region = "us-west-2",
  max = 20000
) %>% 
  as_tibble()

#tier 3
# <variable>.daily.<gcm>.<variant>.<exp_id>.<domain>.<year>.nc
# wspd10max.daily.mpi-esm1-2-lr.r71ip1f1.ssp370.d02.2092.nc

# #variables
# 2-m minimum temperature, t2min
# 2-m maximum temperature, t2max
# 2-m specific humidity, q2
# Soil temperature, soil_t
# Mean 10-m wind speed, wspd10mean
# Net SW flux at the surface (> 0 into sfc), sw_sfc
# Net LW flux at surface (> 0 into atm), lw_sfc
vars= c("t2min", "t2max", "q2", "soil_t", "wspd10mean", "sw_sfc", "lw_sfc")

# time periods
# mpi-esm1-2-hr_r3i1p1f1_historical_bc
# historical 1980-2013
# mpi-esm1-2-hr_r3i1p1f1_ssp370_bc
# 2014 - 2099   
years= c(2015:2020, 2065:2070)

for (year.k in 1:length(years)){
  for (var.k in 1:length(vars)){
    
    filename= paste("downscaled_products/gcm/mpi-esm1-2-hr_r3i1p1f1_ssp370_bc/postprocess/d02/", vars[var.k],".daily.mpi-esm1-2-hr.r3i1p1f1.ssp370.bias-correct.d02.", years[year.k],".nc", sep="") 
    outname= paste("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/WUS-D3/mpi-esm1-2-hr.r3i1p1f1/", vars[var.k],"_", years[year.k],".nc", sep="")
    
    if(var.k==1 & year.k==1)  {filenames= filename; outnames= outname}
    if(var.k>1 | year.k>1)  {filenames= c(filenames, filename); outnames= c(outnames, outname)}
  }} #end loop vars and years

# read coordinates
# https://wrf-cmip6-noversioning.s3.amazonaws.com/index.html#downscaled_products/wrf_coordinates/
coord <- nc_open("wrfinput_d02_coord.nc")
lon <- ncvar_get(coord,"lon2d");# extracts longitude
lat <- ncvar_get(coord,"lat2d");# extracts latitude
lon<- lon[,1]
lat<- lat[1,]

#download nc files
for (file.k in 8:length(filenames)){ #length(filenames)

save_object(
  object = print(filenames[file.k]),
  bucket = "s3://wrf-cmip6-noversioning/", 
  region = "us-west-2",
  file = print(outnames[file.k]))
}  

#Read nc in
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/WUS-D3/mpi-esm1-2-hr.r3i1p1f1/")

year.k=1

t2min= nc_open( paste("t2min_", years[year.k],".nc", sep=""))
t2max= nc_open( paste("t2max_", years[year.k],".nc", sep=""))
q2= nc_open( paste("q2_", years[year.k],".nc", sep=""))
soil_t= nc_open( paste("soil_t_", years[year.k],".nc", sep=""))
wspd10mean= nc_open( paste("wspd10mean_", years[year.k],".nc", sep=""))
sw_sfc= nc_open( paste("sw_sfc_", years[year.k],".nc", sep=""))
lw_sfc= nc_open( paste("lw_sfc_", years[year.k],".nc", sep=""))

time <- ncvar_get(soil_t,"day");# extracts time

lon.ind= which(lon> -109 & lon< -102) 
#lon.ind= which(lon>251 & lon<258) 
lat.ind= which(lat>37 & lat<41)

t2min.sub <- ncvar_get(t2min, "t2min", 
                    start = c(lon.ind[1], lat.ind[1], 1),
                    count = c(length(lon.ind),length(lat.ind),length(time))) 
t2max.sub <- ncvar_get(t2max, "t2max", 
                       start = c(lon.ind[1], lat.ind[1], 1),
                       count = c(length(lon.ind),length(lat.ind),length(time))) 
q2.sub <- ncvar_get(q2, "q2", 
                       start = c(lon.ind[1], lat.ind[1], 1),
                       count = c(length(lon.ind),length(lat.ind),length(time))) 
soil_t.sub <- ncvar_get(soil_t, "soil_t", 
                       start = c(lon.ind[1], lat.ind[1], 1, 1),
                       count = c(length(lon.ind),length(lat.ind),4, length(time))) #get 4 soil layers 
wspd10mean.sub <- ncvar_get(wspd10mean, "wspd10mean", 
                       start = c(lon.ind[1], lat.ind[1], 1),
                       count = c(length(lon.ind),length(lat.ind),length(time))) 
sw_sfc.sub <- ncvar_get(sw_sfc, "sw_sfc", 
                       start = c(lon.ind[1], lat.ind[1], 1),
                       count = c(length(lon.ind),length(lat.ind),length(time))) 
lw_sfc.sub <- ncvar_get(lw_sfc, "lw_sfc", 
                       start = c(lon.ind[1], lat.ind[1], 1),
                       count = c(length(lon.ind),length(lat.ind),length(time))) 

#time series, check out soil temperatures
plot(soil_t.sub[1,1,1,], type="l")
points(soil_t.sub[1,1,2,], type="l", col="red")
points(soil_t.sub[1,1,3,], type="l", col="blue")
points(soil_t.sub[1,1,4,], type="l", col="green")

#make rasters
t2min.rast<- raster(t2min.sub[,,1], xmn=min(lon[lon.ind]), xmx=max(lon[lon.ind]), ymn=min(lat[lat.ind]), ymx=max(lat[lat.ind]))
t2max.rast<- raster(t2max.sub[,,1], xmn=min(lon[lon.ind]), xmx=max(lon[lon.ind]), ymn=min(lat[lat.ind]), ymx=max(lat[lat.ind]))
q2.rast<- raster(q2.sub[,,1], xmn=min(lon[lon.ind]), xmx=max(lon[lon.ind]), ymn=min(lat[lat.ind]), ymx=max(lat[lat.ind]))
soil_t.rast<- raster(soil_t.sub[,,1,1], xmn=min(lon[lon.ind]), xmx=max(lon[lon.ind]), ymn=min(lat[lat.ind]), ymx=max(lat[lat.ind]))
wspd10mean.rast<- raster(wspd10mean.sub[,,1], xmn=min(lon[lon.ind]), xmx=max(lon[lon.ind]), ymn=min(lat[lat.ind]), ymx=max(lat[lat.ind]))
sw_sfc.rast<- raster(sw_sfc.sub[,,1], xmn=min(lon[lon.ind]), xmx=max(lon[lon.ind]), ymn=min(lat[lat.ind]), ymx=max(lat[lat.ind]))
lw_sfc.rast<- raster(lw_sfc.sub[,,1], xmn=min(lon[lon.ind]), xmx=max(lon[lon.ind]), ymn=min(lat[lat.ind]), ymx=max(lat[lat.ind]))

plot(sw_sfc.rast)
plot(lw_sfc.rast)

#----------------------------------------
#microclimate and biophysical

#apply nichemaper

install.packages("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/NicheMapR_3.2.1.tgz", repos = NULL, 
                 type = .Platform$pkgType)
library(NicheMapR)

# FUN_ecto
# https://github.com/mrke/NicheMapR/blob/master/R/FUN_ecto.R

#TrenchR
#scale to organism height
#biophysical model: air, surface temp, windspeed, windspeed, solar radiation








