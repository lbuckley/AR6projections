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
#daily 6km or 4km

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/MACA/")

nc.daily <- nc_open("macav2livneh_tasmax_IPSL-CM5A-MR_r1i1p1_rcp85_2006_2025_CONUS_daily.nc")

lon <- ncvar_get(nc.daily,"lon");# extracts longitude
lat <- ncvar_get(nc.daily,"lat");# extracts latitude
time <- ncvar_get(nc.daily,"time");# extracts time

lon.ind= which(lon>251 & lon<258) 
lat.ind= which(lat>37 & lat<41)

#all sites, one time
nc.sub <- ncvar_get(nc.daily, "air_temperature",
                    start = c(1, 1, 1),
                    count = c(length(lon),length(lat),1)) #abbreviate time

nc.rast<- raster(nc.sub, ymn=min(lon), ymx=max(lon), xmn=min(lat), xmx=max(lat))
nc.rast<- t(flip(nc.rast, 1))

#one site, all times
nc.sub <- ncvar_get(nc.daily, "air_temperature",
                    start = c(200, 200, 1),
                    count = c(1,1,length(time))) #abbreviate time
plot(1:length(time), nc.sub, type="l")
#how is time coded?

#--
nc.sub <- ncvar_get(nc.daily, "air_temperature",
                    start = c(lon.ind[1], lat.ind[1], 1),
                    count = c(length(lon.ind),length(lat.ind),500)) #abbreviate time



nc.rast<- raster(nc.sub, xmn=min(lon[lon.ind]), xmx=max(lon[lon.ind]), ymn=min(lat[lat.ind]), ymx=max(lat[lat.ind]))

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

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/WUS-D3/")
coord <- nc_open("wrfinput_d02_coord.nc")
lon <- ncvar_get(coord,"lon2d");# extracts longitude
lat <- ncvar_get(coord,"lat2d");# extracts latitude

#download nc files
# for (file.k in 8:length(filenames)){ #length(filenames)
# 
# save_object(
#   object = print(filenames[file.k]),
#   bucket = "s3://wrf-cmip6-noversioning/", 
#   region = "us-west-2",
#   file = print(outnames[file.k]))
# }  

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
#soil fields are defined on the 4 Noah-MP soil levels of 5, 25, 70, and 150 centimeters
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

plot(t2min.rast)
plot(sw_sfc.rast)
plot(lw_sfc.rast)

#read in bigger area to check
#map <- raster::getData("GADM",country='USA',level=1)

lon.ind= which(lon[1,]> -125 & lon[1,]< -80) 
#lon.ind= which(lon>251 & lon<258) 
lat.ind= which(lat[,1]>22 & lat[,1]<50)

soil_t.sub <- ncvar_get(soil_t, "soil_t", 
                        start = c(lon.ind[1], lat.ind[1], 1, 1),
                        count = c(length(lon.ind),length(lat.ind),1, length(time))) 
#load all
soil_t.sub <- ncvar_get(soil_t, "soil_t", 
                        start = c(1, 1, 1, 1),
                        count = c(length(lon[,1]),length(lat[1,]),1, 2)) 

#account for lat lon matrices
pts <- cbind(lat=as.vector(lat[lon.ind,lat.ind]), lon=as.vector(lon[lon.ind,lat.ind]), temp=as.vector(soil_t.sub[,,1]))
#load all
pts <- cbind(lat=as.vector(lat), lon=as.vector(lon), temp=as.vector(soil_t.sub[,,1]))

##hangs
# r <- rasterFromXYZ(pts)
# plot(r)

#plot out points with state boundaries
plot1= ggplot(data=as.data.frame(pts), aes(x=lon, y=lat))+
  geom_point(aes(color=temp))+
  borders("state")

# #raster plot
# soil_t.rast<- raster(t(soil_t.sub[,,1]), ymn=min(lat[lon.ind,lat.ind]), ymx=max(lat[lon.ind,lat.ind]), xmn=min(lon[lon.ind,lat.ind]), xmx=max(lon[lon.ind,lat.ind]), crs="+init=epsg:4326" )
# #soil_t.rast <- t(soil_t.rast)
# soil_t.rast <- flip(soil_t.rast, direction='y')
# plot(soil_t.rast)
# plot(map, add=TRUE)

#---------------
#OTHER APPROACHES

#GEB 2019 approach

#We used gridded projections of daily minimum and maximum air temperature (°C, 2 m height) for 1950–2099 from the CMIP5 multi-model ensemble. 
#The projections were statistically downscaled to 1/8° latitude–longitude (c. 12 km by 12 km) resolution using daily bias correction and constructed analogues (bias correction constructed analogs (BCCA) method) 5 archive at http://gdo-dcp.ucllnl.org/downscaled_cmip_projections/). 

#We estimated butterfly body temperatures (see Section 2.4) additionally using hourly microclimate data from the microclim dataset (Kearney, Isaac, & Porter, 2014); we used monthly climate normals (1960–1990) at 15 km resolution for substrate temperatures in 0 and 100% shade (°C), wind speed at the surface (1 cm, m/s) and solar zenith angle (°). We thus assume these conditions do not shift over time.

#---------------
#ISIMIP
#https://www.isimip.org/

#for Butterflies?
#0.5 degree resolution
#https://www.isimip.org/gettingstarted/input-data-bias-adjustment/details/103/
#https://www.isimip.org/gettingstarted/input-data-bias-adjustment/details/102/
#Download from: https://data.isimip.org/search/tree/ISIMIP3b/SecondaryInputData/climate/atmosphere/mpi-esm1-2-hr/

# #variables
# sfcWind: near-surface wind speed
# tasmin, tasmax, tas: daily min, max, ave temp
# rlds: long wave downwelling radiation
# rsds: short wave downwelling radiation 
# pr: tota
#hurs, huss: near surface-relative, specific humidity, 

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/ISIMIP/")

tasmax= nc_open("mpi-esm1-2-hr_r1i1p1f1_w5e5_ssp245_tasmax_global_daily_2021_2030.nc")

time <- ncvar_get(tasmax,"time");# extracts time, days
lon <- ncvar_get(tasmax,"lon");# extracts longitude
lat <- ncvar_get(tasmax,"lat");# extracts latitude

lon.ind= which(lon> -130 & lon< -60) 
lat.ind= which(lat>20 & lat<50)

#1 time all sites
tasmax.us <- ncvar_get(tasmax, "tasmax", 
                       start = c(lon.ind[1], lat.ind[1], 1),
                       count = c(length(lon.ind),length(lat.ind),1)) 

nc.rast<- raster(tasmax.us, ymn=min(lon[lon.ind]), ymx=max(lon[lon.ind]), xmn=min(lat[lat.ind]), xmx=max(lat[lat.ind]))
nc.rast<-t(nc.rast) #(flip(nc.rast, 1))
plot(nc.rast)

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

#FiX gamma constraint on DTR
diurnal_temp_variation_sineexp<- function (T_max, T_min, t, t_r, t_s, alpha = 2.59, beta = 1.55, 
          gamma = 2.2) 
{
  stopifnot(T_max >= T_min, t_s >= 0, t_s <= 24, t_r >= 0, 
            t_r <= 24, t >= 0, t <= 24)
  l <- t_s - t_r
  t_x <- 0.5 * (t_r + t_s) + alpha
  t_n <- t_r + beta
  if (!(t > (t_r + beta) & t < t_s)) {
    T_sn <- T_min + (T_max - T_min) * sin((pi * (t_s - t_r - 
                                                   beta))/(l + 2 * (alpha - beta)))
    if (t <= (t_r + beta)) {
      t_as <- t + 24 - t_s
    }
    if (t >= t_s) {
      t_as <- t - t_s
    }
    Temp <- T_min + (T_sn - T_min) * exp(-(gamma * t_as)/(24 - 
                                                            l + beta))
  }
  else {
    Temp <- T_min + (T_max - T_min) * sin((pi * (t - t_r - 
                                                   beta))/(l + 2 * (alpha - beta)))
  }
  Temp
}





