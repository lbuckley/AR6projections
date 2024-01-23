library(ncdf4)
library(raster)
library(terra)
library(fields)
library(ramify)

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

#---------------------------------
# Subset to CO

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/WUS-D3/mpi-esm1-2-hr.r3i1p1f1/")

#find indices to extract
#make raster brick
r <- brick(paste("t2min_", years[year.k],".nc", sep=""), varname = "t2min")
#single time
r1 <- r[[which(getZ(r)=="20150901")]]

#find indices of nc corresponding to lat lon
lon.ind= which(lon> -113 & lon< -103,arr.ind = T) 
lat.ind= which(lat>34 & lat<40,arr.ind = T)

lon.ind.match<- paste(lon.ind[,1], lon.ind[,2], sep=".")
lat.ind.match<- paste(lat.ind[,1], lat.ind[,2], sep=".")
matched= na.omit(match(lon.ind.match,lat.ind.match))

co.inds= lat.ind[matched,]
#corresponding lats and lons
lats<- lat[inds]
lons<- lon[inds]
co.pts= cbind(co.inds, lats, lons)

#make array to hold data
clim.dat<- array(data = NA, dim = c(366,nrow(co.inds),7), dimnames = NULL)

#retrieve and extract data
for (year.k in 1:1){ #length(years)
  
  t2min <- brick(paste("t2min_", years[year.k],".nc", sep=""), varname = "t2min")
  t2max <- brick(paste("t2max_", years[year.k],".nc", sep=""), varname = "t2max")
  q2 <- brick(paste("q2_", years[year.k],".nc", sep=""), varname = "q2")
  soil_t <- brick(paste("soil_t_", years[year.k],".nc", sep=""), level=1, varname = "soil_t")
  #get 1st soil layers 
  #soil fields are defined on the 4 Noah-MP soil levels of 5, 25, 70, and 150 centimeters
  wspd10mean <- brick(paste("wspd10mean_", years[year.k],".nc", sep=""), varname = "wspd10mean")
  sw_sfc <- brick(paste("sw_sfc_", years[year.k],".nc", sep=""), varname = "sw_sfc")
  lw_sfc <- brick(paste("lw_sfc_", years[year.k],".nc", sep=""), varname = "lw_sfc")
  
  #extract times
  times= getZ(t2min)
  
  #yearly times series for a site
  # t2min[co.inds[1,1],co.inds[1,2]]
  
  #extract one time across sites
  #t2min <- extract(t2min[[which(getZ(t2min)=="20150901")]], inds)
  
  ##combine
  #dat= cbind(lons,lats,t2min)
  
  ##test plot
  #plot.co= ggplot(data=as.data.frame(dat), aes(x=lons, y=lats))+
  #  geom_point(aes(color=t2min), size=4)
  
  #store data across sites
  #just 10 for now
 clim.dat[1:length(times),1:10,1]<- apply(co.inds[1:10,], MARGIN=1, FUN=function(x) t2min[x[1],x[2]])
  clim.dat[1:length(times),1:10,2] <- apply(co.inds[1:10,], MARGIN=1, FUN=function(x) t2max[x[1],x[2]])
  clim.dat[1:length(times),1:10,3] <- apply(co.inds[1:10,], MARGIN=1, FUN=function(x) q2[x[1],x[2]])
  clim.dat[1:length(times),1:10,4] <- apply(co.inds[1:10,], MARGIN=1, FUN=function(x) soil_t[x[1],x[2]])   
  clim.dat[1:length(times),1:10,5] <- apply(co.inds[1:10,], MARGIN=1, FUN=function(x) wspd10mean[x[1],x[2]])
  clim.dat[1:length(times),1:10,6] <- apply(co.inds[1:10,], MARGIN=1, FUN=function(x) sw_sfc[x[1],x[2]])
  clim.dat[1:length(times),1:10,7] <- apply(co.inds[1:10,], MARGIN=1, FUN=function(x) lw_sfc[x[1],x[2]])
  #Check issue with radiation
  
  #slow across sites
  #clim.dat[1:length(times),1:1000,1]<- apply(co.inds[1:1000,], MARGIN=1, FUN=function(x) t2min[x[1],x[2]])
  
} #end loop years

#----------------------------------------
#microclimate and biophysical

#install.packages("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/NicheMapR_3.2.1.tgz", repos = NULL, type = .Platform$pkgType)
library(NicheMapR)
library(TrenchR)

#go through sites
site.k<- 1
doys<- 1:365

#MICROCLIMATE
#NicheMapR microclimate model
#https://github.com/mrke/NicheMapR/blob/master/R/microrun.R

#time of sunrise and sunset
daylengths<- sapply(doys, FUN= daylength, lat= co.pts[site.k,3])
solarnoons<- sapply(doys, FUN= solar_noon, lon= co.pts[site.k,4])

t_r= solarnoons - daylengths/2
t_s= solarnoons + daylengths/2

#diurnal temperature variation
#combine input
cdat<- cbind(clim.dat[1:365,site.k,c(1,2,6)], t_r, t_s)
#integrate lw and sw

clim.dat[1:length(times),1:10,1]

Thr<- sapply(doys, FUN= function(x) diurnal_temp_variation_sineexp(T_max=x[2]+10, T_min=x[1], t=1:24, t_r=x[4], t_s=x[5], alpha = 1.52, beta = 2.00, gamma = -0.18))

diurnal_temp_variation_sineexp(
  T_max,
  T_min,
  t,
  t_r,
  t_s,
  alpha = 1.52,
  beta = 2.00,
  gamma = -0.18
)
# CO values, https://trenchproject.github.io/TrenchR/reference/diurnal_temp_variation_sineexp.html

#diurnal solar variation
diurnal_radiation_variation(doy, S, hour, lon, lat)

#estimate zenith angle
dat$psi <- zenith_angle(doy = dat$J, lat = dat$lat, lon = dat$lon, hour = dat$hour)

z <- 0.001 #specify distance from ground 

# scale temperature 
dat$Tgrass <- air_temp_profile(T_r = dat$Temp, u_r = dat$Wind, zr = 0.5, z0 = 0.02, z = z, T_s = dat$SoilTemp)

# scale wind speed
dat$Ugrass <- wind_speed_profile_neutral(u_r = dat$Wind, zr = 0.5, z0 = 0.02, z = z)

#BIOPHYSICAL
# FUN_ecto
# https://github.com/mrke/NicheMapR/blob/master/R/FUN_ecto.R

#TrenchR
#scale to organism height
#biophysical model: air, surface temp, windspeed, windspeed, solar radiation

dat$Te <- Tb_grasshopper(T_a = dat$Tgrass, T_g = dat$SoilTemp, u = dat$Ugrass, S = dat$Rad, K_t = 1, psi = dat$psi, l = 0.0211, Acondfact = 0.0, abs = 0.9, r_g = 0.5) 







