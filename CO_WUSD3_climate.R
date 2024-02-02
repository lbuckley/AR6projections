library(ncdf4)
library(raster)
library(terra)
library(fields)
library(ramify)

#DOWNLOAD DATA
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
#Process data
# Subset to CO

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/WUS-D3/mpi-esm1-2-hr.r3i1p1f1/")

#find indices to extract
#make raster brick
r <- brick(paste("t2min_", years[1],".nc", sep=""), varname = "t2min")
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
lats<- lat[co.inds]
lons<- lon[co.inds]
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
  
  # #extract times
  times= getZ(t2min)
  
  # #yearly times series for a site
  # t2min1<- t2min[co.inds[1,1],co.inds[1,2]]
  # t2max1<- t2max[co.inds[1,1],co.inds[1,2]]
  # q21<- q2[co.inds[1,1],co.inds[1,2]]
  # soil_t1<- soil_t[co.inds[1,1],co.inds[1,2]]
  # wspd10mean1<- wspd10mean[co.inds[1,1],co.inds[1,2]]
  # sw_sfc1<- extract(sw_sfc, co.inds[1]) #sw_sfc[co.inds[1,1],co.inds[1,2]]
  # lw_sfc1<- extract(lw_sfc, co.inds[1]) #lw_sfc[co.inds[1,1],co.inds[1,2]]
  # 
  # #temps
  # plot(1:366, t2min1, type="l", col="blue")
  # points(1:366, t2max1, type="l", col="red")
  # points(1:366, soil_t1, type="l", col="green")
  # 
  # #humidity and windspeed
  # plot(1:366, q21, type="l", col="blue", ylim=c(2,18))
  # points(1:366, wspd10mean1, type="l", col="red")
  # 
  # #solar radiation
  # plot(1:366, sw_sfc1, type="l", col="blue")
  # points(1:366, lw_sfc1, type="l", col="red")
  
  #---------------------
  
  #extract one time across sites
  #t2min <- extract(t2min[[which(getZ(t2min)=="20150901")]], inds)
  #sw_sfc1 <- extract(sw_sfc[[which(getZ(t2min)=="20150901")]], inds)
  
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
  #clim.dat[1:length(times),1:10,6] <- apply(co.inds[1:10,], MARGIN=1, FUN=function(x) sw_sfc[x[1],x[2]])
  #clim.dat[1:length(times),1:10,7] <- apply(co.inds[1:10,], MARGIN=1, FUN=function(x) lw_sfc[x[1],x[2]])
  #Check issue with radiation
  rad= apply(co.inds[1:10,], MARGIN=1, FUN=function(x) extract(sw_sfc, x) )
  clim.dat[1:length(times),1:10,6] <- apply(co.inds[1:10,], MARGIN=1, FUN=function(x) extract(sw_sfc, x)[1,] )
  clim.dat[1:length(times),1:10,7] <- apply(co.inds[1:10,], MARGIN=1, FUN=function(x) extract(lw_sfc, x)[2,] )
  
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

# scale temperature from two meters
z <- 0.001 #specify distance from ground 
Tminz <- air_temp_profile(T_r = clim.dat[1:365,site.k,c(1)], u_r = clim.dat[1:365,site.k,c(5)], zr = 2, z0 = 0.02, z = z, T_s = clim.dat[1:365,site.k,c(4)])
Tmaxz <- air_temp_profile(T_r = clim.dat[1:365,site.k,c(2)], u_r = clim.dat[1:365,site.k,c(5)], zr = 2, z0 = 0.02, z = z, T_s = clim.dat[1:365,site.k,c(4)])
#check surface roughness

# scale wind speed
Uz <- wind_speed_profile_neutral(u_r = clim.dat[1:365,site.k,c(5)], zr = 2, z0 = 0.02, z = z)

#plot
plot(1:365, clim.dat[1:365,site.k,c(1)], col="blue", type="l")
points(1:365, Tminz, col="blue", lty="dashed", type="l")
points(1:365, clim.dat[1:365,site.k,c(1)], col="red", type="l")
points(1:365, Tminz, col="red", lty="dashed", type="l")
#max and min close

#combine input
cdat<- cbind(clim.dat[1:365,site.k,c(1,2,6)], t_r, t_s, 1:365)
#Tmin, Tmax, sw

#diurnal temperature variation
##use sine exp, need to fix Gamma
#Thr<- apply(cdat, MARGIN=1, FUN= function(x) diurnal_temp_variation_sineexp(T_max=x[2], T_min=x[1], t=1:24, t_r=x[4], t_s=x[5])) #, alpha = 1.52, beta = 2.00, gamma = -0.18
# CO values, https://trenchproject.github.io/TrenchR/reference/diurnal_temp_variation_sineexp.html
#use sine
Thr<- apply(cdat, MARGIN=1, FUN= function(x) diurnal_temp_variation_sine(T_max=x[2], T_min=x[1], t=1:24))

#DTR for soil

#diurnal solar variation
Shr<- apply(cdat, MARGIN=1, FUN= function(x) unlist(diurnal_radiation_variation(doy=x[6], S=x[3], hour=1:24, lon=co.pts[site.k,"lons"], lat=co.pts[site.k,"lats"])))
#need to fix format

#estimate zenith angle
psi<- apply(cdat, MARGIN=1, FUN= function(x) zenith_angle(doy = x[6], lat=co.pts[site.k,"lats"], lon=co.pts[site.k,"lons"], hour = 1:24))

#BIOPHYSICAL
# FUN_ecto
# https://github.com/mrke/NicheMapR/blob/master/R/FUN_ecto.R

#TrenchR
#scale to organism height
#biophysical model: air, surface temp, windspeed, windspeed, solar radiation

#add data to cdat
#loop through hours

#for(hr in 1:24){
hr=10

  #combine data for hours
  gdat= cbind(Thr[hr,], clim.dat[1:365,site.k,c(4)], Uz, clim.dat[1:365,site.k,6], psi[hr,])  
  #fix radiation
  
Te <- apply(gdat, MARGIN=1, FUN= function(x) Tb_grasshopper(T_a = x[1], T_g = x[2], u = x[3], S = x[4], K_t = 1, psi = x[5], l = 0.0211, Acondfact = 0.0, abs = 0.9, r_g = 0.5))

plot(1:365, Te)
            






