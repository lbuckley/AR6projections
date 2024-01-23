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

#find indices to extract
#make raster brick
r <- brick(paste("t2min_", years[year.k],".nc", sep=""), varname = "t2min")
#single time
r1 <- r[[which(getZ(r)=="20150901")]]

#find indices of nc corresponding to lat lon
lon.inds= which(lon> -113 & lon< -103)

lon.ind= which(lon> -113 & lon< -103,arr.ind = T) 
lat.ind= which(lat>34 & lat<40,arr.ind = T)

lon.ind.match<- paste(lon.ind[,1], lon.ind[,2], sep=".")
lat.ind.match<- paste(lat.ind[,1], lat.ind[,2], sep=".")
matched= na.omit(match(lon.ind.match,lat.ind.match))

co.inds= lat.ind[matched,]
#corresponding lats and lons
lats<- lat[inds]
lons<- lon[inds]

#retrieve and extract data
for (year.k in 1:1){ #length(years)
  
  t2min <- brick(paste("t2min_", years[year.k],".nc", sep=""), varname = "t2min")
  t2max <- brick(paste("t2max_", years[year.k],".nc", sep=""), varname = "t2max")
  q2 <- brick(paste("q2_", years[year.k],".nc", sep=""), varname = "q2")
  soil_t <- brick(paste("soil_t_", years[year.k],".nc", sep=""), varname = "soil_t")
  #get 1st soil layers 
  #soil fields are defined on the 4 Noah-MP soil levels of 5, 25, 70, and 150 centimeters
  wspd10mean <- brick(paste("wspd10mean_", years[year.k],".nc", sep=""), varname = "wspd10mean")
  sw_sfc <- brick(paste("sw_sfc_", years[year.k],".nc", sep=""), varname = "sw_sfc")
  lw_sfc <- brick(paste("lw_sfc_", years[year.k],".nc", sep=""), varname = "lw_sfc")
  
  #extract times
  times= getZ(t2min)
  
  #yearly times series for a site
  t1 <- t2min[inds[1,1],inds[1,2]]
  
  #NEED TO EXTRACT ACROSS TIME
  t2min <- extract(t2min[[which(getZ(t2min)=="20150901")]], inds)
  #add other variables
  
  #UPDATE
  #combine
  dat= cbind(lons,lats,t2min)
  
  #test plot
  plot.co= ggplot(data=as.data.frame(dat), aes(x=lons, y=lats))+
    geom_point(aes(color=t2min), size=4)
  
} #end loop years

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








