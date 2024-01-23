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
# Figure out subset for extracting data

year.k=1
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/WUS-D3/mpi-esm1-2-hr.r3i1p1f1/")
t2min= nc_open( paste("t2min_", years[year.k],".nc", sep=""))

#load all
t2min.sub <- ncvar_get(t2min, "t2min", 
                        start = c(1, 1, 1),
                        count = c(length(lon[,1]),length(lat[1,]),1)) 

#plot
pts <- cbind(lat=as.vector(lat), lon=as.vector(lon), temp=as.vector(soil_t.sub[,,1]))

##hangs
# r <- rasterFromXYZ(pts)
# plot(r)

#plot out points with state boundaries
plot.w= ggplot(data=as.data.frame(pts), aes(x=lon, y=lat))+
  geom_point(aes(color=temp))+
  borders("state")

#subset to CO, 3138 cells
pts.ind= which(pts[,"lat"]>37 & pts[,"lat"]<41 & pts[,"lon"]< -102 & pts[,"lon"]> -109)

pts.sub= as.data.frame(pts[pts.ind,])

plot.co= ggplot(data=as.data.frame(pts.sub), aes(x=lon, y=lat))+
  geom_point(aes(color=temp), size=3)+
  borders("state") +xlim(-110,-101) +ylim(36,42)

#plot.co= ggplot(data=as.data.frame(pts.sub), aes(x=lon, y=lat))+
#  geom_tile(aes(fill = temp))+
#  #borders("state")+
#  xlim(-110,-101) +ylim(36,42)

#alternative strategy of subetting read
lon.ind= which(lon[,150]> -113 & lon[,150]< -103) 
lat.ind= which(lat[150,]>34 & lat[150,]<40)

soil_t.sub <- ncvar_get(soil_t, "soil_t", 
                        start = c(lon.ind[1], lat.ind[1], 1, 1),
                        count = c(length(lon.ind),length(lat.ind),1, 1)) 

#account for lat lon matrices
pts.sub <- cbind(lat=as.vector(lat[lon.ind,lat.ind]), lon=as.vector(lon[lon.ind,lat.ind]), temp=as.vector(soil_t.sub[,]))

ggplot(data=as.data.frame(pts.sub), aes(x=lon, y=lat))+
  geom_point(aes(color=temp), size=3)+
  borders("state")

#---------------------------------
#load and process data

for (year.k in 1:1){ #length(years)
  
  t2min= nc_open( paste("t2min_", years[year.k],".nc", sep=""))
  t2max= nc_open( paste("t2max_", years[year.k],".nc", sep=""))
  q2= nc_open( paste("q2_", years[year.k],".nc", sep=""))
  soil_t= nc_open( paste("soil_t_", years[year.k],".nc", sep=""))
  wspd10mean= nc_open( paste("wspd10mean_", years[year.k],".nc", sep=""))
  sw_sfc= nc_open( paste("sw_sfc_", years[year.k],".nc", sep=""))
  lw_sfc= nc_open( paste("lw_sfc_", years[year.k],".nc", sep=""))
  
  time <- ncvar_get(soil_t,"day");# extracts time
  
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
                          start = c(lon.ind[1], lat.ind[1], 1,1),
                          count = c(length(lon.ind),length(lat.ind),1,length(time)))#get 1st soil layers 
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
  
  #convert to points
  #check order
  lat.flat= flatten(lat[lon.ind,lat.ind], across = c("rows"))
  lon.flat= flatten(lon[lon.ind,lat.ind], across = c("rows"))
  pts.t= cbind(lat=lat.flat, lon=lon.flat)
  #with time
  #pts.t= cbind(lat=rep(lat.flat, length(time[1:2])), lon=rep(lon.flat, length(time[1:2])),time=rep(time[1:2], length(lat.flat)))
               
  
  #------------         
  t2min= flatten(t2min.sub[,,1], across = c("rows"))
  
  clim= cbind(pts.t, t2min= t2min)
  
  plot.co= ggplot(data=as.data.frame(clim), aes(x=lon, y=lat))+
    geom_point(aes(color=t2min), size=3)
  
  #-----
  clim= cbind(pts.t, t2min= flatten(t2min.sub[,,1], across = c("rows")),
                    t2max= flatten(t2max.sub[,,1], across = c("rows")),
                    q2= flatten(q2.sub[,,1], across = c("rows")),
                    soil_t= flatten(soil_t.sub[,,1], across = c("rows")),
                    wspd10mean= flatten(wspd10mean.sub[,,1], across = c("rows")),
                    sw_sfc= flatten(sw_sfc.sub[,,1], across = c("rows")),
                    lw_sfc= flatten(lw_sfc.sub[,,1], across = c("rows")) 
                    )
  
  #subset to Colorado points
  clim.ind= which(clim[,"lat"]>37 & clim[,"lat"]<41 & clim[,"lon"]< -102 & clim[,"lon"]> -109)
  
  clim.sub= as.data.frame(clim[clim.ind,])
  
  #order messed up
  plot.co= ggplot(data=as.data.frame(clim.sub[which(clim.sub$time=="20150901"),]), aes(x=lon, y=lat))+
    geom_point(aes(color=t2min), size=3)
   # + borders("state") +xlim(-110,-101) +ylim(36,42)
  
} #end loop years


#extract points in raster
lat.flat= flatten(lat, across = c("rows"))
lon.flat= flatten(lon, across = c("rows"))
pts= cbind(lon=lon.flat, lat=lat.flat)

#alternative strategy of subetting read
pts.ind= which(pts[,1]> -113 & pts[,1]< -103) 
pts.s= pts[pts.ind,]
pts.s= pts.s[which(pts.s[,2]> 34 & pts.s[,2]< 40),]

#-----
#WORKS
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

inds= lat.ind[matched,]

r2 <- extract(r1, inds)
lats<- lat[inds]
lons<- lon[inds]

dat= cbind(lons,lats,r2)
plot.co= ggplot(data=as.data.frame(dat), aes(x=lons, y=lats))+
  geom_point(aes(color=r2), size=4)


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








