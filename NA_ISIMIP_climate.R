library(ncdf4)
library(raster)
library(terra)
library(fields)
#install.packages("PCICt")
library(PCICt)
#install.packages("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/NicheMapR_3.2.1.tgz", repos = NULL, type = .Platform$pkgType)
library(NicheMapR)
library(TrenchR)
library(ggplot2)
library(tidyr)
library(plyr)
#---------------
#Prepares ISIMIP daily data for NA and runs microclimate and biophysical model

#ISIMIP
#https://www.isimip.org/

#for Butterflies?
#0.5 degree resolution
#https://www.isimip.org/gettingstarted/input-data-bias-adjustment/details/103/
#https://www.isimip.org/gettingstarted/input-data-bias-adjustment/details/102/
#Download from: https://data.isimip.org/search/tree/ISIMIP3b/SecondaryInputData/climate/atmosphere/mpi-esm1-2-hr/

# #variables
# sfcWind: near-surface wind speed, m s^-1
# tasmin, tasmax, tas: daily min, max, ave temp, K
# rlds: long wave downwelling radiation, W m^−2
# rsds: short wave downwelling radiation, W m^−2 
# pr: total precipiation, kg m^-2 s^-1
#hurs, huss: near surface-relative, specific humidity, %, kg kg^−1

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/ISIMIP/")

tasmax<- nc_open("mpi-esm1-2-hr_r1i1p1f1_w5e5_ssp245_tasmax_global_daily_2021_2030.nc")

time <- ncvar_get(tasmax,"time");# extracts time, days
lon <- ncvar_get(tasmax,"lon");# extracts longitude
lat <- ncvar_get(tasmax,"lat");# extracts latitude

tasmin<- nc_open("mpi-esm1-2-hr_r1i1p1f1_w5e5_ssp245_tasmin_global_daily_2021_2030.nc")
sfcwind<- nc_open("mpi-esm1-2-hr_r1i1p1f1_w5e5_ssp245_sfcwind_global_daily_2021_2030.nc")
rsds<- nc_open("mpi-esm1-2-hr_r1i1p1f1_w5e5_ssp245_rsds_global_daily_2021_2030.nc")
hurs<- nc_open("mpi-esm1-2-hr_r1i1p1f1_w5e5_ssp245_hurs_global_daily_2021_2030.nc")

lon.ind= which(lon> -130 & lon< -60) 
lat.ind= which(lat>20 & lat<50)

#make list of sites
grids = expand.grid(lon.ind= lon.ind, lat.ind = lat.ind)
grids$lats<- lat[lat.ind]
grids$lons<- lon[lon.ind]

# #time point at all sites
# tasmax.us <- ncvar_get(tasmax, "tasmax", 
#                        start = c(lon.ind[1], lat.ind[1], 1),
#                        count = c(length(lon.ind),length(lat.ind),1)) 
# 
# nc.rast<- raster(tasmax.us, ymn=min(lon[lon.ind]), ymx=max(lon[lon.ind]), xmn=min(lat[lat.ind]), xmx=max(lat[lat.ind]))
# nc.rast<-t(nc.rast) #(flip(nc.rast, 1))
# plot(nc.rast)

#map time
#proleptic_gregorian: A Gregorian calendar extended to dates before 1582-10-15. That is, a year is a leap year if either (i) it is divisible by 4 but not by 100 or (ii) it is divisible by 400
cal <- "proleptic_gregorian"
origin <- "2021-01-01"
seconds.per.day <- 86400
ts.dat.days <- time-0.5
origin.pcict <- as.PCICt(origin, cal)
time.date <- origin.pcict + (ts.dat.days * seconds.per.day)
#doy
doys<- day_of_year(time.date, format = "%Y-%m-%d")

#time series for each site
grid.ind<- 200

tasmax.ts <- ncvar_get(tasmax, "tasmax", 
                       start = c(grids$lon.ind[grid.ind], grids$lon.ind[grid.ind], 1),
                       count = c(1,1,length(time))) 
tasmin.ts <- ncvar_get(tasmin, "tasmin", 
                       start = c(grids$lon.ind[grid.ind], grids$lon.ind[grid.ind], 1),
                       count = c(1,1,length(time))) 
sfcwind.ts <- ncvar_get(sfcwind, "sfcwind", 
                       start = c(grids$lon.ind[grid.ind], grids$lon.ind[grid.ind], 1),
                       count = c(1,1,length(time))) 
rsds.ts <- ncvar_get(rsds, "rsds", 
                        start = c(grids$lon.ind[grid.ind], grids$lon.ind[grid.ind], 1),
                        count = c(1,1,length(time))) 
hurs.ts <- ncvar_get(hurs, "hurs", 
                        start = c(grids$lon.ind[grid.ind], grids$lon.ind[grid.ind], 1),
                        count = c(1,1,length(time))) 
#convert to C
tasmax.ts<- tasmax.ts-273.15
tasmin.ts<- tasmin.ts-273.15

clim.grid<- as.data.frame(cbind(time, doys, tasmax.ts, tasmin.ts, sfcwind.ts, rsds.ts, hurs.ts))

# #plot time series
# plot(1:length(time), tasmax.ts, type="l")
# points(1:length(time), tasmin.ts, type="l", col="blue")
# plot(1:length(time), sfcwind.ts, type="l")
# plot(1:length(time), rsds.ts, type="l")
# plot(1:length(time), hurs.ts, type="l")

#----
#LOAD BIOLOGICAL DATA
#Larval would mean the individual was reared at the treatment temperature from hatching
#pupal means the individuals were reared at 26 degrees until pupation, then switched into their treatments. 

#"PupalTempData" contains data for all three populations where individuals were reared at treatment temperatures only during the pupal stage. 
#"RearingTempPupalData" contains pupal data for the individuals that were reared at the treatment temperatures since hatching. 
#"LarvalCombData" contains all of the larval data for individuals that were kept at treatment temperatures since hatching. 
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/INTBIO/data/")
ldat<- read.csv("LarvalCombData.csv")
pdat<- read.csv("PupalTempData.csv")
rtdat<- read.csv("RearingTempPupalData.csv")

#proportion surviving to egg
#proportion surviving hatching to pupation
ldat$Surv.p<- 0
ldat$Surv.p[ldat$Surv.pupa=="Y"]<- 1
ldat.s<- aggregate(ldat$Surv.p, by = list(Treatment= ldat$Treatment, Pop=ldat$Population), FUN = "mean")

#survival curves by time and temperature
ggplot(data=ldat.s, aes(x=Treatment, y=x, color=Pop))+
  geom_point()+
  geom_smooth()

#pupal development rate
pd.plot= ggplot(data=pdat, aes(x=Treatment, y=1/DaysPtoE, color=Population))+
  geom_point()+
  geom_smooth()

ggplot(data=pdat, aes(x=Treatment, y=1/DaysHtoP, color=Population))+
  geom_point()+
  geom_smooth()

#larval mass
#to long format
ldat2<- ldat[,c(1:3,20:21,5:6,8:9,13:14)]
ldat.l<- ldat2 %>%
  gather("metric", "value", 6:ncol(ldat2))
ldat.l$type<- "mass"
ldat.l$type[ldat.l$metric %in% c("Time.4th","Time.5th","Age.pupa")]<- "time"
ldat.l$instar<-4
ldat.l$instar[ldat.l$metric %in% c("Time.5th", "Mass.5th")]<-5
ldat.l$instar[ldat.l$metric %in% c("Age.pupa", "Mass.pupa")]<-6

ggplot(data=ldat.l, aes(x=instar, y=value, color=Treatment, group=Treatment))+
  geom_point()+
  geom_smooth()+
  facet_grid(type~Population, scales="free_y")

#development time
ggplot(data=ldat.l[which(ldat.l$type=="time"),], aes(x=Treatment, y=1/value, color=Population))+
  geom_point()+
  geom_smooth()+
  facet_wrap(.~instar)

#pupal mass
ggplot(data=ldat.l[which(ldat.l$type=="mass"),], aes(x=Treatment, y=value, color=Population))+
  geom_point()+
  geom_smooth()+
  facet_wrap(.~instar)

#adult mass
ggplot(data=pdat, aes(x=Treatment, y=Mass.adult, color=Population))+
  geom_point()+
  geom_smooth()

#TO MODEL
#Larval development rate: phenology
#Survival using survival curves: survival
#Adult mass: relate to fedundity? Not clear growth TPC

#----
#microclimate and biophysical

#MICROCLIMATE
#NicheMapR microclimate model
#https://github.com/mrke/NicheMapR/blob/master/R/microrun.R

#time of sunrise and sunset
daylengths<- sapply(clim.grid$doys, FUN= daylength, lat= grids$lats[grid.ind])
solarnoons<- sapply(clim.grid$doys, FUN= solar_noon, lon= grids$lons[grid.ind])

clim.grid$t_r= solarnoons - daylengths/2
clim.grid$t_s= solarnoons + daylengths/2
clim.grid$t0<- solarnoons

# scale wind speed
clim.grid$Uz <- wind_speed_profile_neutral(u_r = clim.grid$sfcwind.ts, zr = 2, z0 = 0.02, z = z)

#diurnal temperature variation
##use sine exp, need to fix Gamma
#Thr<- apply(clim.grid, MARGIN=1, FUN= function(x) diurnal_temp_variation_sineexp(T_max=x[3], T_min=x[4], t=1:24, t_r=x[8], t_s=x[9])) 
# CO values, https://trenchproject.github.io/TrenchR/reference/diurnal_temp_variation_sineexp.html
#use sine
Thr<- apply(clim.grid, MARGIN=1, FUN= function(x) diurnal_temp_variation_sine(T_max=x[3], T_min=x[4], t=1:24))

#diurnal solar variation OPTIONS
#*** Currently reestimating soalr radiation below
#Shr<- apply(clim.grid, MARGIN=1, FUN= function(x) diurnal_radiation_variation(doy=x[2], S=x[6]*24, hour=1:24, lon=grids$lons[grid.ind], lat=grids$lats[grid.ind]))
#need to fix format
#Shr<- matrix(unlist(Shr), ncol = 24, byrow = TRUE)
#restimate solar radiation 
#S_sdir= direct_solar_radiation(lat=grids$lats[grid.ind], doy=clim.grid$doys, elev=1800, t=hr, t0=clim.grid$t0, method = "Campbell 1977")
##estimate diffuse fraction of solar radiation
#Sdf<- partition_solar_radiation(method= "Erbs", kt= 0.5, lat= grids$lats[grid.ind])

#----
#Loop through days for soil estimates and scaling
#make matrix to store data
Tsoil= array(NA, dim=c(nrow(clim.grid),2, 24)) #dim 2: sun then shade
Srad= array(NA, dim=c(nrow(clim.grid),3, 24))
Tz= matrix(NA, nrow=nrow(clim.grid), ncol=24)
Te= matrix(NA, nrow=nrow(clim.grid), ncol=24)
  
for(day.id in 1:nrow(clim.grid)){

#estimate solar radiation
Srad[day.id,,]<- sapply(degrees_to_radians(psi[,day.id]), FUN=solar_radiation, doy =clim.grid$doys[day.id], psi=, tau = tau, elev = 1800, rho = 0.6)

#estimate soil temperatures sun then shade
Tsoil[day.id,1,]<- soil_temperature(z_r.intervals = 12,
                         z_r           = 2,
                         z             = 2,
                         T_a           = Thr[,day.id],
                         u             = rep(clim.grid$"Uz"[1],24),
                         Tsoil0        = Thr[1,day.id],
                         z0            = 0.02,
                         SSA           = 0.7,
                         TimeIn        = 1:24,
                         S             = Srad[day.id,1,],
                         water_content = 0.2,
                         air_pressure  = 85,
                         rho_so        = 1620,
                         shade         = FALSE)

Tsoil[day.id,2,]<- soil_temperature(z_r.intervals = 12,
                                     z_r           = 2,
                                     z             = 2,
                                     T_a           = Thr[,day.id],
                                     u             = rep(clim.grid$"Uz"[1],24),
                                     Tsoil0        = Thr[1,day.id],
                                     z0            = 0.02,
                                     SSA           = 0.7,
                                     TimeIn        = 1:24,
                                     S             = Srad[day.id,1,],
                                     water_content = 0.2,
                                     air_pressure  = 85,
                                     rho_so        = 1620,
                                     shade         = TRUE)

## scale temperature from two meters
z <- 0.001 #specify distance from ground 
Tz[day.id,] <- air_temp_profile(T_r = Thr[,day.id], u_r = clim.grid$Uz[day.id], zr = 2, z0 = 0.02, z = z, T_s = Tsoil[day.id,2,])
# clim.grid$Tmaxz <- air_temp_profile(T_r = clim.grid$tasmax.ts, u_r = clim.grid$sfcwind.ts, zr = 2, z0 = 0.02, z = z, T_s = XXXX)
# #check surface roughness

} #end loop grid

#estimate zenith angle across hours
psi<- apply(clim.grid, MARGIN=1, FUN= function(x) zenith_angle(doy = x[2], lat=grids$lats[grid.ind], lon=grids$lons[grid.ind], hour = 1:24))

#------------------------
#BIOPHYSICAL
# FUN_ecto
# https://github.com/mrke/NicheMapR/blob/master/R/FUN_ecto.R

#TrenchR
#scale to organism height
#Tb_butterfly biophysical model: air, surface temp, windspeed, windspeed, direct and diffuse solar radiation

#add data to clim.grid
#loop through hours

#Te matrix
for(hr in 1:24){

#combine data for hour
hr.dat= cbind(Tz=Tz[,hr], Tg=Tsoil[,1,hr], Tsh=Tsoil[,2,hr], Uz=clim.grid$Uz, Sdir=Srad[,1,hr],Sdif=Srad[,2,hr], psi=psi[hr,])  

Te[,hr] <- apply(hr.dat, MARGIN=1, FUN= function(x) Tb_butterfly(T_a = x[1], T_g = x[2], T_sh = x[3], u= x[4], S_sdir = x[5], S_sdif = x[6], z= x[7], 
                                                          D      = 0.036, delta  = 1.46, alpha  = 0.6, r_g    = 0.3))

} #end loop hour

#plot Tes
plot(1:length(time), Te[,12], type="l", ylim=range(10,55))
plot(1:24, Te[1,], type="l", ylim=range(10,55))




