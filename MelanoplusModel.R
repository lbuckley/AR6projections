library(ncdf4)
library(raster)
library(terra)
library(fields)
#install.packages("PCICt")
library(PCICt)
#install.packages("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/NicheMapR_3.2.1.tgz", repos = NULL, type = .Platform$pkgType)
library(NicheMapR)
library(TrenchR)
library(tidyr)
library(dplyr)

#Read species data
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/FitnessContrib_JEB/data/Data/Grasshopper/")
spec.dat=read.csv("SpecData.csv")
site.list= c("Eldorado","A1","B1","C1")
#specs= c("clav","pell","dodg","sang")
#spec.dat=spec.dat[match(spec.dat$SpecID, specs),]

# Jumping TPC parameters
tpc.dat= read.csv("JumpTPCparams.csv")

#sites
#Chautauqua Mesa (1,752 m, 40.00N, 105.28W), A1 (2,195 m, 40.01N, 105.37W), B1 (2,591 m, 40.02N, 105.43W), C1 (3,048 m, 40.03N, 105.55W)

##Calculate mass and length, m and g? 
#dat$mass_dodg= 1-1.64*10^-4*3048 #*dat$Elev
#dat$L_dodg= exp(3.33*0.247*log(dat$mass_dodg))

#==================================================
#read climate data
#solar, air temp, ground temp, wind

#WUS-D3
#https://registry.opendata.aws/wrf-cmip6/
# browse https://wrf-cmip6-noversioning.s3.amazonaws.com/index.html
# products: https://dept.atmos.ucla.edu/alexhall/downscaling-cmip6
# CNRM-ESM2-1 SSP3-7.0 r1i1p1f2 from 2015-2100*

# 2014 - 2099   
years= c(2015:2020, 2065:2070)
doys<- 1:365

# read coordinates
# https://wrf-cmip6-noversioning.s3.amazonaws.com/index.html#downscaled_products/wrf_coordinates/

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/WUS-D3/")
coord <- nc_open("wrfinput_d02_coord.nc")
lon <- ncvar_get(coord,"lon2d");# extracts longitude
lat <- ncvar_get(coord,"lat2d");# extracts latitude

#Process data
# Subset to CO

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/WUS-D3/mpi-esm1-2-hr.r3i1p1f1/")

#find indices to extract
#make raster brick
r <- brick(paste("t2min_", years[1],".nc", sep=""), varname = "t2min")
#single time
r1 <- r[[which(getZ(r)=="20150901")]]

#find indices of nc corresponding to lat lon
lon.ind= which(lon> -106 & lon< -105,arr.ind = T) 
lat.ind= which(lat>39 & lat<40.5,arr.ind = T)

lon.ind.match<- paste(lon.ind[,1], lon.ind[,2], sep=".")
lat.ind.match<- paste(lat.ind[,1], lat.ind[,2], sep=".")
matched= na.omit(match(lon.ind.match,lat.ind.match))

#--------------
co.inds= lat.ind[matched,]

#corresponding lats and lons
lats<- lat[co.inds]
lons<- lon[co.inds]
co.pts= cbind(co.inds, lats, lons)

#make array to hold data
clim.dat<- array(data = NA, dim = c(366,nrow(co.inds),7, length(years)), dimnames = NULL)
Te.dat<- array(data = NA, dim = c(nrow(co.inds),length(years), 366, 24), dimnames = NULL)
Th.dat<- array(data = NA, dim = c(nrow(co.inds),length(years), 366, 24), dimnames = NULL)

#retrieve and extract data
for (year.k in 1:length(years) ){ #length(years)
  
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
  
  # #extract one time across sites
  # t2min <- extract(t2min[[which(getZ(t2min)=="20150901")]], co.inds)
  # #sw_sfc1 <- extract(sw_sfc[[which(getZ(t2min)=="20150901")]], co.inds)
  # 
  # ##combine
  # dat= cbind(lons,lats,t2min)
  # 
  # ##test plot
  # plot.co= ggplot(data=as.data.frame(dat), aes(x=lons, y=lats))+
  #   geom_point(aes(color=t2min), size=4)
  
  #clim.dat<- array(data = NA, dim = c(366,nrow(co.inds),7, length(years)), dimnames = NULL)
  
  #store data across sites
  clim.dat[1:length(times),,1, year.k]<- apply(co.inds[,], MARGIN=1, FUN=function(x) t2min[x[1],x[2]])
  clim.dat[1:length(times),,2, year.k] <- apply(co.inds[,], MARGIN=1, FUN=function(x) t2max[x[1],x[2]])
  clim.dat[1:length(times),,3, year.k] <- apply(co.inds[,], MARGIN=1, FUN=function(x) q2[x[1],x[2]])
  clim.dat[1:length(times),,4, year.k] <- apply(co.inds[,], MARGIN=1, FUN=function(x) soil_t[x[1],x[2]])   
  clim.dat[1:length(times),,5, year.k] <- apply(co.inds[,], MARGIN=1, FUN=function(x) wspd10mean[x[1],x[2]])
  clim.dat[1:length(times),,6, year.k] <- apply(co.inds[,], MARGIN=1, FUN=function(x) sw_sfc[x[1],x[2]])
  clim.dat[1:length(times),,7, year.k] <- apply(co.inds[,], MARGIN=1, FUN=function(x) lw_sfc[x[1],x[2]])
  # #Check issue with radiation
  # rad= apply(co.inds[,], MARGIN=1, FUN=function(x) extract(sw_sfc, x) )
  # clim.dat[1:length(times),,6, year.k] <- apply(co.inds[,], MARGIN=1, FUN=function(x) extract(sw_sfc, x)[1,] )
  # clim.dat[1:length(times),,7, year.k] <- apply(co.inds[,], MARGIN=1, FUN=function(x) extract(lw_sfc, x)[2,] )
  
#-------

#microclimate and biophysical
#loop sites

for(site.k in 1:nrow(co.inds)){

#MICROCLIMATE
#time of sunrise and sunset
daylengths<- sapply(doys, FUN= daylength, lat= co.pts[site.k,3])
solarnoons<- sapply(doys, FUN= solar_noon, lon= co.pts[site.k,4])

t_r= solarnoons - daylengths/2
t_s= solarnoons + daylengths/2

# scale temperature from two meters
z <- 0.001 #specify distance from ground 
Tminz <- air_temp_profile(T_r = clim.dat[1:365,site.k,c(1), year.k], u_r = clim.dat[1:365,site.k,c(5), year.k], zr = 2, z0 = 0.02, z = z, T_s = clim.dat[1:365,site.k,c(4), year.k])
Tmaxz <- air_temp_profile(T_r = clim.dat[1:365,site.k,c(2), year.k], u_r = clim.dat[1:365,site.k,c(5), year.k], zr = 2, z0 = 0.02, z = z, T_s = clim.dat[1:365,site.k,c(4), year.k])
#check surface roughness

# scale wind speed
Uz <- wind_speed_profile_neutral(u_r = clim.dat[1:365,site.k,c(5), year.k], zr = 2, z0 = 0.02, z = z)

# #plot
# plot(1:365, clim.dat[1:365,site.k,c(1), year.k], col="blue", type="l")
# points(1:365, Tminz, col="blue", lty="dashed", type="l")
# points(1:365, clim.dat[1:365,site.k,c(1), year.k], col="red", type="l")
# points(1:365, Tminz, col="red", lty="dashed", type="l")
# #max and min close

#combine input
cdat<- cbind(clim.dat[1:365,site.k,c(1,2,6), year.k], t_r, t_s, 1:365)
#Tmin, Tmax, sw

#diurnal temperature variation
##use sine exp, need to fix Gamma
#Thr<- apply(cdat, MARGIN=1, FUN= function(x) diurnal_temp_variation_sineexp(T_max=x[2], T_min=x[1], t=1:24, t_r=x[4], t_s=x[5])) #, alpha = 1.52, beta = 2.00, gamma = -0.18
# CO values, https://trenchproject.github.io/TrenchR/reference/diurnal_temp_variation_sineexp.html
#use sine
Thr<- apply(cdat, MARGIN=1, FUN= function(x) diurnal_temp_variation_sine(T_max=x[2], T_min=x[1], t=1:24))
Th.dat[site.k, year.k, 1:365,]=t(Thr)

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
for(hr in 1:24){

#combine data for hours
gdat= cbind(Thr[hr,]-273.15, clim.dat[1:365,site.k,c(4), year.k]-273.15, Uz, clim.dat[1:365,site.k,6, year.k], psi[hr,])  
#fix radiation

Te.dat[site.k, year.k, 1:365, hr] <- apply(gdat, MARGIN=1, FUN= function(x) Tb_grasshopper(T_a = x[1], T_g = x[2], u = x[3], S = x[4], K_t = 1, psi = x[5], l = 0.0211, Acondfact = 0.0, abs = 0.9, r_g = 0.5))

} #end loop hours
} #end loop sites
} #end loop years, 1:length(years)

#save data
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperRoL/data/")
saveRDS(clim.dat, file = "climdat.rds")
saveRDS(Te.dat, file = "Tedat.rds")
saveRDS(Th.dat, file = "Thdat.rds")

# Restore the object
# clim.dat<- readRDS(file = "climdat.rds")
# Te.dat<- readRDS(file = "Tedat.rds")
# Th.dat<- readRDS(file = "Thdat.rds")

#==================================================
#development
#Rearing data used to estimate lower developmental temperature and growing degree days required and estimate phenology
#look up our past rearing data
#M. sanuinipes, Hilbert and Logan 1983

#Development as a proportion per day as a function of T (Kelvin)
devel.prop= function(T, alpha=0.0455, beta=8814.36, gamma=-14877.95, delta=298.81, lambda=47258.52, theta=316.695, R=1.987){
  alpha*(T/298)*exp((beta/R)*(1/298 - 1/T))/(1+exp((gamma/R)*(1/delta - 1/T))+exp((lambda/R)*(1/theta - 1/T)))}

plot(1:60, devel.prop(1:60+273.15), type="l", ylab="Development (proportion)", xlab="body temperature (C)", cex.lab=1.5)

# ENERGETICS
# Energy in:
#   Foraging time based on temperature limits for activity, but unlikely limiting
# Assimilation rate [likely limiting?]: can roughly estimate based on current assimilation rate data; Julia improving by weighing wheat grass consumption; and defectation?

#Deutsch et al. TPC
#Performance Curve Function from Deutsch et al. 2008
tpc.perf= function(T,Topt,CTmin, CTmax){
  F=T
  F[]=NA
  sigma= (Topt-CTmin)/4
  F[T<=Topt & !is.na(T)]= exp(-((T[T<=Topt & !is.na(T)]-Topt)/(2*sigma))^2) 
  F[T>Topt & !is.na(T)]= 1- ((T[T>Topt & !is.na(T)]-Topt)/(Topt-CTmax))^2
  #1% performance at and outside CT limits
  F[F<=0]<-0.01
  
  return(F)
}

#plot(1:50, tpc.perf(1:50, Topt=spec.dat[1,"PBT"], CTmin=spec.dat[1,"CTmin"], CTmax=spec.dat[1,"CTmax"]))

#   Energy use:
#   Metabolic rate, accounting for elevation due to activity?
# temperature dependence of MR (ml CO2/hr): Buckley et al. 2013 JAE, https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.12083
# converted into lipid consumed based on a 2L O2 consumed for every 1g of lipid consumed
#assume respiratory quotient of 0.7 for lipids update
#converted to energy use assuming 39 kJ g^−1 lipid

vCO2= function(M, Tb, elev_m, b0, b1, b2, b3, k=8.62*10^-5) exp(b0 +b1*log(M)+b2*(1/(k*(Tb+273.15)))+b3*elev_m)
lipid.g= function(vCO2) (0.7*vCO2/2)

plot(1:60, vCO2(M=spec.dat[1,"Massg_C1"], 1:60, elev_m=3048, b0=spec.dat[1,"b0"], b1=spec.dat[1,"b1"], b2=spec.dat[1,"b2"], b3=spec.dat[1,"b3"]), type="l", ylab="MR (ml CO2/hr)", xlab="body temperature (C)", cex.lab=1.5)

#   Assume portion of energy allocated to maintenance, reproduction
# [Needs: estimation approach; based on lipid content?
#     

# FECUNDITY: 
#   Estimate energy available for reproduction then use an energetic estimate of reproduction to estimate fecundity
# [Needs: energetic cost of reproduction
#   Options: bomb calorimetry of eggs; mass (Levy and Nufio paper) and composition of eggs [preferred?]
#   Infer from eggpods laid in field, but seems better to use data to validate estimates]

#Levy and Nufio
#Mass of eggs
#Number of eggs
#https://github.com/lbuckley/HopperPhenology/blob/master/HopperClinesForGenetics.R
setwd("/Volumes/GoogleDrive/My Drive/GrasshopperData/TPCfield/data/")
#egg= read.csv("LevyNufioData.csv")
#size= read.csv("Levy_Male&FemaleBodySizesGradient.csv")
repro= read.csv("Levy_FemaleGradientDataGrasshopper.csv")

repro2= as.data.frame(repro) %>% group_by(Site,Elevation, Species) %>% summarise(clutch=mean(Clutch.Size, na.rm=TRUE), egg.mass=mean(Mean.Egg.Mass, na.rm=TRUE), Novarioles=mean(Number.Ovarioles, na.rm=TRUE), NFovarioles=mean(Number.Functioing.Ovarioles, na.rm=TRUE), PFovarioles=mean(Proportion.Functioning.Ovarioles, na.rm=TRUE), clutch.mass=mean(clutch.weight.g, na.rm=TRUE) )

#accumulate energy until egg mass and clutch size reached

assim=0.215 #assimilation rate, Fielding 2004
#Assume 20% efficiency in converting energy into eggs, Fielding 2004

#-----------
# SURVIVAL: Use thermal tolerance and stress response data to estimate daily survival. We will model survival as the product of survival to adulthood and daily survival. 

#Functions
surv<- function(T, CTmin, CTmax, td=4.34){ 
  #20 to 80% of CT range
  CTmin1= CTmin+(CTmax-CTmin)*0.2#*0.1
  CTmax1= CTmin+(CTmax-CTmin)*0.8#*0.9
  
  s1= ifelse(T<CTmax1, s<-1, s<- exp(-(T-CTmax1)/td) )
  s2= ifelse(T>CTmin1, s<-1, s<- exp(-(CTmin1-T)/td) )
  s= s1*s2
  return( s*0.8 )
}

#plot(0:70, surv(0:70, 10, 60), type="l", ylab="survival (%)", xlab="body temperature (C)", cex.lab=1.5)
#points(c(10,60),c(0.8,0.8))

#define geometric mean
geo_mean <- function(data) {
  log_data <- log(data)
  gm <- exp(mean(log_data[is.finite(log_data)]))
  return(gm)
}

#=================================
#RUN MODEL

#array for fitness data
fit.dat<- array(data = NA, dim = c(nrow(co.inds),length(years),4), dimnames = NULL)
#3rd dimension in fecundity, survival, fitness, phenology

##loop grid cells
for(grid.k in 1:nrow(co.pts) ){
#grid.k=1

##loop years
for(yr.k in 1:length(years)){    #1:length(years)
#yr.k<- 1
year= years[yr.k]

##loop phenotypes?

#estimate development (proportion per day) starting March 1 through August, doy=60:243

#make time series to estimate development
doys.dev<- 1:243
T.series= cbind(doys.dev, Th.dat[grid.k, yr.k, 1:243,])
colnames(T.series)=c("doy", 1:24)
#to long format
Tl<- as.data.frame(T.series) %>%
  gather("hr", "Thr", 2:ncol(T.series))
Tl$hr<- as.numeric(Tl$hr)
#sort by doy then hr
Tl<- Tl[order(Tl$doy,Tl$hr),]

#estimate development
Tl$dev= devel.prop(T= Tl$Thr)/24
Tl$dev.sum= cumsum(replace_na(Tl$dev, 0))
#Find phenology
doy.ad= Tl$doy[min(which(Tl$dev.sum>1))]
fit.dat[grid.k, yr.k,4]<- doy.ad
#add in delay after adult development until egg laying,  https://www.jstor.org/stable/4219411

#diapause and embryogenesis, Hilbert 1985

#model survival and longevity
#Fielding 2004 https://doi.org/10.1016/j.ecolmodel.2003.10.014
# Daily mortality rate for eggs, juveniles, adults 0.01

#Carter et al 1998 https://doi.org/10.1093/ee/27.4.892
#Adult longevity functions, 
#long lab survival: Tatat et al. 1997, https://link.springer.com/article/10.1007/s004420050246

#------
if(!is.na(doy.ad) & doy.ad<240){

#start accumulating energy to translate into eggs
Te.series= cbind(doy.ad:243, Te.dat[grid.k, yr.k, doy.ad:243,])
colnames(Te.series)=c("doy", 1:24)

#to long format
Tel<- as.data.frame(Te.series) %>%
  gather("hr", "Te", 2:ncol(Te.series))
#limit Te to 10C over CTmax
Tel$Te[Tel$Te> (spec.dat[1,"CTmax"]+5)]<- spec.dat[1,"CTmax"] +5

#estimate survival in response to extremes
Tel$surv= surv(Tel$Te, CTmin=spec.dat[1,"CTmin"], CTmax=spec.dat[1,"CTmax"])

#energy in #NEED TO UPDATE
Tel$Ein<- tpc.perf(Tel$Te,Topt=spec.dat[1,"PBT"], CTmin=spec.dat[1,"CTmin"], CTmax=spec.dat[1,"CTmax"])
#doesn't yet account for assimilation: assim*

#energy use
Tel$vCO2<- vCO2(M=spec.dat[1,"Massg_C1"], Tel$Te, elev_m=3048, b0=spec.dat[1,"b0"], b1=spec.dat[1,"b1"], b2=spec.dat[1,"b2"], b3=spec.dat[1,"b3"])
Tel$lipid.g= lipid.g(Tel$vCO2)
#account for activity?

#net energy to eggs
Tel$Enet<- Tel$Ein - Tel$lipid.g

#account for survival: 0.01 daily mortality rate for eggs, juveniles, and adults Fielding 2004
Tel$dsurv= 0.99^Tel$doy
Tel$Enet.surv= Tel$Enet*Tel$dsurv

#cumulative energetics
#replace NAs with zeros
Tel$Enet.surv[is.na(Tel$Enet.surv)]<-0
Tel$Ecumsum= cumsum(Tel$Enet.surv)*0.2 #assume 20% conversion to eggs #UPDATE

#estimate fecundity discounting by survival
clutch.mass<- as.numeric(repro2[repro2$Site=="C1" & repro2$Species=="sanguinipes","clutch.mass"])

n.eggs<- as.numeric(floor(Tel$Ecumsum[nrow(Tel)]/clutch.mass)*repro2[repro2$Site=="C1" & repro2$Species=="sanguinipes","clutch"])
#limit to max number eggs #Approximate value, https://doi.org/10.4039/Ent98617-6
if(n.eggs>500) n.eggs<-500

#estimate fitness as (sum of fecundity)(product of survival)
#https://github.com/lbuckley/FitnessJEB
fit.dat[grid.k, yr.k,1]<- n.eggs
fit.dat[grid.k, yr.k,2]<- geo_mean(na.omit(Tel$surv))
fit.dat[grid.k, yr.k,3]<- n.eggs * geo_mean(na.omit(Tel$surv))

} #end check can develop
} #end loop years
} #end loop grid cells

#------------
#plot output

#add elevation to pts
library(elevatr)
ll_prj <- "EPSG:4326"
co.pts<- as.data.frame(co.pts)
examp_df<- data.frame(x=co.pts[,4], y=co.pts[,3])
elevs<- get_elev_point(examp_df, prj = ll_prj, src = "aws")
co.pts$elev.m <- elevs$elevation

# number eggs
fit1= cbind(co.pts, fit.dat[,,1])
colnames(fit1)[6:ncol(fit1)]<- years
#to long format
fit1<- as.data.frame(fit1) %>%
  gather("year", "value", 6:ncol(fit1))
fit1$comp<- "eggs"

#survival
fit2= cbind(co.pts, fit.dat[,,2])
colnames(fit2)[6:ncol(fit2)]<- years
#to long format
fit2<- as.data.frame(fit2) %>%
  gather("year", "value", 6:ncol(fit2))
fit2$comp<- "surv"

#fitness
fit3= cbind(co.pts, fit.dat[,,3])
colnames(fit3)[6:ncol(fit3)]<- years
#to long format
fit3<- as.data.frame(fit3) %>%
  gather("year", "value", 6:ncol(fit3))
fit3$comp<- "fit"

#phenology
fit4= cbind(co.pts, fit.dat[,,4])
colnames(fit4)[6:ncol(fit4)]<- years
#to long format
fit4<- as.data.frame(fit4) %>%
  gather("year", "value", 6:ncol(fit4))
fit4$comp<- "phen"

fit<- rbind(fit1, fit2, fit3, fit4)

plot.co= ggplot(data=fit, aes(x=lons, y=lats))+
  geom_point(aes(color=value), size=4)+
  facet_grid(comp~year)

ggplot(data=fit4, aes(x=lons, y=lats))+
  geom_point(aes(color=value), size=5)+
  facet_wrap(.~year, ncol=4)

# #doesn't work due to uneven grid
# ggplot(data=fit4, aes(x=lons, y=lats))+
#   geom_raster(aes(fill=value), interpolate = TRUE)+
#   facet_wrap(.~year, ncol=4)

#plot by elevation
ggplot(data=fit, aes(x=elev.m, y=value, color=year, group=year))+
  geom_point()+
  geom_smooth()+
  facet_wrap(.~comp, scales="free_y")




