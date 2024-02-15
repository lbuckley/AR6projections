library(ncdf4)
library(raster)
library(terra)
library(fields)
#install.packages("PCICt")
library(PCICt)
#install.packages("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/NicheMapR_3.2.1.tgz", repos = NULL, type = .Platform$pkgType)
library(NicheMapR)
library(TrenchR)

#Read species data
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/FitnessContrib_JEB/data/Data/Grasshopper/")
spec.dat=read.csv("SpecData.csv")
site.list= c("Eldorado","A1","B1","C1")
#specs= c("clav","pell","dodg","sang")
#spec.dat=spec.dat[match(spec.dat$SpecID, specs),]

# Jumping TPC parameters
tpc.dat= read.csv("JumpTPCparams.csv")

#Calculate mass and length, m and g? 
dat$mass_dodg= 1-1.64*10^-4*3048 #*dat$Elev
dat$L_dodg= exp(3.33*0.247*log(dat$mass_dodg))

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
for (year.k in 7:length(years) ){ #length(years)
  
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
  #clim.dat[1:length(times),,6, year.k] <- apply(co.inds[,], MARGIN=1, FUN=function(x) sw_sfc[x[1],x[2]])
  #clim.dat[1:length(times),,7, year.k] <- apply(co.inds[,], MARGIN=1, FUN=function(x) lw_sfc[x[1],x[2]])
  #Check issue with radiation
  rad= apply(co.inds[,], MARGIN=1, FUN=function(x) extract(sw_sfc, x) )
  clim.dat[1:length(times),,6, year.k] <- apply(co.inds[,], MARGIN=1, FUN=function(x) extract(sw_sfc, x)[1,] )
  clim.dat[1:length(times),,7, year.k] <- apply(co.inds[,], MARGIN=1, FUN=function(x) extract(lw_sfc, x)[2,] )
  
  #slow across sites
  #clim.dat[1:length(times),1:1000,1]<- apply(co.inds[1:1000,], MARGIN=1, FUN=function(x) t2min[x[1],x[2]])
  
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
clim.dat
Te.dat
Th.dat

# Save an object to a file
saveRDS(clim.dat, file = "climdat.rds")
saveRDS(Te.dat, file = "Tedat.rds")
saveRDS(Th.dat, file = "Thdat.rds")

# # Restore the object
# clim.dat<- readRDS(file = "climdat.rds")
# Te.dat<- readRDS(file = "Tedat.rds")
# Th.dat<- readRDS(file = "Thdat.rds")

#==================================================
#development
#Rearing data used to estimate lower developmental temperature and growing degree days required and estimate phenology
#look up our past rearing data
#M. sanuinipes, Hilbert and Logan 1983

#T is in Kelvin
dr= function(T, alpha=0.0455, beta=8814.36, gamma=-14877.95, delta=298.81, lambda=47258.52, theta=316.695, R=1.987){
  alpha*(T/298)*exp((beta/R)*(1/298 - 1/T))/(1+exp((gamma/R)*(1/delta - 1/T))+exp((lambda/R)*(1/theta - 1/T)))}

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

#plot(1:50, tpc.perf(1:50, Topt=30, CTmin=10, CTmax=40))

TPC.gausgomp= function(T, To, rho=0.9, sigma, Fmax) Fmax*exp(-exp(rho*(T-To)-6)-sigma*(T-To)^2) # rho and sigma determine, the thermal sensitivity of feeding at temperatures above and below Topt, respectively

#   Energy use:
#   Metabolic rate, accounting for elevation due to activity?
# temperature dependence of MR (ml CO2/hr): Buckley et al. 2013 JAE, https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.12083
# converted into lipid consumed based on a 2L O2 consumed for every 1g of lipid consumed
#assume respiratory quotient of 0.7 for lipids update
#converted to energy use assuming 39 kJ g^−1 lipid

bs.ms<- c(0.90, 0.48, 8.68 * 10^{-5})
bs.mb<- c(0.98, 0.48, 1.24 * 10^{-4})

vCO2= function(M, Tb, elev_m, b1, b2, b3, k=8.62*10^-5) exp(b1*log(M)+b2*(1/(k*Tb))+b3*elev_m)
lipid.g= function(vCO) (0.7*vCo2/2)

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
setwd("/Volumes/GoogleDrive/My Drive/AlexanderResurvey/TPCfield/data")
egg= read.csv("LevyNufioData.csv")
size= read.csv("Levy_Male&FemaleBodySizesGradient.csv")
repro= read.csv("Levy_FemaleGradientDataGrasshopper.csv")

repro2= repro %>% group_by(Site,Elevation, Species) %>% summarise(clutch=mean(Clutch.Size, na.rm=TRUE), clutch.N= length(na.omit(Clutch.Size)),clutch.sd=sd(Clutch.Size, na.rm=TRUE), clutch.se=clutch.sd/sqrt(clutch.N),
                                                                  egg.mass=mean(Mean.Egg.Mass, na.rm=TRUE), egg.mass.N= length(na.omit(Mean.Egg.Mass)),egg.mass.sd=sd(Mean.Egg.Mass, na.rm=TRUE), egg.mass.se=egg.mass.sd/sqrt(egg.mass.N),
                                                                  Novarioles=mean(Number.Ovarioles, na.rm=TRUE), Novarioles.N= length(na.omit(Number.Ovarioles)),Novarioles.sd=sd(Number.Ovarioles, na.rm=TRUE), Novarioles.se=Novarioles.sd/sqrt(Novarioles.N),
                                                                  NFovarioles=mean(Number.Functioing.Ovarioles, na.rm=TRUE), NFovarioles.N= length(na.omit(Number.Functioing.Ovarioles)),NFovarioles.sd=sd(Number.Functioing.Ovarioles, na.rm=TRUE), NFovarioles.se=NFovarioles.sd/sqrt(NFovarioles.N),
                                                                  PFovarioles=mean(Proportion.Functioning.Ovarioles, na.rm=TRUE), PFovarioles.N= length(na.omit(Proportion.Functioning.Ovarioles)),PFovarioles.sd=sd(Proportion.Functioning.Ovarioles, na.rm=TRUE), PFovarioles.se=PFovarioles.sd/sqrt(PFovarioles.N), 
                                                                  clutch.mass=mean(clutch.weight.g, na.rm=TRUE), clutch.mass.N= length(na.omit(clutch.weight.g)),clutch.mass.sd=sd(clutch.weight.g, na.rm=TRUE), clutch.mass.se=clutch.mass.sd/sqrt(clutch.mass.N) )

#--------------
sp.col<-c("A. clavatus"="yellow",  "M. boulderensis"="#fecc5c", "C. pellucida"="#fd8d3c","M. sanguinipes"="#e31a1c")  
repro2$Species= factor(repro2$Species, levels=c("A. clavatus","M. boulderensis","C. pellucida","M. sanguinipes"), ordered=TRUE)

ov.plot=ggplot(data=repro2, aes(x=Elevation, y = Novarioles,color=Species))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+ 
  geom_errorbar(aes(ymin=Novarioles-Novarioles.se, ymax=Novarioles+Novarioles.se), width=.1)+
  scale_colour_manual(values = sp.col)

clutch.plot=ggplot(data=repro2, aes(x=Elevation, y = clutch,color=Species))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+ 
  geom_errorbar(aes(ymin=clutch-clutch.se, ymax=clutch+clutch.se), width=.1)+
  scale_colour_manual(values = sp.col)

egg.mass.plot=ggplot(data=repro2, aes(x=Elevation, y = egg.mass*1000,color=Species))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+ xlim(1550,3530)+
  geom_errorbar(aes(ymin=egg.mass-egg.mass.se, ymax=egg.mass+egg.mass.se), width=.1)+ylab("Egg mass (mg)")+ theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  scale_colour_manual(values = sp.col)

clutch.mass.plot=ggplot(data=repro2, aes(x=Elevation, y = clutch.mass,color=Species))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+ xlim(1550,3530)+
  geom_errorbar(aes(ymin=clutch.mass-clutch.mass.se, ymax=clutch.mass+clutch.mass.se), width=.1)+ylab("Clutch mass (g)")+xlab("Elevation (m)")+
  scale_colour_manual(values = sp.col)

Fov.plot=ggplot(data=repro2, aes(x=Elevation, y = NFovarioles,color=Species))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+ xlim(1550,3530)+
  geom_errorbar(aes(ymin=NFovarioles-NFovarioles.se, ymax=NFovarioles+NFovarioles.se), width=.1)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())+ylab("N functional ovarioles")+
  scale_colour_manual(values = sp.col)

PFov.plot=ggplot(data=repro2, aes(x=Elevation, y = PFovarioles,color=Species))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+ 
  geom_errorbar(aes(ymin=PFovarioles-PFovarioles.se, ymax=PFovarioles+PFovarioles.se), width=.1)+
  scale_colour_manual(values = sp.col)


assim=0.215 #assimilation rate, Fielding 2004
#Assume 20% efficiency in converting eggs  into Fielding 2004

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

#plot(0:70, surv(0:70, 10, 60), type="l")
#points(c(10,60),c(0.8,0.8))

#=================================
#FECUNDITY
#Based on CTmin, Topt, CTmax
Te.dat[ind.k,spec.k,4,]<-tpc.perf( ts, Topt=spec.dat[spec.k,"PBT"], CTmin=spec.dat[spec.k,"CTmin"], CTmax=spec.dat[spec.k,"CTmax"])

#Based on hopping TPC
# if(!specs[spec.k]=="clav"){ #NOT CLAVATUS
# ind.elev= dat[ind.k,"Elev"]
# if(ind.elev==1708)ind.elev<-2195
# tpc= tpc.dat[which(tpc.dat$spec==specs[spec.k] & tpc.dat$elev_m==2591),]  ##use mid elevation
# 
# Te.dat[ind.k,spec.k,6,]<-TPC.gausgomp(ts, To=tpc$To, rho=0.7, sigma=tpc$sigma, Fmax=tpc$Pmax)
#} #check for clavatus
#Te.dat[ind.k,spec.k,6,]<-predict(lo,ts)
perf= rezende_2019(ts, q10=2.27, a=0.109, b=9.02, c=0.00116)
perf[which(perf<0)]=0
Te.dat[ind.k,spec.k,6,]<-perf

#Survival  
Te.dat[ind.k,spec.k,5,]<- surv( ts, CTmin= spec.dat[spec.k,"CTmin"], CTmax=spec.dat[spec.k,"CTmax"])  #HEAT STRESS ONLY spec.dat[spec.k,"CTmin"] -20


#-------------------------------
# FITNESS: We will multiply these additive (fecundity) and multiplicative (survival) components to estimate fitness and the fitness estimates will provide estimates of selection on phenotypes / genotypes. The model could potentially be extended to include a quantitative genetic model of evolution


#Fielding 2004 https://doi.org/10.1016/j.ecolmodel.2003.10.014
# *** IBM based on resources
# somatic mass; reproductive mass
# Assimilation rate (food quality) at beginning of season, and on day t 0.215
# Consumption rate (mg food mg−1 grasshopper) 0.530
# Respiratory loss (mg mg−1 grasshopper) 0.034
# Growth rate = (ca − l) 0.080
# Daily mortality rate for eggs, juveniles, adults 0.01
# After individuals attain their target size, they be- come adults and assimilated biomass is no longer allo- cated towards growth, but towards reproduction.
# The same function is used for accumulation of reproductive biomass as somatic growth, except the rate of weight gain was divided by a factor of 5 to account for the greater concentration of energy and protein in the eggs and the extra biomass associated with production of a pod. 
# non-selective background mortal- ity of 0.01 per day was instituted for juveniles and adults. Eggs were subject to a constant mortality rate (default = 0.015 per day). 

#Carter et al 1998 https://doi.org/10.1093/ee/27.4.892
# *** M sanguinipes population model
# A maximum rate of egg production per female per degree day was estimated from data in Pickford (1958), pfadt (1949), Smith (1966), and Smith (1968). 
# egg production factor
# mortality associated with precipitation
# development functions #USE?
# Adult longevity functions #USE?
# rate of embryogenesis as a function of temperature
# diapause development rate

#Hilbert and Logan 1983. development https://doi.org/10.1093/ee/12.1.1
# *** temperature dependence of development for Melanoplus sanguinipes
# development data #USE

#Hilbert and Logan 1983. https://doi.org/10.1016/B978-0-444-42179-1.50041-9 
#     0.69 eggs per degree day greater than 20C

#Berry et al. 1993. Object-oriented Simulation Model of Rangeland Grasshopper Population Dynamics 
#growth parameter food to grasshopper and egg mass
# ***
# cGrowth: parameter, conversion factor used to convert food to grasshopper mass (0.05)
# cEggs: parameter, conversion factor used to convert food to grasshopper egg mass (0.1)

#Rodell 1977. A GRASSHOPPER MODEL FOR A GRASSLAND ECOSYSTEM. https://esajournals.onlinelibrary.wiley.com/doi/pdf/10.2307/1935600
#*** rate of oviposition: proportion of female bodymass going to eggs per day
# egg production function
# low temp and egg viability function
# Other literature values include assimilation efficiencies of 0.274 (Smalley 1960), 0.31-0.48 (Husain et al. 1946), 0.16-0.20 (Hussain 1972), 0.35-0.78 (Davey 1954), 0.161-0.486 (Gyllenberg 1969), and 0.226-0.322 (Mitchell 1973).
#------------------------
