#read climate data
#solar, air temp, ground temp, wind

#microclimate model

#biophysical model 
#Lactin and Johnson 1997  https://doi.org/10.1111/j.1365-3032.1997.tb01150.x
#melanization, using Spectroradiometer like that at UC Berkeley, Fielding and Defoliant https://doi.org/10.1665/1082-6467(2005)14[107:DATMOF]2.0.CO;2

#development
#Rearing data used to estimate lower developmental temperature and growing degree days required and estimate phenology
#look up our past rearing data
#Alaska temp dependence of development: Fielding 2004 https://doi.org/10.1603/0046-225X-33.6.1513


# ENERGETICS
# Energy in:
#   Foraging time based on temperature limits for activity, but unlikely limiting
# Assimilation rate [likely limiting?]: can roughly estimate based on current assimilation rate data; Julia improving by weighing wheat grass consumption; and defectation?
#   
#   Energy use:
#   Metabolic rate, accounting for elevation due to activity?
#   Assume portion of energy allocated to maintenance, reproduction
# [Needs: estimation approach; based on lipid content?
#     

# FECUNDITY: 
#   Estimate energy available for reproduction then use an energetic estimate of reproduction to estimate fecundity
# [Needs: energetic cost of reproduction
#   Options: bomb calorimetry of eggs; mass (Levy and Nufio paper) and composition of eggs [preferred?]
#   Infer from eggpods laid in field, but seems better to use data to validate estimates]

# SURVIVAL: Use thermal tolerance and stress response data to estimate daily survival. We will model survival as the product of survival to adulthood and daily survival. 

# FITNESS: We will multiply these additive (fecundity) and multiplicative (survival) components to estimate fitness and the fitness estimates will provide estimates of selection on phenotypes / genotypes. The model could potentially be extended to include a quantitative genetic model of evolution

#Past models
#Branson 2003 https://doi.org/10.4039/n02-061
#Cage experiment with Melanoplus sanguinipes in Montana
#Regressions of eggs laid as a function of digestible biomass (but weird units)

#Branson 2008 https://doi.org/10.1665/1082-6467-17.2.259
#reproductive correlates in Melanoplus sanguinipes

#Branson 2004 https://doi.org/10.1665/1082-6467(2004)013[0239:RIONAA]2.0.CO;2
#Reproductive allocation in experiments

#Branson 2022 https://doi.org/10.3390/geographies2010003, Wyoming abundance sampling

# *** De Sousa Santos and Begon 1987 https://www.jstor.org/stable/2389423 
#Estimates of reproductive energetics for UK grasshopper 

#I. Filin and O. Ovadia . 2007 
# *** Energetic population model for M. femurrubrum https://doi.org/10.1086/522091 

#Hatle et al 2003 https://doi.org/10.1093/icb/43.5.635
# *** Reproductive allocation through time

#Jonas et al 2015 https://doi.org/10.1016/j.rama.2014.12.011
#population dynamic regressions

#Fielding 2004 https://doi.org/10.1016/j.ecolmodel.2003.10.014
# *** IBM based on resources

#Ackman and Whitman 2008. https://doi.org/10.1665/1082-6467-17.2.249
# size and fecundity

#Carter et al 1998 https://doi.org/10.1093/ee/27.4.892
# *** M sanguinipes population model
# 0.69 eggs produced per degree day greater than 20C
# egg production factor
# mortality associated with precipitation

#Hilbert and Logan 1983. development https://doi.org/10.1093/ee/12.1.1
# *** temperature dependence of development for Melanoplus sanguinipes

#Hilbert and Logan 1983. https://doi.org/10.1016/B978-0-444-42179-1.50041-9 
# Get book

#Berry et al. 1993. Object-oriented Simulation Model of Rangeland Grasshopper Population Dynamics 
#growth parameter food to grasshopper and egg mass
# ***

#Rodell 1977. A GRASSHOPPER MODEL FOR A GRASSLAND ECOSYSTEM. https://esajournals.onlinelibrary.wiley.com/doi/pdf/10.2307/1935600
#*** proportion of female bodymass going to eggs per day
# egg production function
# low temp and egg viability

#------------------------
### LOAD DATA
setwd("/Volumes/GoogleDrive/My\ Drive/Buckley/Work/ColiasBiogeog/Data/")
pts.sel<- read.csv("pts.sel.csv")
Ts_sh_min<- read.csv("Ts_sh_min.csv")
Ts_sh_max<- read.csv("Ts_sh_max.csv")

#load objects
solar <- readRDS("solar.rds")
Ts_sun <- readRDS("Ts_sun.rds")
Ts_sh <- readRDS("Ts_sh.rds")
wind <- readRDS("wind.rds")
zenith <- readRDS("zenith.rds")

### SET UP DATA STRUCTURES
#Make array to store data
#Lambda= survival*fecundity
#Matrix of lambdas
#dims: stats, year, Lambda 
Lambda<-array(NA, dim=c(length(years),nrow(pts.sel),length(seq(0.4,0.7,0.05)),4)) #Last dimension is Lambda, FAT,Egg Viability
dimnames(Lambda)[[1]]<-years

#Matrix for pupual temps
pup.temps<-array(NA, dim=c(12, length(years),nrow(pts.sel),3)) #3 generations 
#Add names
dimnames(pup.temps)[[1]]= c("stat","yr","gen","Jlarv", "Jpup","Jadult","Tlarv","Tpup","Tad","Tlarv_fixed","Tpup_fixed","Tad_fixed") 

#Te
dayk= 152:243
Te.mat.all= array(data=NA, dim=c(nrow(pts.sel), length(6:20),7, length(dayk) ) )

aseq= seq(0.4,0.7,0.05)

#===================================================================
#LOOP YEARS

for(yr.k in 1:length(years) ){
  print(yr.k)
  
  inds= which(time.mat[,2]==years[yr.k])
  
  tmax.yr= tmax[, , inds, proj.k]  
  tmin.yr= tmin[, , inds, proj.k] 
  
  #============================================================================
  #CALCULATE DEVELOPMENT TIMING AND TEMPS
  
  #Temps at plant height
  lon.inds= pts.sel[, "lon.ind"]
  lat.inds= pts.sel[, "lat.ind"]
  lonlat.inds <- cbind(lon.inds,lat.inds)
  
  Ta_plant_min<- foreach(d=60:273, .combine='cbind')  %do% {
    T_mat= tmin.yr[, ,d]
    T_mat= cbind(T_mat[lonlat.inds], Ts_sh_min[mo[d]-2, ])
    apply(T_mat,MARGIN=1, FUN=air_temp_at_height_z_mat, z_0=z_0_1, z_r=2, z=z_0_1)}
  
  Ta_plant_max<- foreach(d=60:273, .combine='cbind')  %do% {
    T_mat= tmax.yr[, ,d]
    T_mat= cbind(T_mat[lonlat.inds], Ts_sh_max[mo[d]-2, ])
    apply(T_mat,MARGIN=1, FUN=air_temp_at_height_z_mat, z_0=z_0_1, z_r=2, z=z_0_1)}
  
  #transpose
  Ta_plant_min= t(Ta_plant_min)
  Ta_plant_max= t(Ta_plant_max)
  
  #----------------------------
  # ESTIMATE DEVELOPMENTAL TIMING
  
  DevZeros= c(9.2176, 11.5, 9.7) #4th and 5th, larval, pupal
  GddReqs= c(117.06, 270.39 ,101.9) 
  
  #Calc GDDs
  GDDs_45<- foreach(d=(60:273)-59, .combine='cbind')  %do% {
    T_mat= cbind(Ta_plant_min[d,], Ta_plant_max[d,])
    apply(T_mat,MARGIN=1, FUN=degree.days.mat, LDT=DevZeros[1])}
  
  GDDs_l<- foreach(d=(60:273)-59, .combine='cbind')  %do% {
    T_mat= cbind(Ta_plant_min[d,], Ta_plant_max[d,])
    apply(T_mat,MARGIN=1, FUN=degree.days.mat, LDT=DevZeros[2])}
  
  GDDs_p<- foreach(d=(60:273)-59, .combine='cbind')  %do% {
    T_mat= cbind(Ta_plant_min[d,], Ta_plant_max[d,])
    apply(T_mat,MARGIN=1, FUN=degree.days.mat, LDT=DevZeros[3])}
  
  GDDs_45= t(GDDs_45)
  GDDs_l= t(GDDs_l)
  GDDs_p= t(GDDs_p)
  
  #-----------------------------------
  #Calculate development timing
  Js=60:273
  
  for(cell.k in 1:nrow(pts.sel)){  #  for(cell.k in 1:nrow(pts.sel))
    
    Ta_plant_mean= rowMeans(cbind(Ta_plant_min[, cell.k], Ta_plant_max[, cell.k]))
    
    for(gen.k in 1:3){
      
      #Assume 7 days from eclosion to eggs laid
      #Hatching ~5days (~70hrs, based on Heidi's heat shock data, Jessica's development data doesn't seem to have hatching time)
      Jlarv= ifelse(gen.k>1, Jadult+12, Js[which.max(GDDs_45[,cell.k]>0)] )  
      if(Jlarv>max(Js)) Jlarv=max(Js)
      
      ##TO PUPATION
      check=ifelse(gen.k==1, gdds<-GDDs_45[,cell.k],gdds<-GDDs_l[,cell.k])
      
      Jpup= Jlarv + which.max(cumsum(gdds[which(Js==Jlarv):length(Js)])> ifelse(gen.k==1, GddReqs[1],GddReqs[2])  )
      if(Jpup>max(Js) | length(Jpup)==0) Jpup=max(Js) 
      
      #PUPATION
      gdds<-GDDs_p[,cell.k]
      Jadult= Jpup + which.max(cumsum(gdds[which(Js==Jpup):length(Js)])> GddReqs[3])
      if(Jadult>max(Js) | length(Jadult)==0) Jadult=max(Js)
      
      #----------------------
      #Calculate temps
      Tlarv= mean(Ta_plant_mean[Js %in% Jlarv:Jpup], na.rm=TRUE)
      Tpup= mean(Ta_plant_mean[Js %in% Jpup:Jadult], na.rm=TRUE)
      Tad= mean(Ta_plant_mean[Js %in% Jadult:(Jadult+7)], na.rm=TRUE)
      ### ADULT TEMP IS AIR TEMP
      #Check if more than 5 NAs
      
      #Write data in array
      pup.temps[3:9,yr.k, cell.k, gen.k]=c(gen.k,Jlarv,Jpup,Jadult,Tlarv,Tpup,Tad)
      
    } #end loop generation
    
  } #end loop cells
  
  #=========================================================================
  #Run Te calculations
  
  #Set up climate data
  hr.dec= 1:24
  
  tmax.yr= tmax[, , inds[152:243], proj.k]  
  tmin.yr= tmin[, , inds[152:243], proj.k]  
  
  clim= as.data.frame(152:243)
  names(clim)="J"
  clim$month=NA
  clim[which(clim$J<182),"month"]=6
  clim[which(clim$J>181 & clim$J<213),"month"]=7
  clim[which(clim$J>212),"month"]=8
  
  #estimate daylength
  Trise.set= suncalc(clim$J, Lat = pts.sel[1,"lat"], Long = pts.sel[1,"lon"])
  clim$set= Trise.set$sunset
  clim$rise= Trise.set$sunrise
  
  #TEMP
  Thr<- foreach(cell.k=1:nrow(pts.sel) ) %do% {
    t(apply( cbind( tmax.yr[pts.sel[cell.k, "lon.ind"], pts.sel[cell.k, "lat.ind"],], tmin.yr[pts.sel[cell.k, "lon.ind"], pts.sel[cell.k, "lat.ind"],], clim[,c("rise","set")]) , FUN=Thours.mat, MARGIN=1))}
  #Collapse list into array
  Thr= array(unlist(Thr), dim = c(nrow(Thr[[1]]), ncol(Thr[[1]]), length(Thr)))
  
  #RADIATION
  Rhr= foreach(cell.k=1:nrow(pts.sel)) %do% {
    t(apply(rbind(matrix(solar[cell.k,,1],nrow = 30,ncol = 24, byrow=TRUE), matrix(solar[cell.k,,2],nrow = 31,ncol = 24, byrow=TRUE), matrix(solar[cell.k,,3],nrow = 31,ncol = 24, byrow=TRUE) ), FUN=Rad.mat, MARGIN=1)) } 
  #Collapse list into array
  Rhr= array(unlist(Rhr), dim = c(nrow(Rhr[[1]]), ncol(Rhr[[1]]), length(Rhr)))
  #columns 1:24 are direct, 25:28 are diffuse
  #reflected is direct * albedo of 0.7
  
  #--------------------------------------------------
  #Calculate Te #cell.k: length(grid.sel)
  
  for(hr.k in 1:15 ){
    
    Thr.d= Thr[,hr.k,]
    Ts_sun.d= Ts_sun[,hr.k,]
    Ts_sh.d= Ts_sh[,hr.k,]
    wind.d= wind[,hr.k,]
    Rhr.d= Rhr[,c(hr.k,hr.k+24),]
    zenith.d= zenith[,hr.k,]
    
    #combine data
    Te.dat=abind(t(Thr.d), Ts_sun.d[,(clim$month-5)], Ts_sh.d[,(clim$month-5)], wind.d[,(clim$month-5)], t(Rhr.d[,1,]), t(Rhr.d[,2,]),zenith.d[,(clim$month-5)],along=3)
    
    for(a.k in 1:7){
      
      a= aseq[a.k]
      Te.mat.all[,hr.k,a.k,]<-  apply(Te.dat,MARGIN=c(1,2), FUN=biophys.var_sh.mat, D, delta, a)  
      
    } # end loop across absorptivities
  } # end loop across hours
  
  #=======================================
  #DEMOGRAPHY
  
  #Flight probability
  FlightProb<- foreach(a= 1:dim(Te.mat.all)[[3]] )  %:% foreach(hr=1:15) %:% foreach(d=1:dim(Te.mat.all)[[4]], .combine='cbind') %do% {  out=sapply(Te.mat.all[,hr, a,d], FUN=fl.ph) }
  #Back to array
  FlightProb.all= array(unlist(FlightProb), dim = c(nrow(FlightProb[[1]][[1]]), ncol(FlightProb[[1]][[1]]), length(FlightProb[[1]]), length(FlightProb)))
  
  #Egg viability
  EggViab<- foreach(a= 1:dim(Te.mat.all)[[3]] )  %:% foreach(hr=1:15) %:% foreach(d=1:dim(Te.mat.all)[[4]], .combine='cbind') %do% {  out=sapply(Te.mat.all[,hr, a,d], FUN=egg.viab) }
  #Back to array        
  EggViab.all= array(unlist(EggViab), dim = c(nrow(EggViab[[1]][[1]]), ncol(EggViab[[1]][[1]]), length(EggViab[[1]]), length(EggViab)))
  
  #-------------------------------------------
  #DEMOGRAPHY
  
  EV1=matrix(NA,2, nrow(SpecDat))
  
  for(gen.k in 1:3 ){ #loop generation
    
    #find cells that can complete generation
    Jfls=pup.temps["Jadult",yr.k,, gen.k] 
    cell.inds= which(Jfls<244 & Jfls>151)
    
    for(cell.k in cell.inds ){ #loop cells
      # if(cell.k/100==round(cell.k/100)) print(cell.k)
      
      for(abs.k in 1:dim(Te.mat.all)[3] ){ #loop absorptivity
        
        Te.mat= Te.mat.all[cell.k,,abs.k,]
        
        #get flight dates
        Jfl=pup.temps["Jadult",yr.k, cell.k, gen.k]   
        
        #average over hours
        FAT= rowSums(FlightProb.all[cell.k,,,abs.k], na.rm=TRUE)
        EggViab= apply(EggViab.all[cell.k,,,abs.k], FUN=geo_mean, MARGIN=1) #Egg viability GEOMETRIC MEAN ACROSS HOURS
        Temps= colMeans(Te.mat[,], na.rm=TRUE)
        
        ##CALCULATE EGG VIABILITY OVER 5 DAY PERIOD (GEOMETRIC MEAN ACROSS HOURS)
        #sample flight day from truncated normal distribution
        Nind=1000 #changed from 100
        f.low= max(Jfl-7,min(clim$J)+2)
        f.up= min(Jfl+7,max(clim$J)-2)
        
        flightday= round(rtruncnorm(Nind, a=f.low, b=f.up, mean = Jfl, sd = 2) )
        f.ind= match(flightday, clim$J)
        #if NA day, use mean
        f.ind[is.na(f.ind)]<-match(Jfl, clim$J)
        
        #calculate geometric mean of egg viability within flight period
        ev.ind=sapply(f.ind, function(x)  geo_mean(EggViab[(x-2):(x+2)]) )
        #AVERAGE FAT OVER DAYS
        FAT.ind= sapply(f.ind, function(x)  mean(FAT[(x-2):(x+2)], na.rm=TRUE) )
        #AVERAGE TEMP
        T.ind= sapply(f.ind, function(x)  mean(Temps[(x-2):(x+2)], na.rm=TRUE) )
        
        Eggs.ind= 60*PropFlight*OviRate*FAT.ind * ev.ind #account for Egg viability
        Eggs.ind_noViab= 60*PropFlight*OviRate*FAT.ind
        
        #Means across individuals
        Eggs= mean(Eggs.ind)
        Eggs_noViab= mean(Eggs.ind_noViab)
        EV1[1,spec.k]= mean(FAT.ind)
        EV1[2,spec.k]= mean(ev.ind)
        
        if(!is.nan(Eggs)){
          MaxDay=5
          Lambda1=0
          for(day in 1:MaxDay){
            Eggs1= min(Eggs, MaxEggs-Eggs*(day-1))  ###LIMIT MAX NUMBER EGGS
            if(Eggs1<0) Eggs1=0
            Lambda1= Lambda1+ SurvMat * SurvDaily^day *Eggs1;                        
          }#end loop days
          
          Lambda[yr.k, cell.k, abs.k, gen.k, ]= c(Lambda1, mean(FAT.ind), mean(ev.ind), mean(T.ind, na.rm=T) ) 
          
        } #Check Eggs
        
      } #end loop absortivity
    } #end loop cells    
  } #end loop generation
  
} #end loop years

#SAVE OBJECT
setwd("/Volumes/GoogleDrive/My\ Drive/Buckley/Work/ColiasBiogeog/OUT/")

filename= paste("lambda1_",projs[proj.k],".rds",sep="")
saveRDS(Lambda, filename)
#Lambda1 <- readRDS("mymodel.rds")

#Write out pupal temps
saveRDS(pup.temps, paste("PupTemps_",projs[proj.k],".rds",sep="") )

#write out points
write.csv(pts.sel, paste("COpoints_",projs[proj.k],".rds",sep="") )

#========================================

#EVOLUTIONARY MODEL
N.ind=1000
a= seq(0.4,0.7,0.05)
#for finding a with max fitness
a.fit= as.data.frame(seq(0.4,0.7,0.01))
names(a.fit)="a"

#Read points
setwd("/Volumes/GoogleDrive/My\ Drive/Buckley/Work/ColiasBiogeog/OUT/") #mac version

pts.sel= read.csv( paste("COpoints.csv", sep="") ) #_",projs[proj.k],"

#Read lambdas and pupal temps
#Lambda[years, sites, abs, gen, metrics: Lambda, FAT,Egg Viability]
Lambda <- readRDS( paste("lambda1_",projs[proj.k],".rds", sep="") )
pup.temps <- readRDS( paste("PupTemps_",projs[proj.k],".rds", sep="") )

#Find years with calculations
counts= rowSums(is.na(pup.temps[6,,,1]))

inds=1:150
years= years[inds]

#==============================================
#Calculate optimal absorptivity

abs.opt= array(NA, dim=c(length(years),nrow(pts.sel), 3))  

for(yr.k in 1:length(years)) {
  
  ##loop through generations in each year
  for(gen.k in 1:ngens) {
    
    Lambda.yr.gen= Lambda[yr.k, , , gen.k, ]
    
    #Extract temperatures
    Tp= pup.temps["Tpup",yr.k, , gen.k]
    
    #--------------------------
    #Fitness models
    
    if(!all(is.na(Lambda.yr.gen[,,1]))){ #check has data
      
      #Estimate fitness functions across cells
      fit= array(unlist(apply(Lambda.yr.gen[,,1], 1, function(x) if(sum(is.na(x))==0) lm(x~a+I(a^2))$coefficients)), dim=c(3, nrow(pts.sel)) )
      #Save model
      fit.mod= apply(Lambda.yr.gen[,,1], 1, function(x) if(sum(is.na(x))==0) lm(x~a+I(a^2)) )
      
      #find maxima lambda
      abs.opt[yr.k,,gen.k]= as.vector(array(unlist(sapply(fit.mod, function(x) if(!is.null(x))a.fit$a[which.max(predict.lm(x, a.fit))] )), dim=c(1, nrow(pts.sel)) ) )
      
    } #end check data
  } #end gen loop
} #end year loop

#save optimal Alphas
setwd("/Volumes/GoogleDrive/My\ Drive/Buckley/Work/ColiasBiogeog/OUT/")

#saveRDS(abs.opt, paste("abs.opt_",projs[proj.k],".rds", sep=""))
abs.opt <- readRDS( paste("abs.opt_",projs[proj.k],".rds", sep="") )

#***************************************
#compute initial AbsMean 
int_elev = 0.4226; slope_elev = 0.06517
Tmid = 20; slope_plast = -0.0083  #if Tmid=22.5, -0.006667;

elev_km= pts.sel$elev/1000
abs.init <- int_elev+ slope_elev*elev_km

## NEED TO CALC ABS.OPT
#initialize with optimum value yrs 1950-1960, across generations
abs.init2 <- rowMeans(colMeans(abs.opt[1:10,, ], na.rm=TRUE))

plot(elev_km, abs.init, ylim=range(0.5, 0.7), type="l")
points(elev_km, abs.init2)

#Use optimal
abs.init<- abs.init2

#-----------------------
#Save values
abs.mean= array(NA, dim=c(length(years),nrow(pts.sel), 3, 5,5))  #dims: yr.k, cell.k, gen.k, scen.k:no plast, plast, only plast, metrics: abssample, absmid, rn, Babsmid, Brn)
abs.mean[1,,1,,2]= abs.init
abs.mean[1,,1,,3]= slope_plast
dimnames(abs.mean)[[5]]= c("abssample", "absmid", "rn", "Babsmid", "Brn") 

lambda.mean= array(NA, dim=c(length(years),nrow(pts.sel), 3, 5)) #dims: yr.k, cell.k, gen.k, scen.k:no plast, plast, only plast)

BetaRN= rep(NA, nrow(pts.sel))
#-------------------------------
scen.mat= rbind(c(0,0,0),c(1,0,0),c(0,1,0),c(1,1,0),c(1,1,1) )
colnames(scen.mat)= c("plast","evol","evolRN"  )

for(yr.k in 1:length(years)) {
  
  ##loop through generations in each year
  for(gen.k in 1:ngens) {
    
    BetaAbsmid=NA
    
    Lambda.yr.gen= Lambda[yr.k, , , gen.k, ]
    
    #determine those completing generations
    comp.gen= which(pup.temps["Jadult",yr.k,,gen.k]<243)
    nocomp.gen= which(pup.temps["Jadult",yr.k,,gen.k]==243)
    #set those not completing generations to NA
    if(length(nocomp.gen)>0) Lambda.yr.gen[nocomp.gen,,]=NA
    
    #account for NA lambdas
    l.no.na= which(!is.na(Lambda.yr.gen[,1,1]))
    
    if(length(l.no.na)>0){ #CHECK LAMBDA DATA EXISTS
      
      #Extract temperatures
      Tp= pup.temps["Tpup",yr.k,l.no.na, gen.k]
      
      #--------------------------
      #Fitness models
      #Estimate fitness functions across cells
      fit= array(unlist(apply(Lambda.yr.gen[l.no.na,,1], 1, function(x) if(sum(is.na(x))==0) lm(x~a+I(a^2))$coefficients)), dim=c(3, nrow(pts.sel)) )
      #Save model
      fit.mod= apply(Lambda.yr.gen[l.no.na,,1], 1, function(x) if(sum(is.na(x))==0) lm(x~a+I(a^2)) )
      
      #find maxima lambda
      abs.max= as.vector(array(unlist(sapply(fit.mod, function(x) if(!is.null(x))a.fit$a[which.max(predict.lm(x, a.fit))] )), dim=c(1, length(l.no.na)) ) )
      
      #-------------------------
      # LOOP PLASTICITY SCENARIOS
      for(scen.k in 1:5){ #plast0evol0, plast1evol0, plast0evol1, plast1evol1, plast1evol1rnevol1
        
        if(scen.mat[scen.k,1]==1) rn.mean1= rep(slope_plast, length(l.no.na) )
        if(scen.mat[scen.k,1]==0) rn.mean1= rep(0, length(l.no.na) )
        if(scen.k==5 & gen.k==1) rn.mean1= abs.mean[yr.k,l.no.na,gen.k,scen.k,"rn"]
        if(scen.k==5 & gen.k>1) rn.mean1= abs.mean[yr.k,l.no.na,gen.k-1,scen.k,"rn"]
        
        if(gen.k==1) abs.mean1= abs.mean[yr.k,l.no.na,gen.k,scen.k,"absmid"]
        if(gen.k>1) abs.mean1= abs.mean[yr.k,l.no.na,gen.k-1,scen.k,"absmid"]
        
        #change NA values to negative values 
        abs.na.inds= abs.mean1[which( is.na(abs.mean1))]
        rn.na.inds= rn.mean1[which( is.na(rn.mean1))]
        
        #   #check abs mean
        #    if(!all(is.na(abs.mean1))){
        
        abs.mean1[which( is.na(abs.mean1))]= -10 
        rn.mean1[which( is.na(rn.mean1))]= -1000
        
        #Choose random sample of abs and rn values from current distribution (truncated normal) 
        abs.sample= sapply(abs.mean1, function(x) rtnorm(N.ind, mean = x, sd = abs.sd, lower=0.400, upper=0.700) )
        rn.sample= sapply(rn.mean1, function(x) rtnorm(N.ind, mean = x, sd = rn.sd, lower=-1, upper=1) )
        if(scen.mat[scen.k,1]==0) rn.sample[]=0
        
        #Add plasticity across sites and sample
        abs.plast <- abs.sample + rn.sample*(Tp-Tmid)
        #abs.mean[yr.k,,gen.k] <- abs.mean[yr.k,,gen.k]+abs.plast
        
        ##calculate fitness
        #use fitness function to predict Lambda for each individual
        #extract coefficients and calculate across abs samples
        fit.sample= foreach(cell.k=1:length(l.no.na), .combine="cbind") %do% {
          sapply(abs.plast[,cell.k], function(x) if( sum(is.na(fit[,cell.k]))==0) fit[1,cell.k]+x*fit[2,cell.k]+x^2*fit[3,cell.k] )
        } 
        #Fit.pred <- eval.fd(Abs.sample,Fitmod.year.gen) ### for spline
        
        #standardize to relative fitness and centered on trait mean
        fit.mean= colMeans(fit.sample)
        lambda.mean[yr.k,l.no.na,gen.k,scen.k]=fit.mean
        rel.fit= fit.sample/fit.mean
        
        absmid.dif= t( apply(abs.sample,1,'-',abs.mean1) )
        rn.dif= t( apply(rn.sample,1,'-',rn.mean1) )
        
        R2selnAbsmid<- rep(0, length(l.no.na) ) #No response to selection if no evolution
        R2selnRN<- rep(0, length(l.no.na) ) 
        #------------
        if(scen.k<5 & scen.mat[scen.k,2]==1){    
          ##selection analysis
          sel.fit= sapply(1:length(l.no.na), function(x) if(sum(is.na(x))==0) lm(rel.fit[,x]~absmid.dif[,x] +I(absmid.dif[,x]^2))$coefficients)
          
          #Save model
          sel.mod= sapply(1:length(l.no.na), function(x) if(sum(is.na(x))==0) lm(rel.fit[,x]~absmid.dif[,x] +I(absmid.dif[,x]^2) ) )
          ## EXTRACT SUMMARY?:   fitr2 <- summary(lm.fitmod.yr)$r.squared
          
          #Response to selection
          BetaAbsmid <-sel.fit[2,]
          R2selnAbsmid <- h2*(abs.sd^2)*BetaAbsmid
        } #end scen.k<5
        #------------
        if(scen.k==5){    
          ##selection analysis
          sel.fit= sapply(1:length(l.no.na), function(x) if(sum(is.na(x))==0) lm(rel.fit[,x]~absmid.dif[,x] + rn.dif[,x] +I(absmid.dif[,x]^2) +I(rn.dif[,x]^2)+ rn.dif[,x]*absmid.dif[,x])$coefficients)
          
          #Save model
          sel.mod= sapply(1:length(l.no.na), function(x) if(sum(is.na(x))==0) lm(rel.fit[,x]~absmid.dif[,x] + rn.dif[,x] +I(absmid.dif[,x]^2) +I(rn.dif[,x]^2)+ rn.dif[,x]*absmid.dif[,x]) )
          ## EXTRACT SUMMARY?:   fitr2 <- summary(lm.fitmod.yr)$r.squared
          
          #Response to selection
          BetaAbsmid <-sel.fit[2,]
          R2selnAbsmid <- h2*(abs.sd^2)*BetaAbsmid
          
          BetaRN <- sel.fit[3,] 
          R2selnRN <- h2*(rn.sd^2)*BetaRN
        } #end scen.k==5
        #-------------
        
        #Response to selection
        if(gen.k<3) {
          abs.mean[yr.k,l.no.na,gen.k+1,scen.k,"absmid"]= abs.mean[yr.k,l.no.na,gen.k,scen.k,"absmid"] + R2selnAbsmid
          #Constain abs
          abs.mean[yr.k,which(abs.mean[yr.k,l.no.na,gen.k+1,scen.k,"absmid"]>0.7),gen.k+1,scen.k,"absmid"]=0.7
          abs.mean[yr.k,which(abs.mean[yr.k,l.no.na,gen.k+1,scen.k,"absmid"]<0.4),gen.k+1,scen.k,"absmid"]=0.4   
          
          #rn evolution
          if(scen.k==5){
            abs.mean[yr.k,l.no.na,gen.k+1,scen.k,"rn"]= abs.mean[yr.k,l.no.na,gen.k,scen.k,"rn"] + R2selnRN
            #Constain abs
            abs.mean[yr.k,which(abs.mean[yr.k,l.no.na,gen.k+1,scen.k,"rn"]>1),gen.k+1,scen.k,"rn"]= 1
            abs.mean[yr.k,which(abs.mean[yr.k,l.no.na,gen.k+1,scen.k,"rn"]< -1),gen.k+1,scen.k,"rn"]= -1
          }
        } #end evolutionary scenarios
        
        #Account for missing lambdas
        if(length(abs.na.inds)>0)  R2selnAbsmid[abs.na.inds]=NA    
        if(length(rn.na.inds)>0)   R2selnRN[rn.na.inds]=NA  
        
        #also put in next year's slot
        abs.mean[yr.k+1,l.no.na,1,scen.k,"absmid"]= abs.mean[yr.k,l.no.na,gen.k,scen.k,"absmid"] + R2selnAbsmid
        #Constain abs
        abs.mean[yr.k+1,which(abs.mean[yr.k+1,l.no.na,1,scen.k,"absmid"]>0.7),1,scen.k,"absmid"]=0.7
        abs.mean[yr.k+1,which(abs.mean[yr.k+1,l.no.na,1,scen.k,"absmid"]<0.4),1,scen.k,"absmid"]=0.4 
        
        if(scen.k==5) abs.mean[yr.k+1,l.no.na,1,scen.k,"rn"]= abs.mean[yr.k,l.no.na,gen.k,scen.k,"rn"] + R2selnRN
        #Constain abs
        abs.mean[yr.k+1,which(abs.mean[yr.k+1,l.no.na,gen.k,scen.k,"rn"]>1),1,scen.k,"rn"]= 1
        abs.mean[yr.k+1,which(abs.mean[yr.k+1,l.no.na,gen.k,scen.k,"rn"]< -1),1,scen.k,"rn"]= -1
        
        #Store other metrics
        abs.mean[yr.k,l.no.na,gen.k,scen.k,"abssample"]= colMeans(abs.plast)
        abs.mean[yr.k,l.no.na,gen.k,scen.k,"Babsmid"]= BetaAbsmid
        if(scen.k==5) abs.mean[yr.k,l.no.na,gen.k,scen.k,"Brn"]= BetaRN
        
      } #end scen loop
      
    } #Check lambda values exist
    
  } #end generation
  print(yr.k)
} #end year 

#=====================================
#Save output

setwd("/Volumes/GoogleDrive/My\ Drive/Buckley/Work/ColiasBiogeog/OUT/") #mac version

#saveRDS(abs.mean, "absmean.abs")
#saveRDS(lambda.mean, "lambdamean.abs")

abs.mean <- readRDS("absmean.abs")
lambda.mean <- readRDS("lambdamean.abs")

#================================