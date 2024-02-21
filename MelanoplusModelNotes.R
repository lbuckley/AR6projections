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
# temperature dependence of MR (ml CO2/hr): https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.12083
# converted into lipid consumed based on a 2L O2 consumed for every 1g of lipid consumed
#assume respiratory quotient of 0.7 for lipids update
#converted to energy use assuming 39 kJ g^−1 lipid

#bs.ms<- c(0.90, 0.48, 8.68 * 10^{-5})
#bs.md<- c(0.98, 0.48, 1.24 * 10^{-4})

#vCO2= function(M, Tb, elev_m) exp(b1*log(M)+b2*(1/(k*Tb))+b3*elev_m)
#lipid.g= function(vCO) (0.7*vCo2/2)

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
#Sibly, R.M., Monk, K., 1987. A theory of grasshopper life cycles. Oikos 48, 186–194.
# k is the rate at which resources destined for egg biomass are accu- mulated

#Branson 2003 https://doi.org/10.4039/n02-061
#Cage experiment with Melanoplus sanguinipes in Montana
#Regressions of eggs laid as a function of digestible biomass (but weird units)

#Branson 2008 https://doi.org/10.1665/1082-6467-17.2.259
#reproductive correlates in Melanoplus sanguinipes

#Branson 2004 https://doi.org/10.1665/1082-6467(2004)013[0239:RIONAA]2.0.CO;2
#Reproductive allocation in experiments

#Branson 2022 https://doi.org/10.3390/geographies2010003, Wyoming abundance sampling

#I. Filin and O. Ovadia . 2007 
# *** Energetic population model for M. femurrubrum https://doi.org/10.1086/522091 

#Hatle et al 2003 https://doi.org/10.1093/icb/43.5.635
# Reproductive allocation through time

#Jonas et al 2015 https://doi.org/10.1016/j.rama.2014.12.011
#population dynamic regressions

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

#Ackman and Whitman 2008. https://doi.org/10.1665/1082-6467-17.2.249
# size and fecundity
# Romalea microptera
# Fecundity (= reproductive) variables included clutch size Pod 1, clutch size Pod 2, total eggs Pods 1 2, mass Pod 1, mass Pod 2, time between adult eclosion and Pod 1, time between Pod 1 and Pod 2, and time between eclosion and Pod 2. 

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

#ENERGETICS
#Jackson et al 1968, https://doi.org/10.1016/0022-1910(68)90208-4
#Energetics of M. bivittatus eggs

#Kohler et al 1987, https://www.jstor.org/stable/pdf/4218439.pdf
#Energy partitioning across grasshoppers species
#proportion of production used for egg production: 28-44%
#In C. parallelus, with an average fecundity of 17 eggs/ female, this loss is 391 J per average individual, i.e. without distinguishinh between males and females
#G. rufus: 643J (3.1% of energy ingested)

#De Sousa Santos and Begon 1987 https://www.jstor.org/stable/2389423 
#Estimates of reproductive energetics for UK grasshopper 
#Energy values of eggs etc

#Hook 1971, https://esajournals.onlinelibrary.wiley.com/doi/abs/10.2307/1942433 
#Energy and nutrient dynamics of spider and orthopteran populations in a grassland ecosystem
#Energy budget but not for reproduction
#------------------------
