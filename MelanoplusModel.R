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

bs.ms<- c(0.90, 0.48, 8.68 * 10^{-5})
bs.md<- c(0.98, 0.48, 1.24 * 10^{-4})

vCO2= function(M, Tb, elev_m) exp(b1*log(M)+b2*(1/(k*Tb))+b3*elev_m)
lipid.g= function(vCO) (0.7*vCo2/2)


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
