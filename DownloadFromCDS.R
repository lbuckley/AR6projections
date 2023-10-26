library(ncdf4)
library(raster)
library(terra)

#climate factory

#https://esgf-node.ipsl.upmc.fr/esgf-idp/openid/lbuckley

#c3s-cmip5-adjust.bias-adjusted-output.BNU.BNU-ESM.rcp85.mon.atmos.Amon.r1i1p1.gr025.TCDF-CDFT23-ERA5-1981-2010.tasmaxAdjust
#THREDDS Catalog ]   [ WGET Script

#https://cran.r-project.org/web/packages/thredds/index.html

#https://esgf-node.ipsl.upmc.fr/search/cmip6-adjust/

#https://github.com/ideas-lab-nus/epwshiftr/blob/master/R/esgf.R

#Expert users may also want to use the ESGF Search RESTful API from http://esgf-node.ipsl.upmc.fr/esg-search/search?
# https://esgf.github.io/esg-search/ESGF_Search_RESTful_API.html
# https://rdrr.io/github/ideas-lab-nus/epwshiftr/man/EsgfQuery.html

#Data syntax: https://eartharxiv.org/repository/view/2722/
#bias adjusted file
# tasAdjust_day_IPSL-CM6A-LR_ssp585_r1i1p1f1_gr010_TCDF-CDFT23-ERA5Land-1981-2010_20160101-20251231.nc

#https://oceanhealthindex.org/news/cmip_2_what_is_this/
#https://docs.ropensci.org/tidync/
#https://www.climatologylab.org/maca.html

#can obtain and combine climate models using R packages climate futures toolbox (CFT) (https://github.com/earthlab/cft), or climateR (https://mikejohnson51.github.io/climateR-intro/#123)

#---
#CMIP6-Adjust @ IPSL

esgf_query(variable = "rss", experiment = "ssp126", resolution = "100 km", limit = 1)

#R epwshiftr
#https://github.com/ideas-lab-nus/epwshiftr

#install.packages("epwshiftr", repos = "https://hongyuanjia.r-universe.dev")
#install.packages("epwshiftr")
library(epwshiftr)

#Try query
# query1= esgf_query(
#   activity = "CMIP6-Adjust",
#   variable = c("tasAdjust"),
#   frequency = "day",
#   experiment = c("rcp45"),
#   source = c("BCC.bcc-csm1-1"),
#   variant = "r1i1p1",
#   replica = FALSE,
#   latest = TRUE,
#   type = "Dataset",
#   limit = 10000L,
#   data_node = "esgf-node.ipsl.upmc.fr"
# )

#------------
# set directory to store files
options(epwshiftr.dir = tempdir())
options(epwshiftr.verbose = TRUE)

(nodes <- get_data_node())

#https://esgf-node.ipsl.upmc.fr/

# create a CMIP6 output file index
idx <- init_cmip6_index(
  # only consider ScenarioMIP activity
  activity = "ScenarioMIP",
  
  # specify dry-bulb temperature and relative humidity
  variable = c("tas", "hurs"),
  
  # specify report frequent
  frequency = "day",
  
  # specify experiment name
  experiment = c("ssp585"),
  
  # specify GCM name
  source = "AWI-CM-1-1-MR",
  
  # specify variant,
  variant = "r1i1p1f1",
  
  # specify years of interest
  years = c(2050, 2080),
  
  # save to data dictionary
  save = TRUE
)


#----------------
#ISIMIP
#https://www.isimip.org/gettingstarted/input-data-bias-adjustment/
#https://www.isimip.org/gettingstarted/input-data-bias-adjustment/details/100/

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/")
nc.daily <- nc_open("gfdl-esm4_r1i1p1f1_w5e5_hist-nat_tasmax_global_daily_1850_1850.nc")

t <- ncvar_get(nc.daily, "time")
#convert the hours into date + hour
#as_datetime() function of the lubridate package needs seconds
#timestamp <- as_datetime(c(t*60*60),origin="1900-01-01")

lon <- ncvar_get(nc.daily, "lon")
lat <- ncvar_get(nc.daily, "lat")

#import the data
tasmax <- ncvar_get(nc.daily,"tasmax")

#plot
r <- raster(t(tasmax[,,1]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
plot(r)

#----------------------
#Climate factory
#https://esgf-node.ipsl.upmc.fr/projects/cmip6-adjust/
#https://esgf-node.ipsl.upmc.fr/search/cmip6-adjust/
  
#raster strategy
library(raster)
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/")
nc.daily <- nc_open("tasmaxAdjust_day_GFDL-ESM4_ssp585_r1i1p1f1_gr010_TCDF-CDFT23-ERA5Land-1981-2010_19510101-19651231.nc")

lon <- ncvar_get(nc.daily,"lon");# extracts longitude
lat <- ncvar_get(nc.daily,"lat");# extracts latitude
time <- ncvar_get(nc.daily,"time");# extracts time

lon.ind= which(lon>251 & lon<258) 
lat.ind= which(lat>37 & lat<41)

nc.sub <- ncvar_get(nc.daily, "tasmaxAdjust",
                     start = c(lon.ind[1], lat.ind[1], 1),
                     count = c(length(lon.ind),length(lat.ind),length(time)))

nc.rast<- raster(nc.sub[,,1])
plot(nc.rast)

#-------
# f <- "tasmaxAdjust_day_GFDL-ESM4_ssp585_r1i1p1f1_gr010_TCDF-CDFT23-ERA5Land-1981-2010_19510101-19651231.nc"
# b <- brick(f)  
# 
# #Colorado
# e <- extent(251, 258, 37, 41)
# x <- crop(b, e)
# 
# r <- raster(f, band = "tasmaxAdjust", xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat) )
# plot(r)

#-------
#processNC approach
#https://github.com/RS-eco/processNC

# #if(!"remotes" %in% installed.packages()[,"Package"]) install.packages("remotes")
# #if(!"processNC" %in% installed.packages()[,"Package"]) remotes::install_github("RS-eco/processNC", build_vignettes=T)
# library(processNC)
# 
# tas_files <- list.files("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/", 
#                         pattern="tasmax.*\\.nc", full.names=T)
# 
# nc.sub2<- subsetNC(files=tas_files, ext=c(251, 258, 37, 41), 
#          startdate=1960, enddate=1965, varid="tasmaxAdjust")

#-------
# #terra approach
# library(terra)
# 
# setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/")
# nc.daily <- rast("tasmaxAdjust_day_GFDL-ESM4_ssp585_r1i1p1f1_gr010_TCDF-CDFT23-ERA5Land-1981-2010_19510101-19651231.nc")
# 
# # define extent
# e <- c(251, 258, 37, 41) |> ext()
# 
# # crop raster to extent
# nc_crop <- crop(nc.daily, e)

#===================================
#NCCS Thredds
#CATALOG: https://ds.nccs.nasa.gov/thredds/catalog/AMES/NEX/GDDP-CMIP6/catalog.html
#DESCRIPTION: https://www.nccs.nasa.gov/services/data-collections/land-based-products/nex-gddp-cmip6

#R access: https://rdrr.io/github/bocinsky/thredds/f/README.md

#https://pjbartlein.github.io/REarthSysSci/netCDF.html




