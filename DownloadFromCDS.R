library(ncdf4)
library(raster)

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

# set directory to store files
options(epwshiftr.dir = tempdir())
options(epwshiftr.verbose = TRUE)

(nodes <- get_data_node())

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

#https://pjbartlein.github.io/REarthSysSci/netCDF.html

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/AR6projections/")
nc.daily <- nc_open("tasmaxAdjust_day_GFDL-ESM4_ssp585_r1i1p1f1_gr010_TCDF-CDFT23-ERA5Land-1981-2010_19510101-19651231.nc")

lon <- ncvar_get(nc.daily,"lon");# extracts longitude
lat <- ncvar_get(nc.daily,"lat");# extracts latitude
time <- ncvar_get(nc.daily,"time");# extracts time

f <- "tasmaxAdjust_day_GFDL-ESM4_ssp585_r1i1p1f1_gr010_TCDF-CDFT23-ERA5Land-1981-2010_19510101-19651231.nc"
b <- brick(f)  

#Colorado
e <- extent(251, 258, 37, 41)
x <- crop(b, e)

r <- raster(f, band = "tasmaxAdjust", xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat) )
plot(r)


#import the data
tasmax <- ncvar_get(nc.daily,"tasmaxAdjust")


#see https://github.com/lbuckley/HopperBodysize/blob/main/AssembleClimateData.R
#https://github.com/bluegreen-labs/ecmwfr
# non downscaled CMIP6, https://ecmwf-projects.github.io/copernicus-training-c3s/projections-cmip6.html

#Other CDS options
#https://cds.climate.copernicus.eu/cdsapp#!/dataset/sis-agroclimatic-indicators?tab=overview


setwd("/Volumes/GoogleDrive/Shared drives/RoL_FitnessConstraints/projects/BodySize/data/ClimateData/")

# assign your credentials (register here: https://cds.climate.copernicus.eu/user/register)
uid <- "78176"
cds_api_key <- "062c3e77-bcc8-4c56-8e72-4872e7a92be6"

ecmwfr::wf_set_key(user = uid, key = cds_api_key, service = "cds")

k=1

#----
request.month <- list(
  "dataset_short_name" = "reanalysis-era5-single-levels-monthly-means",
  "product_type" = "monthly_averaged_reanalysis",
  "variable" = c("2m_temperature","skin_temperature","snow_depth","snowfall","total_precipitation"),
  "year" = c('1959', '1960', '1961', '1962', '1963', '1964', '1979', '1981', '2006',
             '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2022'),
  "month" = c("03","04","05","06","07","08"),
  "time" = "00:00",
  "area" = "40.06/-105.59/39.99/-105.28",
  "format" = "netcdf",
  "target" = "ERA5_month.nc"
)

request.hr <- list(
  "dataset_short_name" = "reanalysis-era5-single-levels",
  "product_type" = "reanalysis",
  "variable" = c("2m_temperature"), #,"skin_temperature","snow_depth","snowfall","total_precipitation"
  "year" = c('1959', '1960', '1961', '1962', '1963', '1964', '1979', '1981', '2006',
             '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2022'),
  "month" = c("03","04","05","06","07","08"),
  "day" = c('01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30',
            '31'),
  "time" = c('00:00', '01:00', '02:00',
             '03:00', '04:00', '05:00',
             '06:00', '07:00', '08:00',
             '09:00', '10:00', '11:00',
             '12:00', '13:00', '14:00',
             '15:00', '16:00', '17:00',
             '18:00', '19:00', '20:00',
             '21:00', '22:00', '23:00'),
  "area" = "40.06/-105.59/39.99/-105.28",
  "format" = "netcdf",
  "target" = "ERA5_hr.nc"
)

#request
file <- wf_request(
  user     = "78176",   # user ID (for authentification)
  request  = request.hr,  # the request
  transfer = TRUE,     # download the file
  path     = "."       # store data in current working directory
)


#----
#process data

setwd("/Volumes/GoogleDrive/Shared drives/RoL_FitnessConstraints/projects/BodySize/data/ClimateData/")

nc.hr <- nc_open("Evans_ERA5_hr.nc")

#---
#hr data

#extract the time
t <- ncvar_get(nc.hr, "time")
#convert the hours into date + hour
#as_datetime() function of the lubridate package needs seconds
timestamp <- as_datetime(c(t*60*60),origin="1900-01-01")

#import the data
t2m <- ncvar_get(nc.hr,"t2m")
#add soil temp, wind for microclimate?

#combine data
if(k==1) {
  nc.hr.data= cbind(year(timestamp),month(timestamp), day(timestamp), hour(timestamp), yday(timestamp), t2m[1,])
  nc.hr.data= as.data.frame(nc.hr.data)
  names(nc.hr.data)=c("year","month","day","hour","doy","t2m")
  nc.hr.data$site= "Boulder1"
  
  nc.hr.data2= cbind(year(timestamp),month(timestamp), day(timestamp), hour(timestamp), yday(timestamp), t2m[2,])
  nc.hr.data2= as.data.frame(nc.hr.data2)
  names(nc.hr.data2)=c("year","month","day","hour","doy","t2m")
  nc.hr.data2$site= "Boulder2"
  
  nc.hr.all= rbind(nc.hr.data,nc.hr.data2)
}

if(k>1){
  nc.hr.data= cbind(year(timestamp),month(timestamp), day(timestamp), hour(timestamp), yday(timestamp), t2m)
  nc.hr.data= as.data.frame(nc.hr.data)
  names(nc.hr.data)=c("year","month","day","hour","doy","t2m")
  
  #add site info
  if(k==2) nc.hr.data$site= "Rollins"
  if(k==3) nc.hr.data$site= "Evans"
  
  nc.hr.all= rbind(nc.hr.all, nc.hr.data)
}

#---
#close the connection with the ncdf file
nc_close(nc.m)
nc_close(nc.hr)

#write out climate data
saveRDS(nc.hr.all, "nchr.RDS") 

#===================================
#NCCS Thredds
#CATALOG: https://ds.nccs.nasa.gov/thredds/catalog/AMES/NEX/GDDP-CMIP6/catalog.html
#DESCRIPTION: https://www.nccs.nasa.gov/services/data-collections/land-based-products/nex-gddp-cmip6

#R access: https://rdrr.io/github/bocinsky/thredds/f/README.md





