#https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html
library(ncdf4)
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting

#Grab the sst data
nc_data <- nc_open("sst.mnmean.nc")

#Get the lng data
lon <- ncvar_get(nc_data, "lon")
lon[lon > 180] <- lon[lon > 180] - 360

#Get the lat data
lat <- ncvar_get(nc_data, "lat", verbose = F)
t <- ncvar_get(nc_data, "time")
