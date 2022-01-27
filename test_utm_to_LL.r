# https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html
# https://ggplot2.tidyverse.org/reference/coord_map.html
# https://rdrr.io/github/James-Thorson-NOAA/FishStatsUtils/src/R/make_extrapolation_info.R
# https://www.earthdatascience.org/courses/earth-analytics/spatial-data-r/geographic-vs-projected-coordinate-reference-systems-UTM/


#Jim has someting called convert_shapefile to change coordinate system from LL to UTM northing and easting.

library("rnaturalearth")
library("rnaturalearthdata")
library("sf")
library(rgdal)
library(ggplot2)
#Change northing and easting back into lat and long
utm <- data.frame(Easting=cog_est[,1], Northing=cog_est[,2])
utm <- utm[complete.cases(utm),]
utm1 <- data.frame(y=utm$Northing,x=utm$Easting) 
coordinates(utm1) <- ~x+y
# class(utm1)
proj4string(utm1) <- CRS("+proj=utm +zone=10 +datum=WGS84 +units=km")
utm2 <- spTransform(utm1,CRS("+proj=longlat +datum=WGS84"))


usamap <- ggplot(states, aes(long, lat)) +
  geom_polygon(fill = "white", colour = "white") +
  ylim(34,49) + 
  # coord_sf(xlim = c(-127, -115.00), ylim = c(34.0, 55.), expand = FALSE) +
  geom_point(data=as.data.frame(utm2@coords),aes(x=x,y=y, group=NA))

# Use cartesian coordinates
usamap

print(p)
