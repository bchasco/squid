cog_plot <- function(fit, 
                     plot=TRUE, 
                     col_arg=1,
                     myYears=c(2014:2016,2019),
                     elipse_alpha = 0.05,
                     text_alpha = 1){
  library(shape)
  #Name for derived variable
  CogName <- "mean_Z_ctm"
  Sdreport <- fit$parameter_estimates$SD
  SD <- TMB::summary.sdreport(Sdreport)
  COG <- SD[which(rownames(SD)==CogName),]
  cog_est <- array(COG[,c('Estimate')],c(nrow(COG)/2,2))
  cog_sd <- array(COG[,c('Std. Error')],c(nrow(COG)/2,2))
  
  source("func_ellipse.r")
  
  if(plot==TRUE){
    plot(cog_est,
         # xlim=c(min(cog_est[,1]-cog_sd[,1]),max(cog_est[,1]+cog_sd[,1])),
         # ylim=c(min(cog_est[,2]-cog_sd[,2]),max(cog_est[,2]+cog_sd[,2])),
         las=1,
         ylab="",
         xlab="Easting (km)",
         xlim=c(375,520),
         ylim=c(4000,4600),
         xaxs="i",
         yaxs="i",
         type="n")
    mtext("Northing (km)",side=2,line=4)
  }
  
  
  for(i in 2:nrow(cog_est)){
    Arrows(cog_est[i-1,1],
           cog_est[i-1,2],
           cog_est[i,1],
           cog_est[i,2],
           arr.type="triangle",
           arr.width = 0.2,
           arr.length = 0.2,
           lwd=0.5,
           col="lightgrey")
  }
  
  years <- sort(unique(raw$Year))
  for(i in 1:nrow(cog_est)){
  #   # i <- 2
    cog_ellipse <- ellipse_func(cog_est[i,1],
                                cog_est[i,2],
                                cog_sd[i,1],
                                cog_sd[i,2],
                                step=10)

    if(col_arg==1){
      mycol <- rgb(0.1,0.9,0.2,elipse_alpha)
      textcol <- rgb(0.1,0.9,0.2,text_alpha)
    }
    if(col_arg==2){
      mycol <- rgb(0.2,0.1,0.9,elipse_alpha)
      textcol <- rgb(0.2,0.1,0.9,text_alpha)
    }
    if(col_arg==3){
      mycol <- rgb(0.9,0.1,0.2,elipse_alpha)
      textcol <- rgb(0.9,0.1,0.2,text_alpha)
    }
    polygon(cog_ellipse$x,
            cog_ellipse$y,
            col=mycol,
            border = FALSE)
    if(years[i]%in%myYears){
      text(cog_est[i,1],
           cog_est[i,2],
           years[i],
           cex=1.2,
           pos = 4,
           col = "black")
           # col = textcol)
    }
  }
  # 
  if(plot==TRUE){
    library(mapproj)
    library(sp)
    library(rgdal)
    xy <- data.frame(ID = 1, 
                     X = c(-122.4194,-123.8191551,-124.1839227,-124.4084, -124.09344), 
                     Y = c(37.7749,39.4447541,40.7850507,43.1190, 46.246922))
    coordinates(xy) <- c("X", "Y")
    proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
    res <- spTransform(xy, CRS("+proj=utm +zone=10 ellps=WGS84"))
    # points(res@coords/1000, cex=10)
    # text(510, res@coords[,c('Y')]/1000,
    #      labels=c("San \n Francisco", "Fort \n Bragg", 'Eureka', 'Bandon','Columbia \n River'),
    #      pos=1,
    #      cex=0.8,
    #      srt=90
    # )
    
    op <- par()
    # par(mar=c(6,7,4,2)+0.1)
    axis(4,
         at=res@coords[,c('Y')]/1000,
         labels=c("San \n Francisco", "Fort \n Bragg", 'Eureka', 'Bandon','Columbia \n River'),
    )
    # axis(4, res@coords[,c('Y')]/1000,
    #      labels=rep("",5),
    #      line=3,
    #      # labels=c("San /n Francisco", "Fort Bragg", 'Eureka'),
    #      tck = 0.025,
    #      las=0,
    #      srt=180
    # )
    par(op)
  }
}
# 
# # https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html
# # https://ggplot2.tidyverse.org/reference/coord_map.html
# # https://rdrr.io/github/James-Thorson-NOAA/FishStatsUtils/src/R/make_extrapolation_info.R
# # https://www.earthdatascience.org/courses/earth-analytics/spatial-data-r/geographic-vs-projected-coordinate-reference-systems-UTM/
# 
# 
# #Jim has someting called convert_shapefile to change coordinate system from LL to UTM northing and easting.
# 
# library("rnaturalearth")
# library("rnaturalearthdata")
# library("sf")
# library(rgdal)
# library(ggplot2)
# #Change northing and easting back into lat and long
# utm <- data.frame(Easting=cog_est[,1], Northing=cog_est[,2])
# utm <- utm[complete.cases(utm),]
# utm1 <- data.frame(y=utm$Northing,x=utm$Easting) 
# coordinates(utm1) <- ~x+y
# # class(utm1)
# proj4string(utm1) <- CRS("+proj=utm +zone=10 +datum=WGS84 +units=km")
# utm2 <- spTransform(utm1,CRS("+proj=longlat +datum=WGS84"))
# 
# 
# usamap <- ggplot(states, aes(long, lat)) +
#   geom_polygon(fill = "white", colour = "white") +
#   ylim(34,49) + 
#   # coord_sf(xlim = c(-127, -115.00), ylim = c(34.0, 55.), expand = FALSE) +
#   geom_point(data=as.data.frame(utm2@coords),aes(x=x,y=y, group=NA))
# 
# # Use cartesian coordinates
# usamap
# 
# print(p)
