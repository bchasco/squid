
library(maps)
library(mapdata)
library(mapproj)

# col_seq <- seq(floor(min(rep$pred_c_ef)),ceiling(max(rep$pred_c_ef)),abs(floor(min(rep$pred_c_ef))-ceiling(max(rep$pred_c_ef)))/50)
col_seq <- exp(seq(floor(min(rep$beta_c[1]+rep$pred_c_ef)),ceiling(max(rep$beta_c[1]+rep$pred_c_ef)),
               abs(ceiling(max(rep$beta_c[1]+rep$pred_c_ef))-floor(min(rep$beta_c[1]+rep$pred_c_ef)))/50))
par(mfrow=c(5,5))

obs_mat <- matrix(NA,nrow(loi_lu),nrow(lati_lu))
for(loi in c(loi_lu$lo_i)){
  for(lati in c(lati_lu$lat_i)){
    xlo <- loi_lu$lo[Data$lo_i]==loi_lu$lo[loi]
    xla <- lati_lu$lat[Data$lat_i]==lati_lu$lat[lati]
    if(max(xla + xlo)==2){
      obs_mat[loi,lati] <- 1
    }
  }
}
    
for(i in 1:22){
  par(mai=c(0.4,0.6,0.0,0.0))
  x <- exp(rep$beta_c[1]+rep$pred_c_ef[,,i])
  # map("usa", xlim=c(-126,-122), ylim=c(36,50), col=grey(0.9), fill=TRUE, yaxs="i", xaxs="i", add=FALSE)  
  # x[is.na(obs_mat)] <- NA
  image(loi_lu$lo
        , lati_lu$lat
        , x
        # ,main=i+1997
        ,col=rainbow(length(col_seq)-1)
        ,breaks=col_seq
        ,las=1
        ,xlim=c(-125.5,-122)
        ,ylab=""
        # , ylim=c(36,49.5) 
        # ,add=TRUE
  )
  map("usa"
      , xlim=c(-125.5,-122)
      # , ylim=c(36,49.5) 
      , col=grey(0.9)
      , fill=TRUE
      , yaxs="i"
      , xaxs="i"
      , add=TRUE)  
  
  points(survey$Long[survey$Year==(i-1)]
         ,survey$Lat[survey$Year==(i-1)]
         , pch=16
         , col=rgb(0,0,0,0.3)
         , cex=0.75)
  
  text(x=-122.5,y=46,i+1997, cex=1.5)
  box()
}

