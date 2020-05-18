library(INLA); library(sp); library(fields)
library(geoR)
library(viridisLite)
library(TMB)
# library(TMBhelper)
library(ggthemes)
library(maps)
require(sp)
require(maptools)
require(ggplot2)

survey <- read.csv("comb_catches.txt", sep="\t", header = TRUE)
stations <- unique(survey$Station_Code)
survey$Year <- survey$Year - min(survey$Year)

round.to <- function(x, b) {
  round(x/b)*b
}
lat_bin <- 0.5
lo_bin <- 0.2

n_i <- nrow(survey)
c_i <- survey$CPUE
s_i <- data.frame(Station=survey$Station)
si_lu <- data.frame(s_i=1:length(unique(survey$Station)),
                    Station=unique(survey$Station))
s_i <- merge(s_i,si_lu)$s_i-1

lat <- data.frame(lat=as.numeric(round.to(survey$Lat,lat_bin)))

lat_range <- round.to(seq(min(as.integer(floor(lat$lat))),
                       max(as.integer(ceiling(lat$lat))),lat_bin),
                   lat_bin)
lati_lu <- data.frame(lat_i=1:length(lat_range),
                      lat=lat_range)

lat_i_tmp <- NA
for(i in 1:dim(lat)[1]){
  lat_i_tmp[i] <- lati_lu$lat_i[lat$lat[i]==lati_lu$lat]
}

lo <- data.frame(lo=as.numeric(round.to(survey$Long,lo_bin)))
lo_range <- seq(min(as.integer(floor(lo$lo))),
                max(as.integer(ceiling(lo$lo))),lo_bin)
loi_lu <- data.frame(lo_i=1:length(lo_range),
                      lo=lo_range)
lo_i_tmp <- NA
for(i in 1:dim(lo)[1]){
  lo_i_tmp[i] <- loi_lu$lo_i[round(lo$lo[i]*100,0)==round(loi_lu$lo*100,0)]
}
length(lo_i_tmp)

lat_dist <- lati_lu$lat

t_i <- survey[,'Year']-min(survey[,'Year']) 
n_x <- nrow(survey)
myVars <- c("X3m_Temp","X3m_Salinity","X3m_Chl","Station_Depth_m")
X_xp <- cbind(as.matrix(rep(1,n_x),1),survey[,myVars])
# X_xp <- cbind(survey[,c("Year","X3m_Temp","X3m_Salinity","X3m_Chl","survey$Distance_Offshore_nmi",survey$Dep)])
for(i in 2:ncol(X_xp)){
  X_xp[is.na(X_xp[,i]),i]<-mean(na.omit(X_xp[,i]))
  X_xp[,i] <- X_xp[,i] - mean(X_xp[,i])
} 

#Switches 
use_p_y <- 1
use_c_y <- 1
use_p_s <- 1
use_c_s <- 1
use_p_lat <- 1
use_c_lat <- 1
use_p_lo <- 1
use_c_lo <- 1
use_p_yll <- 1
use_c_yll <- 1

try(dyn.unload("squid_v1_delta_glm_AR1yr_AR2Dspace"))
# compile("squid_v1_delta_glm_AR1yr_AR2Dspace.cpp")
dyn.load("squid_v1_delta_glm_AR1yr_AR2Dspace")

Data = list( scaler = 0.1,
             n_i=n_i, 
             c_i=c_i, 
             s_i=s_i, 
             t_i=t_i, 
             lat_i = lat_i_tmp,
             lo_i = lo_i_tmp,
             lat_dist = lat_dist,
             X_xp=as.matrix(X_xp),
             use_p_y = use_p_y,
             use_c_y = use_c_y,
             use_p_s = use_p_s,
             use_c_s = use_c_s,
             use_p_lat = use_p_lat,
             use_c_lat = use_c_lat,
             use_p_lo = use_p_lo,
             use_c_lo = use_c_lo,
             use_p_yll = use_p_yll,
             use_c_yll = use_c_yll
)

Parameters = list(beta_c = rep(0,ncol(Data$X_xp)),
                  beta_p = rep(0,ncol(Data$X_xp)),
                  eps_c_s = unique(Data$s_i)*0,
                  eps_p_s = unique(Data$s_i)*0,
                  eps_c_y = unique(Data$t_i)*0,
                  eps_p_y = unique(Data$t_i)*0,
                  eps_p_lat = lati_lu$lat_i*0,
                  eps_c_lat = lati_lu$lat_i*0,
                  eps_p_lo = loi_lu$lo_i*0,
                  eps_c_lo = loi_lu$lo_i*0,
                  eps_p_yll = array(0,c(length(loi_lu$lo_i),
                                        length(lati_lu$lat_i),
                                        length(unique(Data$t_i)))),
                  eps_c_yll = array(0,c(length(loi_lu$lo_i),
                                        length(lati_lu$lat_i),
                                        length(unique(Data$t_i)))),
                  fc_y_sd = 0,
                  fp_y_sd = 0,
                  fc_y_rho = 0,
                  fp_y_rho = 0,
                  fc_yll_rho_y = 0,
                  fp_yll_rho_y = 0,
                  fc_yll_rho_lat = 0,
                  fp_yll_rho_lat = 0,
                  fc_yll_rho_lo = 0,
                  fp_yll_rho_lo = 0,
                  fc_yll_sd = 0,
                  fp_yll_sd = 0,
                  fc_lat_sd = 0,
                  fp_lat_sd = 0,
                  fc_lo_sd = 0,
                  fp_lo_sd = 0,
                  fc_lat_rho = 0,
                  fp_lat_rho = 0,
                  fc_lo_rho = 0,
                  fp_lo_rho = 0,
                  fs_c_sd = 0,
                  fs_p_sd = 0,
                  fsd = 0)

Random = c("eps_c_s","eps_p_s", "eps_c_y", "eps_p_y", 
           "eps_p_lat", "eps_c_lat", 
           "eps_p_lo", "eps_c_lo", 
           "eps_p_yll", "eps_c_yll")

myMap <- list()
if(!use_p_s){
  myMap <- append(myMap,
                  list(eps_p_s=as.factor(rep(NA,length(Parameters$eps_p_s)))
                       ,fs_p_sd = as.factor(NA)
                  )
  )
}

if(!use_c_s){
  myMap <- append(myMap,
                  list(eps_c_s=as.factor(rep(NA,length(Parameters$eps_c_s)))
                       ,fs_c_sd = as.factor(NA)
                  )
  )
}

if(!use_p_y){
  myMap <- append(myMap,
                  list(eps_p_y=as.factor(rep(NA,length(Parameters$eps_p_y)))
                       ,fp_y_rho = as.factor(NA)
                       ,fp_y_sd = as.factor(NA)
                  )
  )
}

if(!use_c_y){
  myMap <- append(myMap,
                  list(eps_c_y=as.factor(rep(NA,length(Parameters$eps_c_y)))
                       ,fc_y_rho = as.factor(NA)
                       ,fc_y_sd = as.factor(NA)
                  )
  )
}

if(!use_p_lat){
  myMap <- append(myMap,
                  list(eps_p_lat=as.factor(rep(NA,length(Parameters$eps_p_lat)))
                       ,fp_lat_rho = as.factor(NA)
                       ,fp_lat_sd = as.factor(NA)
                  )
  )
}

if(!use_c_lat){
  myMap <- append(myMap,
                  list(eps_c_lat=as.factor(rep(NA,length(Parameters$eps_c_lat)))
                       ,fc_lat_rho = as.factor(NA)
                       ,fc_lat_sd = as.factor(NA)
                  )
  )
}

if(!use_p_lo){
  myMap <- append(myMap,
                  list(eps_p_lo=as.factor(rep(NA,length(Parameters$eps_p_lo)))
                       ,fp_lo_rho = as.factor(NA)
                       ,fp_lo_sd = as.factor(NA)
                  )
  )
}

if(!use_c_lo){
  myMap <- append(myMap,
                  list(eps_c_lo=as.factor(rep(NA,length(Parameters$eps_c_lo)))
                       ,fc_lo_rho = as.factor(NA)
                       ,fc_lo_sd = as.factor(NA)
                  )
  )
}

if(!use_c_yll){
  myMap <- append(myMap,
                  list(eps_c_yll=as.factor(array(NA,c(dim(Parameters$eps_c_yll))))
                       ,fc_yll_rho_y = as.factor(NA)
                       ,fc_yll_rho_lat = as.factor(NA)
                       ,fc_yll_rho_lo = as.factor(NA)
                       ,fc_yll_sd = as.factor(NA)
                  )
  )
}

if(!use_p_yll){
  myMap <- append(myMap,
                  list(eps_p_yll=as.factor(array(NA,c(dim(Parameters$eps_p_yll))))
                       ,fp_yll_rho_y = as.factor(NA)
                       ,fp_yll_rho_lat = as.factor(NA)
                       ,fp_yll_rho_lo = as.factor(NA)
                       ,fp_yll_sd = as.factor(NA)
                  )
  )
}

# # Make object
Obj = MakeADFun(data=Data,
                parameters=Parameters,
                random=Random,
                map = myMap,
                hessian=FALSE,
                DLL="squid_v1_delta_glm_AR1yr_AR2Dspace")

# 
out <- nlminb(Obj$par,Obj$fn,Obj$gr)
# #This is the output of the model
rep <- Obj$report()

rm(SD)
SD <- sdreport(Obj)
print(SD)


par(mfrow=c(4,2))
plot(lati_lu$lat, rep$eps_p_lat, type="l")
plot(lati_lu$lat, rep$eps_c_lat, type="l")
plot(loi_lu$lo, rep$eps_p_lo, type="l")
plot(loi_lu$lo, rep$eps_c_lo, type="l")
plot(sort(unique(survey$Year))+1998,rep$eps_p_y, type="l")
plot(sort(unique(survey$Year))+1998,rep$eps_c_y, type="l")
hist(rep$eps_p_s, breaks=20)
hist(rep$eps_c_s, breaks=20)
if(use_p_yll){
  par(mfrow=c(5,5))
  for(i in 1:22){
    image((rep$eps_p_yll[,,i]))
  }
}
# if(use_c_yll){
source("fig_ggplot_maps.r")
print(paste("AIC\n",out$objective+2*length(out$par)))
write.csv(file="year effect.csv",cbind(SD$value[names(SD$value)=="eps_c_y"],SD$sd[names(SD$value)=="eps_c_y"]))      
      