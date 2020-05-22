library(INLA); library(sp); library(fields)
library(geoR)
library(viridisLite)
library(TMB)
library(ggthemes)
library(maps)
require(sp)
require(maptools)
require(ggplot2)

survey <- read.csv("comb_catches.txt", sep="\t", header = TRUE)
stations <- unique(survey$Station_Code)
survey$Year <- survey$Year - min(survey$Year)


n_i <- nrow(survey)
c_i <- survey$CPUE
s_i <- data.frame(Station=survey$Station)
si_lu <- data.frame(s_i=1:length(unique(survey$Station)),
                    Station=unique(survey$Station))
s_i <- merge(s_i,si_lu)$s_i-1

lat_i <- data.frame(lat=round(survey$Lat,1))
lati_lu <- data.frame(lat_i=1:length(unique(lat_i)[,1]),
                    lat=sort(unique(lat_i)[,1]))
lat_i <- merge(lat_i,lati_lu)$lat_i-1
lat_dist <- lati_lu$lat

t_i <- survey[,'Year']-min(survey[,'Year']) 
n_x <- nrow(survey)
myVars <- c("X3m_Temp")#,"X3m_Salinity","X3m_Chl","Distance_Offshore_nmi","Station_Depth_m")
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
use_p_lat <- 0
use_c_lat <- 0
use_p_yl <- 0
use_c_yl <- 1

try(dyn.unload("squid_v1_delta_glm_AR1yr_MV1Dspace"))
# compile("squid_v1_delta_glm_AR1yr_MV1Dspace.cpp")
dyn.load("squid_v1_delta_glm_AR1yr_MV1Dspace")

Data = list( scaler = 0.1,
             n_i=n_i, 
             c_i=c_i, 
             s_i=s_i, 
             t_i=t_i, 
             lat_i = lat_i,
             lat_dist = lat_dist,
             X_xp=as.matrix(X_xp),
             use_p_y = use_p_y,
             use_c_y = use_c_y,
             use_p_s = use_p_s,
             use_c_s = use_c_s,
             use_p_lat = use_p_lat,
             use_c_lat = use_c_lat,
             use_p_yl = use_p_yl,
             use_c_yl = use_c_yl
)

Parameters = list(beta_c = rep(0,ncol(Data$X_xp)),
                  beta_p = rep(0,ncol(Data$X_xp)),
                  eps_c_s = unique(Data$s_i)*0,
                  eps_p_s = unique(Data$s_i)*0,
                  eps_c_y = unique(Data$t_i)*0,
                  eps_p_y = unique(Data$t_i)*0,
                  eps_p_lat = Data$lat_dist*0,
                  eps_c_lat = Data$lat_dist*0,
                  eps_p_yl = array(0,c(length(Data$lat_dist),length(unique(Data$t_i)))),
                  eps_c_yl = array(0,c(length(Data$lat_dist),length(unique(Data$t_i)))),
                  fc_y_sd = 0,
                  fp_y_sd = 0,
                  fc_y_rho = 0,
                  fp_y_rho = 0,
                  fc_yl_rho_y = 0,
                  fp_yl_rho_y = 0,
                  fc_yl_rho_lat = 0,
                  fp_yl_rho_lat = 0,
                  fc_yl_sd = 0,
                  fp_yl_sd = 0,
                  fc_lat_sd = 0,
                  fp_lat_sd = 0,
                  fc_lat_rho = 0,
                  fp_lat_rho = 0,
                  fs_c_sd = 0,
                  fs_p_sd = 0,
                  fsd = 0)

Random = c("eps_c_s","eps_p_s", "eps_c_y", "eps_p_y", "eps_p_lat", "eps_c_lat", "eps_p_yl", "eps_c_yl")

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

if(!use_c_yl){
  myMap <- append(myMap,
                  list(eps_c_yl=as.factor(array(NA,c(nrow(Parameters$eps_c_yl),ncol(Parameters$eps_c_yl))))
                       ,fc_yl_rho_y = as.factor(NA)
                       ,fc_yl_rho_lat = as.factor(NA)
                       ,fc_yl_sd = as.factor(NA)
                  )
  )
}

if(!use_p_yl){
  myMap <- append(myMap,
                  list(eps_p_yl=as.factor(array(NA,c(nrow(Parameters$eps_p_yl),ncol(Parameters$eps_p_yl))))
                       ,fp_yl_rho_y = as.factor(NA)
                       ,fp_yl_rho_lat = as.factor(NA)
                       ,fp_yl_sd = as.factor(NA)
                  )
  )
}

# Make object
Obj = MakeADFun(data=Data,
                parameters=Parameters,
                random=Random,
                map = myMap,
                hessian=FALSE,
                DLL="squid_v1_delta_glm_AR1yr_MV1Dspace")


out <- nlminb(Obj$par,Obj$fn,Obj$gr)
#This is the output of the model
rep <- Obj$report()

par(mfrow=c(4,2))
image(t(rep$eps_p_yl))
image(t(rep$eps_c_yl))
plot(rep$eps_p_lat, type="l")
plot(rep$eps_c_lat, type="l")
plot(rep$eps_p_y, type="l")
plot(rep$eps_c_y, type="l")
hist(rep$eps_p_s, breaks=20)
hist(rep$eps_c_s, breaks=20)

rm(SD)
SD <- sdreport(Obj)
print(SD)

library(ggplot2)
library(ggmap)
ggimage((rep$eps_c_yl), fullpage = TRUE)
