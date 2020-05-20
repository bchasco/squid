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

try(dyn.unload("squid_v1_delta_glm_AR1yr_MVspace"))
compile("squid_v1_delta_glm_AR1yr_MVspace.cpp")
dyn.load("squid_v1_delta_glm_AR1yr_MVspace")

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
X_xp <- cbind(as.matrix(rep(1,n_x),1),survey[,c("X3m_Temp","X3m_Salinity","X3m_Chl")])
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

Data = list( n_i=n_i, 
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
             use_c_lat = use_c_lat
)

Parameters = list(beta_c = rep(0,ncol(Data$X_xp)),
                  beta_p = rep(0,ncol(Data$X_xp)),
                  eps_c_s = unique(Data$s_i)*0,
                  eps_p_s = unique(Data$s_i)*0,
                  eps_c_y = unique(Data$t_i)*0,
                  eps_p_y = unique(Data$t_i)*0,
                  eps_p_lat = Data$lat_dist*0,
                  eps_c_lat = Data$lat_dist*0,
                  fc_y_sd = 0,
                  fp_y_sd = 0,
                  fc_y_rho = 0,
                  fp_y_rho = 0,
                  fc_lat_sd = 0,
                  fp_lat_sd = 0,
                  fc_lat_rho = 1,
                  fp_lat_rho = 1,
                  fs_c_sd = 0,
                  fs_p_sd = 0,
                  fsd = 0)

Random = c("eps_c_s","eps_p_s", "eps_c_y", "eps_p_y", "eps_p_lat", "eps_c_lat")

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

# Make object
Obj = MakeADFun(data=Data,
                parameters=Parameters,
                random=Random,
                map = myMap,
                hessian=FALSE,
                DLL="squid_v1_delta_glm_MVspace")


out <- nlminb(Obj$par,Obj$fn,Obj$gr)
SD <- sdreport(Obj)
# 
#This is the output of the model
rep <- Obj$report()

print(SD)
