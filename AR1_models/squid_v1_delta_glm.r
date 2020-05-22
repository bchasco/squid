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

try(dyn.unload("squid_v1_delta_glm"))
compile("squid_v1_delta_glm.cpp")
dyn.load("squid_v1_delta_glm")

n_i <- nrow(survey)
c_i <- survey$CPUE
s_i <- data.frame(Station=survey$Station)
si_lu <- data.frame(s_i=1:length(unique(survey$Station)),
                    Station=unique(survey$Station))
s_i <- merge(s_i,si_lu)$s_i-1

t_i <- survey[,'Year']-min(survey[,'Year']) 
n_x <- nrow(survey)
# X_xp <- cbind(as.matrix(rep(1,n_x),1),survey[,c("Year","X3m_Temp","X3m_Salinity","X3m_Chl")])
X_xp <- cbind(survey[,c("Year","X3m_Temp","X3m_Salinity","X3m_Chl")])
for(i in 1:ncol(X_xp)) 
  X_xp[is.na(X_xp[,i]),i]<-mean(na.omit(X_xp[,i]))

X_xp <- t(t(X_xp) - colMeans(X_xp))

Data = list( n_i=n_i, 
             c_i=c_i, 
             s_i=s_i, 
             t_i=t_i, 
             X_xp=as.matrix(X_xp))

Parameters = list(beta_c = rep(0,ncol(Data$X_xp)),
                  beta_p = rep(0,ncol(Data$X_xp)),
                  eps_c_s = unique(Data$s_i)*0,
                  eps_p_s = unique(Data$s_i)*0,
                  fs_c_sd = 0,
                  fs_p_sd = 0,
                  fsd = 10)

Random = c("eps_c_s","eps_p_s")

myMap <- list()
UseP_Map <- FALSE
if(UseP_Map){
  myMap <- list(eps_p_s=as.factor(rep(NA,length(Parameters$eps_p_s)))
                ,eps_c_s=as.factor(rep(NA,length(Parameters$eps_c_s)))
                ,fs_p_sd = as.factor(NA)
                ,fs_c_sd = as.factor(NA)
  )
}

# Make object
Obj = MakeADFun(data=Data,
                parameters=Parameters,
                random=Random,
                map = myMap,
                hessian=FALSE,
                DLL="squid_v1_delta_glm")


out <- nlminb(Obj$par,Obj$fn,Obj$gr)
# SD <- sdreport(Obj)
# 
#This is the output of the model
rep <- Obj$report()

