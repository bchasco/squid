library(INLA); library(sp); library(fields)
library(geoR)
library(viridisLite)
library(TMB)
library(ggthemes)
library(maps)
require(sp)
require(maptools)
require(ggplot2)

survey <- read.csv("SquidData.csv", sep=",", header = TRUE)
survey <- survey[!is.na(survey$Distance_Offshore_nmi),]
survey <- survey[survey$Study_Type=="Regular",]
stations <- unique(survey$Station)
# for(i in stations){
#   survey$Lat[survey$Station==i] <- mean(survey$Lat[survey$Station==i])
#   survey$Lon[survey$Station==i] <- mean(survey$Lon[survey$Station==i])
# 
# }

df = data.frame(y = log(survey$California_market_squid_), 
                locx = survey$Lon, 
                locy = survey$Lat)
# spatial.scaling = 1
# df$locx = (df$locx - min(df$locx))/spatial.scaling
# df$locy = (df$locy - min(df$locy))/spatial.scaling
df$dist = survey$Distance_Offshore_nmi - mean(survey$Distance_Offshore_nmi)
df$y = df$y-min(df$y)

# quilt.plot(x=df$locx,y=df$locy,z=df$y,nx=40,ny=40, 
#            col = plasma(101),
#            main = "Data")

max.edge = 0.99
mesh = inla.mesh.2d(loc=cbind(df$locx, df$locy),
                     max.edge = c(0.9,0.9), cutoff = 0.1)
try(dyn.unload("squid"))
# compile("squid.cpp")
dyn.load("squid")

n_i <- nrow(survey)
n_x <- mesh$n
n_t <- length(unique(survey$Year))
n_p <- 1 #number of covariates

x_s <- mesh$idx$loc - 1
c_i <- survey$California_market_squid_
s_i <- data.frame(Station=survey$Station)
si_lu <- data.frame(s_i=1:length(unique(survey$Station)),Station=unique(survey$Station))
s_i <- merge(s_i,si_lu)$s_i-1

t_i <- survey[,'Year']-min(survey[,'Year']) 
X_xp <- as.matrix(rep(1,n_x),n_x,1)

# prior.median.sd = 10; prior.median.range = 7
# - diff(range(mesh$loc[, 1]))/2
# - sd(df$y)/10
# spde = inla.spde2.pcmatern(mesh,
#                            prior.range = c(prior.median.range, .5),
#                            prior.sigma = c(prior.median.sd, .5),
#                            constr = T)

spde = inla.spde2.matern(mesh)

Data = list( n_i=n_i, 
             n_x=n_x, 
             n_t=n_t, 
             n_p=ncol(X_xp),
             x_s=x_s, 
             c_i=c_i, 
             s_i=s_i, 
             t_i=t_i, 
             X_xp=X_xp, 
             G0=spde$param.inla$M0, 
             G1=spde$param.inla$M1, 
             G2=spde$param.inla$M2)

Parameters = list(fp= 0,
                  alpha=rep(0.0), 
                  fsd = 5.0,
                  phi=0.0, 
                  log_tau_E=1.0, 
                  log_tau_O=1.0, 
                  log_kappa=0.0,	
                  rho=0.5, 
                  fs_sd = 0,
                  eps_s = unique(Data$s_i)*0,
                  eps_p = unique(Data$s_i)*0,
                  fp_sd = 0,
                  Epsilon_input=matrix(0,nrow=mesh$n,ncol=Data$n_t), 
                  Omega_input=rep(0,mesh$n))

Random = c("Epsilon_input","Omega_input", "eps_s","eps_p")

# library(TMBhelper)
myMap <- list()
UseMap <- FALSE
if(UseMap){
  myMap <- list(Epsilon_input=as.factor(matrix(NA,nrow=mesh$n,ncol=Data$n_t))
                ,Omega_input=as.factor(rep(NA,mesh$n))
                ,phi = as.factor(NA)
                ,log_tau_E = as.factor(NA)
                ,log_tau_O = as.factor(NA)
                ,log_kappa = as.factor(NA)
                ,rho = as.factor(NA)
  )
}

UseP_Map <- TRUE
if(UseP_Map){
  myMap <- list(eps_p=as.factor(rep(NA,length(Parameters$eps_p)))
                ,eps_s=as.factor(rep(NA,length(Parameters$eps_s)))
                ,fp_sd = as.factor(NA)
                ,fs_sd = as.factor(NA)
  )
}

# Make object
Obj = MakeADFun(data=Data,
                parameters=Parameters,
                random=Random,
                map = myMap,
                hessian=FALSE,
                DLL="squid")


out <- nlminb(Obj$par,Obj$fn,Obj$gr)
# SD <- sdreport(Obj)
# 
#This is the output of the model
rep <- Obj$report()

