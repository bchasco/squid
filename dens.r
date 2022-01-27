library(TMB)
load("dens.rData")
try(dyn.unload(dynlib("dens")))
compile("dens.cpp")
dyn.load(dynlib("dens"))

data <- list(ef = ef[1:22],
             idx=idx[1:22],
             ef_sd=ef_sd[1:22],
             idx_sd=idx_sd[1:22])

parameters <- list(u=rep(0,length(data$ef)),
                   a=0, 
                   b=0)


obj <- MakeADFun(data, 
                 parameters, 
                 random=c('u'),
                 DLL="dens")
system.time(opt<-nlminb(obj$par,obj$fn,obj$gr))

opt
sd <- sdreport(obj)
