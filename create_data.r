#Starting with the goolge drive
#https://drive.google.com/drive/folders/0BxEB7Ne107iyX1JvZGpoc2Rnemc
#This data is made using the 2010-2015 STD matrix format, and the ocean trawl stations file
#I collected the data from the matrix file and made the following fields
#cruise,	waypoint,	state	line,	station,	depth, order,	haul,	Lopa

#Lopa is the standardized catch. It has a zero for trawl where no squid were observed
#and a number greater than zero when catches were observed

library(TMB)
lat_index <- 1
x <- read.csv("SW_Data_Combined.csv")
x$Year <- as.numeric(as.character(x$Year))
x <- x[!is.na(x$Year),]

x$TransectLat_i <- as.integer(as.factor(round(x$Latitude,lat_index)))-1
x$TransectLat <- round(x$Latitude,lat_index)

transDist <- sort(unique(round(x$Latitude,lat_index)))

c_i <- x$Lopa
s_i <- as.integer(x$station)-1
ag <- aggregate(c_i, by=list(s_i),sum)
ag$s_i_re <- -1
ii <- 0
for(i in 1:nrow(ag)){
  if(ag$x[i]>0){
    ag$s_i_re[i] <- ii
    ii <- ii +1 
  }
}


Data <- list(c_i = c_i,
             s_i = s_i,
             s_i_re = ag$s_i_re,
             y_i = x$Year-min(x$Year),
             transDist = transDist,
             trans_i = x$TransectLat_i)

Parameters <- list(fp = -0,
                   ln_lam = log(300),
                   ln_theta = log(1),
                   eps_s = rep(0,max(Data$s_i_re)+1),
                   lnsig_eps_s = 0,
                   eps_y = rep(0,length(unique(Data$y))),
                   frho_y = 0,
                   lnsig_eps_y = 0,
                   eps_tr = rep(0,length(unique(Data$trans_i))),
                   frho_tr = 0,
                   lnsig_eps_tr = 0,
                   eps_yt = matrix(0,length(unique(Data$trans_i)),length(unique(Data$y))),
                   frho_yt_1 = 0,
                   frho_yt_2 = 0,
                   lnsig_eps_yt = 0)



try(dyn.unload("squid2"))
compile("squid2.cpp")
dyn.load("squid2")

map_s <- list(eps_s = as.factor(rep(NA,max(Data$s_i_re)+1)),
              lnsig_eps_s = as.factor(NA))

map_tr <- list(eps_tr = as.factor(rep(NA,length(unique(x$TransectLat_i)))),
              lnsig_eps_tr = as.factor(NA),
              frho_tr = as.factor(NA))

map_y <- list(eps_y = as.factor(rep(NA,max(Data$y_i)+1)),
              lnsig_eps_y = as.factor(NA),
              frho_y = as.factor(NA))

map_yt <- list(eps_yt = as.factor(matrix(NA,length(unique(x$TransectLat_i)),max(Data$y_i)+1)),
              lnsig_eps_yt = as.factor(NA),
              frho_yt_1 = as.factor(NA),
              frho_yt_2 = as.factor(NA))

myMap <- list()
est_tr <- 1
est_s <- 1
est_y <- 1
est_yt <- 1

sim_tr <- 1
sim_s <- 1
sim_y <- 1
sim_yt <- 1

runSim <- TRUE

myMap <- list()

if(!est_tr){
  myMap <- append(myMap
                  ,map_tr)
}
if(!est_s){
  myMap <- append(myMap
                  ,map_s)
}
if(!est_y){
  myMap <- append(myMap
                  ,map_y)
}
if(!est_yt){
  myMap <- append(myMap
                  ,map_yt)
}

# myMap <- append(myMap
#                   ,list(ln_lam=as.factor(NA),
#                   ln_theta=as.factor(NA)))

# myMap <- append(myMap
#                 ,list(fp=as.factor(NA)))

Data <- append(Data,
               list(est_tr = est_tr,
               est_y = est_y,
               est_s = est_s,
               est_yt = est_yt))
# Make object
Obj = MakeADFun(data=Data,
                parameters=Parameters,
                random=c("eps_s", #Station
                         "eps_tr",#Latitude
                         "eps_y", #Year
                         "eps_yt")#Spatio-temporal
                ,
                hessian=TRUE,
                map=myMap,
                DLL="squid2")

# Check_Identifiable(Obj)

out <- nlminb(Obj$par,Obj$fn,Obj$gr)
rep <- Obj$report()


i <- 1
if(runSim){
  myMap <- list()
  
  if(!sim_tr){
    myMap <- append(myMap
                    ,map_tr)
  }
  if(!sim_s){
    myMap <- append(myMap
                    ,map_s)
  }
  if(!sim_y){
    myMap <- append(myMap
                    ,map_y)
  }
  if(!sim_yt){
    myMap <- append(myMap
                    ,map_yt)
  }
  
  simout <- replicate(500,{
    simData <- Obj$simulate(par = Obj$env$last.par.best, complete = TRUE)
    simData$est_s <- sim_s
    simData$est_y <- sim_y
    simData$est_tr <- sim_tr
    simData$est_yt <- sim_yt
    
    obj2 <- MakeADFun(simData,
                      Parameters,
                      random=c("eps_s", #Station
                               "eps_tr",#Latitude
                               "eps_y", #Year
                               "eps_yt")#Spatio-temporal
                      ,
                      hessian=TRUE,
                      silent = TRUE,
                      map=myMap,
                      DLL="squid2")
    print(i)
    i <- i+1
    sim_i <- nlminb(obj2$par,obj2$fn,obj2$gr)$par
    # rep <- obj2$report()
    sim_i})
  
  simMelt <- melt(t(simout),measure.vars = colnames(simout),variable.name="pars")
  library(ggplot2)
  
  myLines <- data.frame(Var2=names(Obj$env$last.par.best),
                        value=Obj$env$last.par.best)
  myLines <- myLines[myLines$Var2%in%rownames(simout),]
  
  p <- ggplot(simMelt, aes(value)) +
    facet_wrap(~Var2, scales = "free") +
    geom_histogram()+
    geom_vline(aes(xintercept=value), colour="blue",myLines)
  
  print(p)
}
