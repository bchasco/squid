library(VAST)
library(TMB)
library(dplyr)

#Data
raw <- read.csv("Update_Comb_Catch_wTrawlDist_flat.csv")
# raw <- read.csv("Update_Comb_Catch_wTrawlDist.csv")

#Adjust by wainright paper
raw$catch[raw$Gear=="264NRT+MMED_Down" & raw$Extension=="No"] <- raw$catch[raw$Gear=="264NRT+MMED_Down" & raw$Extension=="No"]/0.48
raw$catch[raw$Gear=="264NRT+MMED_Up" & raw$Extension=="No"] <- raw$catch[raw$Gear=="264NRT+MMED_Up" & raw$Extension=="No"]/0.88

#Get rid of any blanks
raw <- raw[!apply(raw,1,function(x)return(sum(is.na(x)))),]

#Decisions for VAST analysis. This follows the structure of Thorson (2019)
#1)A multi region input looks like this
strata.limits <- data.frame(
  'STRATA' = c("Coastwide","CA","OR","WA"),
  'north_border' = c(49.0, 42.0, 46.0, 49.0),
  'south_border' = c(37.0, 37.0, 42.0, 46.0)
)

#2) Single size class for now
c_iz <- rep(0,dim(raw)[1]) #This needs to be numeric starting at 0
# c_iz <- as.integer(as.factor(raw$size))-1 #This needs to be numeric starting at 0

#3) CPUE for now must change it to numbers
b_i <- raw$catch

#4) Spatial and Spatio-temporal for both encounter and positive catches, see #6)
FieldConfig = c("Omega1" = 0
                ,"Epsilon1" = 0
                ,"Omega2" = 0
                ,"Epsilon2" = 0)

#5) Choosing the spatial
Mesh.Method <- "samples" #mesh
grid_size_km <- 100
n_x <- 175 #number of knots.  This is really important. To few and the model won't converge, too many and it will take forever.
Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e1 )
Aniso <- FALSE #isotropic

#6) Choosing the number of spatial and temporal factors
#We have a single factor to the encounter probabilites (Omega1 and Epsilon1)
#We have a single factor to the positive encounters (Omega2 and Epsilon2)
#1 is the encounter, 2 is the positive catches
#Omega is spatial, Epsilon is spatio-temporal
FieldConfig <-c(Omega1 = 1 #Spatial corr. encounter probability
                ,Epsilon1 = 1 #Spatial corr for positive catch probability.
                ,Omega2 = 1 #Spatiotemporal corr. for encounter probability 
                ,Epsilon2 = 1)#Spatiotemporal corr. for positive catch probability

#7) Temporal correlation
RhoConfig= c("Beta1" = 3 #Temporal corr. encounter covariate intercepts
             ,"Beta2" = 3 #Temporal corr. for positive catch covariate intercepts
             ,"Epsilon1"= 4 #Temporal corr. for encounter probability intercepts
             ,"Epsilon2" = 4) #Temporal corr. for positive catch intercepts

#8) Density dependent covariates
#None for now
# install.packages("splines")
library(splines)
formula = ~ bs( log(X3m_Temp), knots=3, intercept=FALSE)

# covariate_data <- data.frame(Lat=raw$Lat,Lon=raw$Lat,Year=raw$Year,Station_Depth_m=raw$Station_Depth_m/100)
covariate_data <- data.frame(Lat=raw$Lat,Lon=raw$Lat,Year=NA,X3m_Temp=(raw$X3m_Temp-mean(raw$X3m_Temp))/sd(raw$X3m_Temp))
covariate_data$Year <- NA

#9) Catchability associated with vessel
Q_ik <- raw[,c('Top20m_Temp','Top20m_Salinity')] #rep(1,nrow(raw))
for(i in 1:ncol(Q_ik)){
  Q_ik[is.na(Q_ik[,i]),i] <- mean(na.omit(Q_ik[,i]))
  Q_ik[,i] <- (Q_ik[,i]-mean(Q_ik[,i]))/sd(Q_ik[,i])
}

#10) Area offsets
a_i <- raw$TrawlDist_km*0.02 #distance times net width, see Cheryl's email from Aug. 10th.

#11) Over dispersion vessel effects
#This treats the vessel and MMED effects as fixed
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)


#12)The observation model is
ObsModel = c(2,0) # Distribution for data, and link-function for linear predictors

#13) Options for derived quantities
Options = c(SD_site_density = 0 
            ,SD_site_logdensity = 0
            ,Calculate_Range = 1 #Center of gravity 
            ,Calculate_evenness = 0 
            ,Calculate_effective_area = 0
            ,Calculate_Cov_SE = 0 
            ,Calculate_Synchrony = 0
            ,Calculate_Coherence = 0)
#14)

#15)

settings <- make_settings(
  n_x = n_x
  ,Region = "california_current"
  ,purpose = "index"
  ,strata.limits = strata.limits
  ,FieldConfig = FieldConfig
  ,RhoConfig = RhoConfig
  ,ObsModel = ObsModel
  ,knot_method = "samples"
  ,bias.correct = FALSE
  ,fine_scale = FALSE
)

FieldConfig <- rep(1,4)

settings <- make_settings(
  n_x = n_x
  ,Region = "california_current"
  ,purpose = "index"
  ,strata.limits = strata.limits  
  ,FieldConfig = FieldConfig
  ,RhoConfig = RhoConfig
  ,ObsModel = ObsModel
  ,knot_method = "samples"
  ,bias.correct = FALSE
)

#Initial the parameter lists for the models you are testing
ParHatList <- list()
yrAIC <- NA
icnt <- 1


for(mySurvey in c(0,1)){ #If 1 then make and offset for survey type
  for(q_d in c(1)){

      Q_ik <- raw[,c('Top20m_Temp','Top20m_Salinity')] #rep(1,nrow(raw))
      for(i in 1:ncol(Q_ik)){
        Q_ik[is.na(Q_ik[,i]),i] <- mean(na.omit(Q_ik[,i]))
        Q_ik[,i] <- (Q_ik[,i]-mean(Q_ik[,i]))/sd(Q_ik[,i])
      }
      if(mySurvey==1){
        if(q_d==1){
          Q_ik[,ncol(Q_ik)+1] <- 0
          Q_ik[raw$Survey=="JSOES",ncol(Q_ik)] <- 1
          
        }
        if(q_d==0){
          Q_ik <- rep(0,nrow(raw))
          Q_ik[raw$Survey=="JSOES"] <- 1
        }
      }

    for(k in c(2)){
      for(k2 in c(1)){
        # k <- 
        settings$RhoConfig <- c("Beta1" = k #Temporal corr. encounter covariate intercepts
                                ,"Beta2" = k #Temporal corr. for positive catch covariate intercepts
                                ,"Epsilon1"= k2 #Temporal corr. for encounter probability intercepts
                                ,"Epsilon2" = k2) #Temporal corr. for positive catch intercepts
        
        if(q_d==0 & mySurvey==0){
          fit <- tryCatch(fit_model(settings = settings
                                    ,Lat_i = raw$Lat
                                    ,Lon_i = raw$Lon
                                    ,t_i = raw$Year
                                    # ,c_i = c_iz
                                    ,b_i =  b_i #Number of squid captured.
                                    ,a_i = a_i
                                    # ,v_i = rep(0,nrow(raw))
                                    ,silent = TRUE
                                    # ,Q_ik = as.matrix(Q_ik)
                                    # ,working_dir = getwd()
                                    # ,formula=formula
                                    # ,covariate_data=covariate_data
          ),error=function(e) NULL) 
        }
        
        if(q_d==1 | mySurvey==1){
          fit <- tryCatch(fit_model(settings = settings
                                    ,Lat_i = raw$Lat
                                    ,Lon_i = raw$Lon
                                    ,t_i = raw$Year
                                    # ,c_i = c_iz
                                    ,b_i =  b_i #Number of squid captured.
                                    ,a_i = a_i
                                    # ,v_i = rep(0,nrow(raw))
                                    ,silent = TRUE
                                    ,Q_ik = as.matrix(Q_ik)
                                    # ,working_dir = getwd()
                                    # ,formula=formula
                                    # ,covariate_data=covariate_data
          ),error=function(e) NULL) 
        }
        
        ParHatList[[icnt]] <- fit$ParHat
        if(!is.null(fit$parameter_estimates$AIC)){
          yrAIC[icnt] <- fit$parameter_estimates$AIC
        }
        icnt <- icnt + 1
      }
    }
  }
}

#Run mulitple cross-validations for each mode using 10-folds.  
nsim <- 50
n_fold <- 10
prednll_f <- array(NA, c(nsim,2,n_fold))
nullCnt <- 1
for(ii in 1:nsim){ #number of simulations
  # Generate partitions in data
  for( fI in 1:n_fold){#n_fold #create the folds
    Partition_i <- sample( 1:n_fold, size=nrow(raw), replace=TRUE )
    PredTF_i = ifelse(Partition_i==fI, TRUE, FALSE )
    icnt <- 1
    for(mySurvey in c(0,1)){ #If 1 then make and offset for survey type
      rm(fit_new)
      for(q_d in c(1)){ #change the covariates
        
        Q_ik <- raw[,c('Top20m_Temp','Top20m_Salinity')] #rep(1,nrow(raw))
        for(i in 1:ncol(Q_ik)){
          Q_ik[is.na(Q_ik[,i]),i] <- mean(na.omit(Q_ik[,i]))
          Q_ik[,i] <- (Q_ik[,i]-mean(Q_ik[,i]))/sd(Q_ik[,i])
        }
        if(mySurvey==1){
          if(q_d==1){
            Q_ik[,ncol(Q_ik)+1] <- 0
            Q_ik[raw$Survey=="JSOES",ncol(Q_ik)] <- 1
            
          }
          if(q_d==0){
            Q_ik <- rep(0,nrow(raw))
            Q_ik[raw$Survey=="JSOES"] <- 1
          }
        }
        
        for(k in c(2)){ #with or without autocorrelation
          for(k2 in c(1)){ #with or without autocorrelation
            settings$RhoConfig <- c("Beta1" = k #Temporal corr. encounter covariate intercepts
                                    ,"Beta2" = k #Temporal corr. for positive catch covariate intercepts
                                    ,"Epsilon1"= k2 #Temporal corr. for encounter probability intercepts
                                    ,"Epsilon2" = k2) #Temporal corr. for positive catch intercepts
            if(q_d==0 & mySurvey==0){
              # Loop through partitions, refitting each time with a different PredTF_i
              fit_new <- tryCatch(fit_model(settings = settings
                                            ,Lat_i = raw$Lat
                                            ,Lon_i = raw$Lon
                                            ,t_i = raw$Year
                                            # ,c_i = c_iz
                                            ,b_i =  b_i #Number of squid captured.
                                            ,a_i = a_i
                                            # ,v_i = rep(0,nrow(raw))
                                            ,silent = TRUE
                                            # ,Q_ik = as.matrix(Q_ik)
                                            ,"PredTF_i"=PredTF_i
                                            ,"Parameters"=ParHatList[[icnt]]
                                            ,"getsd"=FALSE
              ),error=function(e) NULL)
            }
            if(q_d==1 | mySurvey==1){
              # Loop through partitions, refitting each time with a different PredTF_i
              fit_new <- tryCatch(fit_model(settings = settings
                                            ,Lat_i = raw$Lat
                                            ,Lon_i = raw$Lon
                                            ,t_i = raw$Year
                                            # ,c_i = c_iz
                                            ,b_i =  b_i #Number of squid captured.
                                            ,a_i = a_i
                                            # ,v_i = rep(0,nrow(raw))
                                            ,silent = TRUE
                                            ,Q_ik = as.matrix(Q_ik)
                                            ,"PredTF_i"=PredTF_i
                                            ,"Parameters"=ParHatList[[icnt]]
                                            ,"getsd"=FALSE
              ),error=function(e) NULL) 
              
            }
            # Save fit to out-of-bag data
            if(!is.null(fit_new$Report$pred_jnll)){
              # prednll_f[ii,icnt,fI-min(raw$Year)+1] = fit_new$Report$pred_jnll
              prednll_f[ii,icnt,fI] = fit_new$Report$pred_jnll
              nullCnt <- nullCnt + 1
            }
            print(paste("############## You are on number ", ii, icnt, k, k2, fI,nullCnt))
            print(prednll_f[ii,,fI])
            icnt <- icnt + 1
          }
        }
      }#q_d
    }#mySurvey
  }#fI
}#sim






  
