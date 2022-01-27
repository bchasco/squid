library(VAST)
library(TMB)
library(dplyr)

#Data
raw <- read.csv("Update_Comb_Catch_wTrawlDist_flat.csv")
# raw <- read.csv("Update_Comb_Catch_wTrawlDist.csv")

#Adjust by wainright paper
raw$catch[raw$Gear=="264NRT+MMED_Down" & raw$Extension=="No"] <- raw$catch[raw$Gear=="264NRT+MMED_Down" & raw$Extension=="No"]/0.48
raw$catch[raw$Gear=="264NRT+MMED_Up" & raw$Extension=="No"] <- raw$catch[raw$Gear=="264NRT+MMED_Up" & raw$Extension=="No"]/0.88
  
#Adjustments for SWFSC - REad Cheryl's email from 9/15/2020
raw$catch[raw$Gear=="264NRT+MMED modified" & raw$Extension=="No"] <- raw$catch[raw$Gear=="264NRT+MMED modified" & raw$Extension=="No"]/0.48
raw$catch[raw$Gear=="264NRT+MMED" & raw$Extension=="No"] <- raw$catch[raw$Gear=="264NRT+MMED" & raw$Extension=="No"]/0.88

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
RhoConfig= c("Beta1" = 0 #Temporal corr. encounter covariate intercepts
               ,"Beta2" = 0 #Temporal corr. for positive catch covariate intercepts
               ,"Epsilon1"= 0 #Temporal corr. for encounter probability intercepts
               ,"Epsilon2" = 0) #Temporal corr. for positive catch intercepts

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
ObsModel = c(1,0) # Distribution for data, and link-function for linear predictors

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
  # ,n_g = Inf
  ,strata.limits = strata.limits
  ,FieldConfig = FieldConfig
  ,RhoConfig = RhoConfig
  # ,OverdispersionConfig = OverdispersionConfig
  ,ObsModel = ObsModel
  ,knot_method = "samples"
  ,bias.correct = FALSE
  ,fine_scale = FALSE
  # ,Options = Options
)

############Putting is all together
# try(rm(fit))

AICTable <- data.frame(Model=NA,AIC=NA)
ii <- 1
for(k_sub in 1:1){
  for(q in 0:1){#0:1include environmental covariates
    for(k in c(2)){#c(0,1,2,4) temporal correlation
      for(i in c(0)){#1:4 number of spatial and spatio-temporal flags turned on
        myComb <- combn(1:4,i)
        for(F_ar in c(1)){#c(1,"AR1") #AR1 does not appear to converge for any of the models
          for(j in 1:ncol(myComb)){ #model combinations of spatial and spatiotemporal flags
            FieldConfig[1:4] <- rep(F_ar,4)
            FieldConfig[myComb[,j]] <- 0
            
            settings$FieldConfig <- FieldConfig
            if(k_sub==1){
              settings$RhoConfig <- c("Beta1" = 1 #Temporal corr. encounter covariate intercepts
                                      ,"Beta2" = 1 #Temporal corr. for positive catch covariate intercepts
                                      ,"Epsilon1"= 2 #Temporal corr. for encounter probability intercepts
                                      ,"Epsilon2" = 2) #Temporal corr. for positive catch intercepts
              
            }
            if(k_sub==2){
              settings$RhoConfig <- c("Beta1" = 0 #Temporal corr. encounter covariate intercepts
                                      ,"Beta2" = 0 #Temporal corr. for positive catch covariate intercepts
                                      ,"Epsilon1"= k #Temporal corr. for encounter probability intercepts
                                      ,"Epsilon2" = k) #Temporal corr. for positive catch intercepts
            }
            if(k_sub==3){
              settings$RhoConfig <- c("Beta1" = k #Temporal corr. encounter covariate intercepts
                                      ,"Beta2" = k #Temporal corr. for positive catch covariate intercepts
                                      ,"Epsilon1"= k #Temporal corr. for encounter probability intercepts
                                      ,"Epsilon2" = k) #Temporal corr. for positive catch intercepts
            }
            
            if(q==1){
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
                                        # ,getsd = TRUE
                                        # ,working_dir = getwd()
                                        # ,formula=formula
                                        # ,covariate_data=covariate_data
              ),error=function(e) NULL) 
            }
            if(q==0){
              fit <- tryCatch(fit_model(settings = settings
                                        ,Lat_i = raw$Lat
                                        ,Lon_i = raw$Lon
                                        ,t_i = raw$Year
                                        # ,c_i = c_iz
                                        ,b_i =  b_i #Number of squid captured.
                                        ,a_i = a_i
                                        # ,v_i = rep(0,nrow(raw))
                                        ,silent = TRUE
                                        # ,getsd = TRUE
                                        # ,Q_ik = as.matrix(Q_ik)
                                        # ,formula=formula
                                        # ,covariate_data=covariate_data
              ),error=function(e) NULL) 
            }
            
            print(paste(c(q,k_sub,k,FieldConfig),collapse=""))
            
            AICTable[ii,] <- c(paste(c(q,k_sub,k,FieldConfig),collapse=""),
                               fit$parameter_estimates$AIC)
            ii <- ii + 1
            
          }
        }
      }
    }  
  }
}


# write.csv(AICTable,file="AICTable.csv")
# source("modelvalidation.r")