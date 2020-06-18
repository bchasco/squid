library(VAST)
library(TMB)
library(dplyr)

#Data
# raw <- read.csv("Update_Comb_Catch_wTrawlDist_flat.csv")
raw <- read.csv("Update_Comb_Catch_wTrawlDist.csv")


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
# b_i <- raw$Total_CPUE
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
Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )
Aniso <- TRUE #isotropic

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
               ,"Epsilon1"= 4 #Temporal corr. for encounter probability intercepts
               ,"Epsilon2" = 4) #Temporal corr. for positive catch intercepts

#8) Density dependent covariates
#None for now
# install.packages("splines")
library(splines)
formula = ~ bs( log(Station_Depth_m), knots=3, intercept=FALSE)

# covariate_data <- data.frame(Lat=raw$Lat,Lon=raw$Lat,Year=raw$Year,Station_Depth_m=raw$Station_Depth_m/100)
covariate_data <- data.frame(Lat=raw$Lat,Lon=raw$Lat,Year=NA,Station_Depth_m=raw$Station_Depth_m/100)
covariate_data$Year <- NA

#9) Catchability associated with vessel
Q_ik <- raw[,c('X3m_Temp','X3m_Salinity','X3m_Chl')] #rep(1,nrow(raw))
for(i in 1:ncol(Q_ik)){
  Q_ik[is.na(Q_ik[,i]),i] <- mean(na.omit(Q_ik[,i]))
}

#10) Area offsets
a_i <- raw$TrawlDist_km*0.085 #distance times net width

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
            ,Calculate_effective_area = 1
            ,Calculate_Cov_SE = 0 
            ,Calculate_Synchrony = 0
            ,Calculate_Coherence = 0)
#14)

#15)

############Putting is all together

AICTable <- data.frame(Model=NA,AIC=NA)
ii <- 1
# for(k in 0:1){
#   for(i in 1:4){
#     myComb <- combn(1:4,i)
#     for(j in 1:ncol(myComb)){
      FieldConfig[1:4] <- rep(1,4)
      FieldConfig[myComb[,j]] <- 0
      # if(sum(FieldConfig[1:2]>0)){
      settings <- make_settings(
        n_x = n_x
        # ,grid_size_km = grid_size_km
        # ,randomseed = Kmeans_Config[["randomseed"]]
        # ,nstart = Kmeans_Config[["nstart"]]
        # ,iter.max = Kmeans_Config[["iter.max"]]
        ,Region = "california_current"
        ,purpose = "index2"
        ,strata.limits = strata.limits  
        ,FieldConfig = FieldConfig
        ,RhoConfig = RhoConfig
        ,OverdispersionConfig = OverdispersionConfig
        ,ObsModel = ObsModel
        ,knot_method = "samples"
        ,bias.correct = FALSE
        ,Options = Options
      )
      
      
      n_i <- nrow(raw)
      if(k==0){
        fit <- fit_model(settings = settings
                         ,Lat_i = raw$Lat
                         ,Lon_i = raw$Lon
                         ,t_iz = raw$Year
                         ,c_iz = c_iz
                         ,e_i = c_iz
                         ,b_i =  b_i #Number of squid captured.
                         ,a_i = a_i
                         ,v_i = rep(0,n_i)
                         ,formula=formula
                         ,covariate_data=covariate_data
        )
      }
      if(k==1){
        fit <- fit_model(settings = settings
                         ,Lat_i = raw$Lat
                         ,Lon_i = raw$Lon
                         ,t_iz = raw$Year
                         ,c_iz = c_iz
                         ,e_i = c_iz
                         ,b_i =  b_i #Number of squid captured.
                         ,a_i = a_i
                         ,v_i = rep(0,n_i)
                         ,Q_ik = as.matrix(Q_ik)
                         ,formula=formula
                         ,covariate_data=covariate_data
        )
      }
      
      # }
      
      AICTable[ii,] <- c(paste(c(k,FieldConfig),collapse=""),
                         fit$parameter_estimates$AIC[1])
      ii <- ii + 1
      
#     }
#   }
# }


plot_dir <- paste0(getwd(),"/VAST_plots/")
plot_ids <- c(1,3,6,7,11,13,14)
for(i in plot_ids){
  plot_results(fit=fit
               # ,years_to_plot = c(3,4)
               ,working_dir = plot_dir
               ,plot_set = i)
}
# 
# Cov_List = summarize_covariance(Report = fit$Report, 
#                                 ParHat = fit$tmb_list$Obj$env$parList(),
#                                 Data = fit$data_list, 
#                                 SD = fit$Opt$SD, 
#                                 plot_cor = TRUE,
#                                 # category_names = c("1","2"),
#                                 plotTF = FieldConfig,
#                                 mgp = c(2,0.5, 0), 
#                                 tck = -0.02, 
#                                 oma = c(0, 5, 2, 2))
# 
