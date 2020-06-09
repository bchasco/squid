library(VAST)
library(TMB)
#Data
raw <- read.csv("comb_catches.csv")

#Decisions for VAST analysis. This follows the structure of Thorson (2019)
#1)A multi region input looks like this
strata.limits <- data.frame(
  'STRATA' = c("Coastwide","CA","OR","WA"),
  'north_border' = c(49.0, 42.0, 46.0, 49.0),
  'south_border' = c(37.0, 37.0, 42.0, 46.0)
)

# strata.limits <- data.frame(
#   'STRATA' = c("Coastwide"),
#   'north_border' = c(49.0),
#   'south_border' = c(32.0)
# )


#2) Single size class for now
c_iz <- rep(0,nrow(raw))

#3) CPUE for now must change it to numbers
b_i <- raw$CPUE

#4) Spatial and Spatio-temporal for both encounter and positive catches
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
FieldConfig <-c(Omega1 = 1
                ,Epsilon1 = 1
                ,Omega2 = 1
                ,Epsilon2 = 1)

#7) Temporal correlation
RhoConfig= c("Beta1" = 0
             ,"Beta2" = 0
             ,"Epsilon1"= 1
             ,"Epsilon2" = 1)

#8) Density dependent covariates
#None for now

#9) Catchability associated with vessel
Q_ik <- raw[,c('X3m_Temp','X3m_Salinity','X3m_Chl')] #rep(1,nrow(raw))
for(i in 1:ncol(Q_ik)){
  Q_ik[is.na(Q_ik[,i]),i] <- mean(na.omit(Q_ik[,i]))
}

#10) Area offsets
a_i <- rep(1,nrow(raw))

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
fit <- fit_model(settings = settings
                 ,Lat_i = raw$Lat
                 ,Lon_i = raw$Long
                 ,t_iz = raw$Year
                 ,c_iz = c_iz
                 ,b_i = raw$CPUE #Number of squid captured.
                 ,a_i = rep(0.01,n_i)
                 ,v_i = rep(1,n_i)
                 ,Q_ik = as.matrix(Q_ik)
)

# plot_dir <- paste0(getwd(),"/VAST_plots/")
# plot_ids <- c(1,3,6,7,11,13,14) 
# for(i in plot_ids){
#   plot_results(fit=fit
#                ,working_dir = plot_dir
#                ,plot_set = i)
# }
# 
