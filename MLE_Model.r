library(VAST)
library(TMB)
library(dplyr)


mySurveys <- c("All")
for(subDat in mySurveys){
  #Data
  raw <- read.csv("Update_Comb_Catch_wTrawlDist_flat.csv")

  file <- paste0(subDat,"MLE.Rdata")
  
  
  #Adjust by wainright paper for JSOES
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

  #9) Catchability associated with surveys
  Q_ik <- raw[,c('Top20m_Temp','Top20m_Salinity')] #rep(1,nrow(raw))
  for(i in 1:ncol(Q_ik)){
    Q_ik[is.na(Q_ik[,i]),i] <- mean(na.omit(Q_ik[,i]))
    Q_ik[,i] <- (Q_ik[,i]-mean(Q_ik[,i]))/sd(Q_ik[,i])
  }

  FieldConfig = c("Omega1" = 1 #spatial detection
                  ,"Epsilon1" = 1 #spatiotemporal detection
                  ,"Omega2" = 1 #spatial positive
                  ,"Epsilon2" = 1) #spatiotemporal positive
  
  #3) CPUE for now must change it to numbers
  b_i <- raw$catch
  
  
  Mesh.Method <- "samples" #mesh
  grid_size_km <- 100
  n_x <- 175 #number of knots.  This is really important. To few and the model won't converge, too many and it will take forever.
  Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e1 )

  RhoConfig= c("Beta1" = 2#Temporal corr. encounter covariate intercepts, 0 SWFSC, , 0 or 1 All
               ,"Beta2" = 2 #Temporal corr. for positive catch covariate intercepts, 0 SWFSC, , 0 or 1 All
               ,"Epsilon1"= 1 #Temporal corr. for encounter probability intercepts, 1 SWFSC, , 0 or 1 All
               ,"Epsilon2" = 1) #Temporal corr. for positive catch intercepts, 1 SWFSC, , 0 or 1 All
  


  #10) Area offsets
  a_i <- raw$TrawlDist_km*0.02 #distance times net width, see Cheryl's email from Aug. 10th.
  
  #11) Over dispersion vessel effects
  #This treats the vessel and MMED effects as fixed
  OverdispersionConfig = c("Eta1"=0, 
                           "Eta2"=0)
  
  
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
    ,use_anisotropy = TRUE #The default is true
    ,Options = Options
  )
  
  fit <- #tryCatch(
    fit_model(settings = settings
                            ,Lat_i = raw$Lat
                            ,Lon_i = raw$Lon
                            ,t_i = raw$Year
                            ,b_i =  b_i #Number of squid captured.
                            ,c_i = c_iz
                            ,a_i = a_i
                            ,silent = TRUE
                            ,Q_ik = as.matrix(Q_ik)
                            ,getReportCovariance = TRUE
                            # ,newtonsteps = 2
                            ,getsd = TRUE #TRUE is the default
              )
  assign(paste0(subDat,"_fit"),fit)
  save(list=ls()[ls()==paste0(subDat,"_fit")],file=file)
}
load("AllMLE.rData")
# for(i in c(1,2,3,4,5,6,7)){
#   plot_results(All_fit,plot_set = i, cex.main=1.5, cex.lab=1.5)
# }
# 
# 
# for(i in c(5)){
#   plot_results(All_fit,plot_set = i, years_to_plot = c(1,22)
#   )
# }
# 

