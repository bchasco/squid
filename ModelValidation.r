ParHat <- fit$ParHat

#7) Temporal correlation
RhoConfig= c("Beta1" = 0 #Temporal corr. encounter covariate intercepts
             ,"Beta2" = 0 #Temporal corr. for positive catch covariate intercepts
             ,"Epsilon1"= 0 #Temporal corr. for encounter probability intercepts
             ,"Epsilon2" = 0) #Temporal corr. for positive catch intercepts

settings <- make_settings(
  n_x = n_x
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



# Generate partitions in data
n_fold <- 10
Partition_i <- sample( 1:n_fold, size=nrow(raw), replace=TRUE )
prednll_f <- rep(NA, n_fold )

# Loop through partitions, refitting each time with a different PredTF_i
n_i <- nrow(raw)
for( fI in 1:n_fold ){
  rm(fit_new)
  
  PredTF_i = ifelse( Partition_i==fI, TRUE, FALSE )
  
  fit_new <- fit_model(settings = settings
                   ,Lat_i = raw$Lat
                   ,Lon_i = raw$Lon
                   ,t_i = raw$Year
                   ,c_i = c_iz
                   ,b_i =  b_i #Number of squid captured.
                   ,a_i = a_i
                   ,v_i = rep(0,n_i)
                   ,"PredTF_i"=PredTF_i
                   ,"Parameters"=ParHat
                   ,"getsd"=FALSE
  )
  
  # Save fit to out-of-bag data
  prednll_f[fI] = fit_new$Report$pred_jnll
}
