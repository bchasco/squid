library(TMB)
library(VAST)

AICoutput <- data.frame(AIC=NA,survey=NA,q=NA,k=NA,k2=NA)

#Initial the parameter lists for the models you are testing
ParHatList <- list()
icnt <- 1
for(mySurvey in c(0,1)){ #If 1 then make and offset for survey type
  for(q_d in c(0,1)){
    
    for(k in c(0,1,2,3)){
      for(k2 in c(1)){

        #Data
        raw <- read.csv("Update_Comb_Catch_wTrawlDist_flat.csv")
        # if(subDat!="All"){
        #   raw <- raw[raw$Survey==subDat,]
        # }
        
        
        #Adjust by wainright paper for JSOES
        raw$catch[raw$Gear=="264NRT+MMED_Down" & raw$Extension=="No"] <- raw$catch[raw$Gear=="264NRT+MMED_Down" & raw$Extension=="No"]/0.48
        raw$catch[raw$Gear=="264NRT+MMED_Up" & raw$Extension=="No"] <- raw$catch[raw$Gear=="264NRT+MMED_Up" & raw$Extension=="No"]/0.88
        
        #Adjustments for SWFSC - REad Cheryl's email from 9/15/2020
        raw$catch[raw$Gear=="264NRT+MMED modified" & raw$Extension=="No"] <- raw$catch[raw$Gear=="264NRT+MMED modified" & raw$Extension=="No"]/0.48
        raw$catch[raw$Gear=="264NRT+MMED" & raw$Extension=="No"] <- raw$catch[raw$Gear=="264NRT+MMED" & raw$Extension=="No"]/0.88
        
        
        #Get rid of any blanks
        raw <- raw[!apply(raw,1,function(x)return(sum(is.na(x)))),]
        

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
        
        #Decisions for VAST analysis. This follows the structure of Thorson (2019)
        #1)A multi region input looks like this
        strata.limits <- data.frame(
          'STRATA' = c("Coastwide","CA","OR","WA"),
          'north_border' = c(49.0, 42.0, 46.0, 49.0),
          'south_border' = c(37.0, 37.0, 42.0, 46.0)
        )
        
        #2) Single size class for now
        c_iz <- rep(0,dim(raw)[1]) #This needs to be numeric starting at 0
        # c_iz <- as.numeric(raw$Survey)-1
        
        #3) CPUE for now must change it to numbers
        b_i <- raw$catch
        
        
        Mesh.Method <- "samples" #mesh
        grid_size_km <- 100
        n_x <- 175 #number of knots.  This is really important. To few and the model won't converge, too many and it will take forever.
        Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e1 )
        Aniso <- FALSE #isotropic
        
        #10) Area offsets
        a_i <- raw$TrawlDist_km*0.02 #distance times net width, see Cheryl's email from Aug. 10th.
        
        #12)The observation model is
        ObsModel = c(2,0) # Distribution for data, and link-function for linear predictors
        
        
        FieldConfig <- rep(1,4)
        
        RhoConfig <- c("Beta1" = k #Temporal corr. encounter covariate intercepts
                       ,"Beta2" = k #Temporal corr. for positive catch covariate intercepts
                       ,"Epsilon1"= k2 #Temporal corr. for encounter probability intercepts
                       ,"Epsilon2" = k2) #Temporal corr. for positive catch intercepts
        
        settings <- make_settings(
          n_x = n_x
          ,Region = "california_current"
          ,purpose = "index"
          ,strata.limits = strata.limits  
          ,FieldConfig = FieldConfig
          ,RhoConfig = RhoConfig
          # ,OverdispersionConfig = OverdispersionConfig
          ,ObsModel = ObsModel
          ,knot_method = "samples"
          ,bias.correct = FALSE
          # ,Options = Options
        )
        
        
        if(q_d==0 & mySurvey==0){
          fit <- tryCatch(fit_model(settings = settings
                                    ,Lat_i = raw$Lat
                                    ,Lon_i = raw$Lon
                                    ,t_i = raw$Year
                                    ,b_i =  b_i #Number of squid captured.
                                    ,a_i = a_i
                                    ,silent = FALSE
                                    ,getsd=FALSE
          ),error=function(e) NULL) 
        }
        if(q_d==1 | mySurvey==1){
          fit <- tryCatch(fit_model(settings = settings
                                    ,Lat_i = raw$Lat
                                    ,Lon_i = raw$Lon
                                    ,t_i = raw$Year
                                    ,b_i =  b_i #Number of squid captured.
                                    ,a_i = a_i
                                    ,silent = FALSE
                                    ,Q_ik = as.matrix(Q_ik)
                                    ,getsd = FALSE
          ),error=function(e) NULL) 
        }
        
        if(!is.null(fit$parameter_estimates$AIC))
          AICoutput[icnt,] <- c(fit$parameter_estimates$AIC,mySurvey,q_d,k,k2)
        if(is.null(fit$parameter_estimates$AIC))
          AICoutput[icnt,] <- c(NA,mySurvey,q_d,k,k2)
        icnt <- icnt + 1
      }
      
    }
  }
}
  

# save(AICoutput,file="log_AICoutput.rData")



  
