
#Start by bringing in the flat data.
raw <- read.csv("comb_catches.csv")

#Decisions for VAST analysis. This follows the structure of Thorson (2019)
#1) Define the grid
#To begin with we will start with one region
#We can expand this to three regions, WA, OR, CA

#-A single region input looks like this

strata.limits <- data.frame(STRATA = c("All_areas"))

#-A multi region input looks like this
strata.limits <- data.frame(
  'STRATA' = c("Coastwide","CA","OR","WA"),
  'north_border' = c(49.0, 42.0, 46.0, 49.0),
  'south_border' = c(32.0, 32.0, 42.0, 46.0)
)


#2) Species/size/age/stage
#To begin with we assume that we only have a single size class for squid.

#For the time being the is only one category,
#In the near future, we will change this include, small and large squid categories

c_i <- rep(0,nrow(raw))

#3) Sampling type, encounter, abundance, or biomass

#For now we're going to use cpueN
#The vector for include the catches is b_i

b_i <- raw$CPUE

#4) Spatial and spatiotemporal variation
#There are two linear parts to the VAST model: the encounter probability and the 
#the positive catches. To determine which of those linear predictors 
#have spatial and spatio-temporal components we modify the FieldConfig
#The 1 refers to the encounter and the 2 refers to the positive catches.

#This field configuration means that there is both spatial and spatiotemporal variation
#for the positive catches (2), but neither for the encounter (1). 
FieldConfig = c("Omega1" = 0, "Epsilon1" = 0, "Omega2" = 1, "Epsilon2" = 1)

#5) Choosing the spatial
#This is where VAST chooses how to create the grid that is laid over the
#top of the spatial data. The machinery behind this is the INLA package
#but there are several arguments that the high level VAST wrapper passes to the
#the INLA package. These input include, mesh configuration, number of knots,
#and the isotropy (equal decorrelation in every direction) or anisotropy (unequal
#decorrelation). 

#Let's start with the default settings for squid: INLA SPDE mesh, 
#200 knots to start with, and isotropy.

Mesh.Method <- "Mesh" #mesh
grid_size_km <- 100
n_x <- 200 #number of knots
Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )
Aniso <- FALSE #isotropic

#6) Choosing the number of spatial and temporal factors
#If you have multiple factors included in the analysis, you need to decide how those factors
#will be aggregated. This is done with the loadings matrix and is similar to what 
#people know as principle components analysis.

#In this example, the linear predictor for postive catches would have two factors for the spatial
#probability, and a single factor for the spatiotemporal probability 

#We have a single factor to the encounter probabilites (Omega1 and Epsilon1)
#We have a single factor to the positive encounters (Omega2 and Epsilon2)
#1 is the encounter, 2 is the positive catches
#Omega is spatial, Epsilon is spatio-temporal
FieldConfig <-c(Omega1 =1,  Epsilon1 =1,  Omega2 =1,Epsilon2 =1)

#7) Temporal correlation
#Again, '1' is encounter and '2' is for positive catches
#Beta is for the yearly intercepts. If =0 then the intercepts are just fixed effects. 
#If = AR1' and there is temporal correlation auto-regressive process
#Epsilon is for the spatio-temporal part. As a default, each spatial mesh at time t is not temporal correlated with
#the next spatial mesh at time t if set = 0.
RhoConfig= c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0)

#9) Catchability associated with vessel
#This catchability is a little different that you may be used to. This has the dimension of 
#n_i by n_k (number of samples by the number of covariates). These are things that affects
#catchability associated with the vessel or example weather, vessel length, 
Q_ik <- rep(1,nrow(raw))


#10) Area offsets
#This is for the effective area sampled. In the case of the squid data, we already have CPUE
#so we would just set the offset to 1
a_i <- rep(1,nrow(raw))

#Question to ask Jim, should effort offset be set in a_i or Q_ik

#11)

#The catchability associated features not linked to habitat are 
#important for reducing sampling bias. In the case of the squid surveys, 
#these would include changes in things like vessel effects, marine mammal excluder 
#device (MMED),

#There are a number of ways of implementing the effects of catchability
#If the data are categorical (e.g., vessel, MMED), the we can use the following argument
#In this instance, treats the vessel and MMED effects a random 
#OverdispersionConfig = c("Eta1"=1, "Eta2"=1)

#This treats the vessel and MMED effects as fixed
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
#Because the catergories for the things affects catchability are so small (MMED is binary),
#and the number of vessel used is vessel small, we will use the later 
#configuration


#12)
#The observation model is
#The first number is for the postive catches: 2 is for log-normal
#The second number is for the encounter probability: 0 is for conventional delta glmm
ObsModel = c(2,0) # Distribution for data, and link-function for linear predictors

#13)
#Options for derived quantities
#This is where we pick the derived variables we are interested in
#Center of gravity is the obvious options, it looks for distribution shifts
#based on an unbalanced design in the surveys
Options = c(SD_site_density = 0 
            ,SD_site_logdensity = 0
            ,Calculate_Range = 1 
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
  ,Region = "california_current"
  ,purpose = "index2"
  ,strata.limits = strata.limits  #"All_areas"#strata.limits
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
                 ,t_i = raw$Year
                 ,c_i = rep(0,n_i)
                 ,b_i = raw$CPUE
                 ,a_i = rep(0.01,n_i)
                 ,v_i = rep(1,n_i)
)

plot_results(fit=fit
             ,plot_set = 3)
