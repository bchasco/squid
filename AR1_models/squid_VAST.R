
#Load the two libraries that you need, TMB for the libraries, and VAST for the high level wrappers
library(TMB)
library(VAST)

Data_Set = c("Chatham_rise_hake"
             , "Iceland_cod"
             , "WCGBTS_canary"
             , "GSL_american_plaice"
             , "BC_pacific_cod"
             , "EBS_pollock"
             , "GOA_Pcod"
             , "GOA_pollock"
             , "GB_spring_haddock"
             , "GB_fall_haddock"
             , "SAWC_jacopever"
             , "Aleutian_islands_POP")[3]

#This is for making the INLA mesh
Method = c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km = 25
n_x = 100   # Specify number of stations (a.k.a. "knots")

#These are the configuration for the VAST model
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1) 
RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0) 
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
ObsModel = c(2,0)  

#These are the derived variables
Options =  c("SD_site_density"=FALSE, "SD_site_logdensity"=FALSE, "Calculate_Range"=TRUE, "Calculate_effective_area"=TRUE)



# Default
if( Data_Set %in% c("GSL_american_plaice","BC_pacific_cod","EBS_pollock","SAWC_jacopever","Chatham_rise_hake","Aleutian_islands_POP")){
  strata.limits <- data.frame('STRATA'="All_areas")
}
# Specific (useful as examples)
if( Data_Set %in% c("WCGBTS_canary","Sim")){
  # In this case, it will calculate a coastwide index, and also a separate index for each state (although the state lines are approximate)
  strata.limits <- data.frame(
    'STRATA' = c("Coastwide","CA","OR","WA"),
    'north_border' = c(49.0, 42.0, 46.0, 49.0),
    'south_border' = c(32.0, 32.0, 42.0, 46.0),
    'shallow_border' = c(55, 55, 55, 55),
    'deep_border' = c(1280, 1280, 1280, 1280)
  )
  # Override default settings for vessels
  OverdispersionConfig = c("Delta1"=1, "Delta2"=1)
}
if( Data_Set %in% c("GOA_Pcod","GOA_pollock")){
  # In this case, will calculating an unrestricted index and a separate index restricted to west of -140W
  strata.limits <- data.frame(
    'STRATA' = c("All_areas", "west_of_140W"),
    'west_border' = c(-Inf, -Inf),
    'east_border' = c(Inf, -140)
  )
}
if( Data_Set %in% c("GB_spring_haddock","GB_fall_haddock")){
  # For NEFSC indices, strata must be specified as a named list of area codes
  strata.limits = list( 'Georges_Bank'=c(1130, 1140, 1150, 1160, 1170, 1180, 1190, 1200, 1210, 1220, 1230, 1240, 1250, 1290, 1300) )
}
if( Data_Set %in% c("Iceland_cod")){
  strata.limits = data.frame( 'STRATA'="All_areas" )
  # Turn off all spatial, temporal, and spatio-temporal variation in probability of occurrence, because they occur almost everywhere
  FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=1, "Epsilon2"=1)
  # Use an observation model that fixes encounter probability at 100% for years where its encoutered everywhere.
  ObsModel = c(2,3)
}

Region = switch( Data_Set, 
                 "Chatham_rise_hake"="New_Zealand", 
                 "WCGBTS_canary"="California_current", 
                 "GSL_american_plaice"="Gulf_of_St_Lawrence", 
                 "BC_pacific_cod"="British_Columbia", 
                 "EBS_pollock"="Eastern_Bering_Sea", 
                 "GOA_Pcod"="Gulf_of_Alaska", 
                 "GOA_pollock"="Gulf_of_Alaska", 
                 "GB_spring_haddock"="Northwest_Atlantic", 
                 "GB_fall_haddock"="Northwest_Atlantic", 
                 "SAWC_jacopever"="South_Africa", 
                 "Aleutian_islands_POP"="Aleutian_Islands",
                 "Other")

DateFile = paste0(getwd(),'/VAST_output/')
dir.create(DateFile)

Version = get_latest_version( package="VAST" )
Record = ThorsonUtilities::bundlelist( c("Data_Set","Version","Method","grid_size_km","n_x","FieldConfig","RhoConfig","OverdispersionConfig","ObsModel","Options") )
save( Record, file=file.path(DateFile,"Record.RData"))
capture.output( Record, file=paste0(DateFile,"Record.txt"))

Method = c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km = 25
n_x = 1000   # Specify number of stations (a.k.a. "knots")


FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1) 
RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0) 
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
ObsModel = c(2,0)  

Options =  c("SD_site_density"=FALSE, "SD_site_logdensity"=FALSE, "Calculate_Range"=TRUE, "Calculate_effective_area"=TRUE)

# Default
if( Data_Set %in% c("GSL_american_plaice","BC_pacific_cod","EBS_pollock","SAWC_jacopever","Chatham_rise_hake","Aleutian_islands_POP")){
  strata.limits <- data.frame('STRATA'="All_areas")
}
# Specific (useful as examples)
if( Data_Set %in% c("WCGBTS_canary","Sim")){
  # In this case, it will calculate a coastwide index, and also a separate index for each state (although the state lines are approximate)
  strata.limits <- data.frame(
    'STRATA' = c("Coastwide","CA","OR","WA"),
    'north_border' = c(49.0, 42.0, 46.0, 49.0),
    'south_border' = c(32.0, 32.0, 42.0, 46.0),
    'shallow_border' = c(55, 55, 55, 55),
    'deep_border' = c(1280, 1280, 1280, 1280)
  )
  # Override default settings for vessels
  OverdispersionConfig = c("Delta1"=1, "Delta2"=1)
}
if( Data_Set %in% c("GOA_Pcod","GOA_pollock")){
  # In this case, will calculating an unrestricted index and a separate index restricted to west of -140W
  strata.limits <- data.frame(
    'STRATA' = c("All_areas", "west_of_140W"),
    'west_border' = c(-Inf, -Inf),
    'east_border' = c(Inf, -140)
  )
}
if( Data_Set %in% c("GB_spring_haddock","GB_fall_haddock")){
  # For NEFSC indices, strata must be specified as a named list of area codes
  strata.limits = list( 'Georges_Bank'=c(1130, 1140, 1150, 1160, 1170, 1180, 1190, 1200, 1210, 1220, 1230, 1240, 1250, 1290, 1300) )
}
if( Data_Set %in% c("Iceland_cod")){
  strata.limits = data.frame( 'STRATA'="All_areas" )
  # Turn off all spatial, temporal, and spatio-temporal variation in probability of occurrence, because they occur almost everywhere
  FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=1, "Epsilon2"=1)
  # Use an observation model that fixes encounter probability at 100% for years where its encoutered everywhere.
  ObsModel = c(2,3)
}


Region = switch( Data_Set, "Chatham_rise_hake"="New_Zealand", 
                 "WCGBTS_canary"="California_current", 
                 "GSL_american_plaice"="Gulf_of_St_Lawrence", 
                 "BC_pacific_cod"="British_Columbia", 
                 "EBS_pollock"="Eastern_Bering_Sea", 
                 "GOA_Pcod"="Gulf_of_Alaska", 
                 "GOA_pollock"="Gulf_of_Alaska", 
                 "GB_spring_haddock"="Northwest_Atlantic", 
                 "GB_fall_haddock"="Northwest_Atlantic", 
                 "SAWC_jacopever"="South_Africa", 
                 "Aleutian_islands_POP"="Aleutian_Islands",
                 "Other")

if(Data_Set=="WCGBTS_canary"){
  data( WCGBTS_Canary_example, package="FishStatsUtils" )
  Year = as.numeric(sapply(WCGBTS_Canary_example[,'PROJECT_CYCLE'], FUN=function(Char){strsplit(as.character(Char)," ")[[1]][2]}))
  Data_Geostat = data.frame( "Catch_KG"=WCGBTS_Canary_example[,'HAUL_WT_KG'], "Year"=Year, "Vessel"=paste(WCGBTS_Canary_example[,"VESSEL"],Year,sep="_"), "AreaSwept_km2"=WCGBTS_Canary_example[,"AREA_SWEPT_HA"]/1e2, "Lat"=WCGBTS_Canary_example[,'BEST_LAT_DD'], "Lon"=WCGBTS_Canary_example[,'BEST_LON_DD'], "Pass"=WCGBTS_Canary_example[,'PASS']-1.5)
}
if( Data_Set %in% c("BC_pacific_cod")){
  data( BC_pacific_cod_example, package="FishStatsUtils" )
  Data_Geostat = data.frame( "Catch_KG"=BC_pacific_cod_example[,'PCOD_WEIGHT'], "Year"=BC_pacific_cod_example[,'Year'], "Vessel"="missing", "AreaSwept_km2"=BC_pacific_cod_example[,'TOW.LENGTH..KM.']/100, "Lat"=BC_pacific_cod_example[,'LAT'], "Lon"=BC_pacific_cod_example[,'LON'], "Pass"=0)
}
if( Data_Set %in% c("GSL_american_plaice")){
  data( GSL_american_plaice, package="FishStatsUtils" )
  Print_Message( "GSL_american_plaice" )
  Data_Geostat = data.frame( "Year"=GSL_american_plaice[,'year'], "Lat"=GSL_american_plaice[,'latitude'], "Lon"=GSL_american_plaice[,'longitude'], "Vessel"="missing", "AreaSwept_km2"=GSL_american_plaice[,'swept'], "Catch_KG"=GSL_american_plaice[,'biomass']*GSL_american_plaice[,'vstd'] )
}
if(Data_Set=="EBS_pollock"){
  data( EBS_pollock_data, package="FishStatsUtils" )
  Data_Geostat = data.frame( "Catch_KG"=EBS_pollock_data[,'catch'], "Year"=EBS_pollock_data[,'year'], "Vessel"="missing", "AreaSwept_km2"=0.01, "Lat"=EBS_pollock_data[,'lat'], "Lon"=EBS_pollock_data[,'long'], "Pass"=0)
}
if(Data_Set=="GOA_Pcod"){
  data( GOA_pacific_cod , package="FishStatsUtils")
  Data_Geostat = data.frame( "Catch_KG"=GOA_pacific_cod[,'catch'], "Year"=GOA_pacific_cod[,'year'], "Vessel"="missing", "AreaSwept_km2"=0.01, "Lat"=GOA_pacific_cod[,'lat'], "Lon"=GOA_pacific_cod[,'lon'], "Pass"=0)
}
if(Data_Set=="GOA_pollock"){
  data( GOA_walleye_pollock, package="FishStatsUtils" )
  Data_Geostat = data.frame( "Catch_KG"=GOA_walleye_pollock[,'catch'], "Year"=GOA_walleye_pollock[,'year'], "Vessel"="missing", "AreaSwept_km2"=0.01, "Lat"=GOA_walleye_pollock[,'lat'], "Lon"=GOA_walleye_pollock[,'lon'], "Pass"=0)
}
if(Data_Set=="Aleutian_islands_POP"){
  data( AI_pacific_ocean_perch, package="FishStatsUtils" )
  Data_Geostat = data.frame( "Catch_KG"=AI_pacific_ocean_perch[,'cpue..kg.km.2.'], "Year"=AI_pacific_ocean_perch[,'year'], "Vessel"="missing", "AreaSwept_km2"=1, "Lat"=AI_pacific_ocean_perch[,'start.latitude'], "Lon"=AI_pacific_ocean_perch[,'start.longitude'], "Pass"=0)
}
if( Data_Set=="GB_spring_haddock"){
  data( georges_bank_haddock_spring, package="FishStatsUtils" )         
  Print_Message( "GB_haddock" )
  Data_Geostat = data.frame( "Catch_KG"=georges_bank_haddock_spring[,'CATCH_WT_CAL'], "Year"=georges_bank_haddock_spring[,'YEAR'], "Vessel"="missing", "AreaSwept_km2"=0.0112*1.852^2, "Lat"=georges_bank_haddock_spring[,'LATITUDE'], "Lon"=georges_bank_haddock_spring[,'LONGITUDE'])
}
if( Data_Set=="GB_fall_haddock"){
  data( georges_bank_haddock_fall, package="FishStatsUtils" )         
  Print_Message( "GB_haddock" )
  Data_Geostat = data.frame( "Catch_KG"=georges_bank_haddock_fall[,'CATCH_WT_CAL'], "Year"=georges_bank_haddock_fall[,'YEAR'], "Vessel"="missing", "AreaSwept_km2"=0.0112*1.852^2, "Lat"=georges_bank_haddock_fall[,'LATITUDE'], "Lon"=georges_bank_haddock_fall[,'LONGITUDE'])
}
if( Data_Set=="SAWC_jacopever"){
  data( south_africa_westcoast_jacopever, package="FishStatsUtils" )         
  Data_Geostat = data.frame( "Catch_KG"=south_africa_westcoast_jacopever[,'HELDAC'], "Year"=south_africa_westcoast_jacopever[,'Year'], "Vessel"="missing", "AreaSwept_km2"=south_africa_westcoast_jacopever[,'area_swept_nm2']*1.852^2, "Lat"=south_africa_westcoast_jacopever[,'cen_lat'], "Lon"=south_africa_westcoast_jacopever[,'cen_long'])
}
if( Data_Set %in% c("Iceland_cod")){
  # WARNING:  This data set has not undergone much evaluation for spatio-temporal analysis
  data( iceland_cod, package="FishStatsUtils" )
  Data_Geostat = data.frame( "Catch_KG"=iceland_cod[,'Catch_b'], "Year"=iceland_cod[,'year'], "Vessel"=1, "AreaSwept_km2"=iceland_cod[,'towlength'], "Lat"=iceland_cod[,'lat1'], "Lon"=iceland_cod[,'lon1'])
}
if( Data_Set %in% c("Chatham_rise_hake")){
  data( chatham_rise_hake, package="FishStatsUtils" )
  Data_Geostat = data.frame( "Catch_KG"=chatham_rise_hake[,'Hake_kg_per_km2'], "Year"=chatham_rise_hake[,'Year'], "Vessel"=1, "AreaSwept_km2"=1, "Lat"=chatham_rise_hake[,'Lat'], "Lon"=chatham_rise_hake[,'Lon'])
}
Data_Geostat = na.omit( Data_Geostat )


if( Region %in% c("California_current","Eastern_Bering_Sea","Gulf_of_Alaska","Aleutian_Islands","Northwest_Atlantic","Gulf_of_St_Lawrence","New_Zealand") ){
  Extrapolation_List = make_extrapolation_info( Region=Region, strata.limits=strata.limits )
}
if( Region == "British_Columbia" ){
  Extrapolation_List = make_extrapolation_info( Region=Region, strata.limits=strata.limits, strata_to_use=c("HS","QCS") )
}
if( Region == "South_Africa" ){
  Extrapolation_List = make_extrapolation_info( Region=Region, strata.limits=strata.limits, region="west_coast" )
}
if( Region == "Other" ){
  Extrapolation_List = make_extrapolation_info( Region=Region, strata.limits=strata.limits, observations_LL=Data_Geostat[,c('Lat','Lon')], maximum_distance_from_sample=15 )
}

Spatial_List = make_spatial_info( grid_size_km=grid_size_km, n_x=n_x, Method=Method, Lon_i=Data_Geostat[,'Lon'], Lat_i=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, DirPath=DateFile, Save_Results=FALSE, fine_scale=TRUE )
# Add knots to Data_Geostat
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

TmbData = make_data("Version"=Version, "FieldConfig"=FieldConfig, "OverdispersionConfig"=OverdispersionConfig, "RhoConfig"=RhoConfig, "ObsModel"=ObsModel, "c_i"=rep(0,nrow(Data_Geostat)), "b_i"=Data_Geostat[,'Catch_KG'], "a_i"=Data_Geostat[,'AreaSwept_km2'], "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1, "s_i"=Data_Geostat[,'knot_i']-1, "t_i"=Data_Geostat[,'Year'], "spatial_list"=Spatial_List, "Options"=Options )

TmbList = make_model("TmbData"=TmbData, "RunDir"=DateFile, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Spatial_List$Method)
Obj = TmbList[["Obj"]]

Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=DateFile, bias.correct=TRUE, newtonsteps=1, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

