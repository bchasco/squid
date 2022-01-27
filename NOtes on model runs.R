Notes

Model runs:
  
  -Changing the spatial to zero for enc and pos catches
    +This led to convergence of all models; however, the shift in the COG seems wonky
      =This is likely due to the fact that rho_config was c(3,3,4,4)
      =Try changing rho to c(1,1,1,1) with field still equal to c(0,1,0,1)
        *The c(1,1,1,1) did not work for the SWFSC only data
      =Try rho c(0,0,1,1) and see if that works
        *This does seem to work with field c(0,1,0,1). You could use this to compare the 
        *The only problem is the COG doesn't seem to make sense. THe shift is doesn't show
          any impact of the heat wave.
      =Try rho c(3,3,1,1) and field c(1,1,1,1)
        *Does not work for the JSOES data
      =Try rho c(3,3,1,1) and field c(0,1,0,1)
        *This works for the JSOES data
      =Now add a flag, different JSOES (rho c(3,3,1,1) and field c(0,1,0,1)) and SWFSC and All (rho c(3,3,4,4) and field c(1,1,1,1))
          


With regard to Q_ik
  -Try a COG plot with model inputs with and without Q_ik. Uses the settings that work for all JSOES, SWFSC and All


With regard to the net effect
  -Need to look and see if we can estimate a net effect based on the coupl eyears that they over lap. COG plots might imply that the SWFSC is better at catching squid.

2012 had the highest catches off SF Bay, thus we would expect this year to have the highest index
