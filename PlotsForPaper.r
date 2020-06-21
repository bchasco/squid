#Figure 1
df_data <- data.frame(E_km = fit$spatial_list$MeshList$loc_x[fit$spatial_list$knot_i,1],
                      N_km = fit$spatial_list$MeshList$loc_x[fit$spatial_list$knot_i,2],
                      Year = raw$Year,
                      Lon = fit$spatial_list$latlon_i[,2],
                      Lat = fit$spatial_list$latlon_i[,1],
                      knot_i = fit$spatial_list$knot_i)

plot_data(Extrapolation_List = fit$extrapolation_list,
          Spatial_List = fit$spatial_list,
          Data_Geostat = df_data,
          PlotDir = paste0(getwd(),"/VAST_plots/"),
          col = as.numeric(as.factor(raw$Survey)),
          cex=0.7, 
          las=1)

#Figures 2 
years_to_plot <- seq(1,22,2)
plot_results(fit, 
             working_dir = paste0(getwd(),"/VAST_plots/"),
             # years_to_plot = years_to_plot,
             plot_set=1)

#Figures 3, 4, 5 
plot_results(fit, 
             working_dir = paste0(getwd(),"/VAST_plots/"),
             strata_names = order(strata.limits$STRATA),
             # years_to_plot = years_to_plot,
             plot_set=2)


#Figures 4 
plot_results(fit, 
             working_dir = paste0(getwd(),"/VAST_plots/"),
             # years_to_plot = years_to_plot,
             strata_names = strata.limits$STRATA, #You have to sort these alphabetically to get the plots to make sense
             plot_set=7)

#Figure 5

mySD <- matrix(fit$parameter_estimates$SD$value[names(fit$parameter_estimates$SD$value)=="effective_area_cyl"],22,4)
myIndex <- fit$Report$effective_area_cyl[1,,]/1000
upr <- myIndex+0.674*mySD/1000
lwr <- myIndex-0.674*mySD/1000
DF <- data.frame(strata = rep(strata.limits$STRATA,each=22), 
                 year=rep(1998:2019,4), 
                 index=c(myIndex),
                 upr = c(upr), 
                 lwr = c(lwr))
p1 <- ggplot(DF, aes(x=year, y = index, fill=strata, color=strata)) + 
  geom_line(size=2)+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.2, colour = NA) +
  ylab("Effective area occupied (1000 km2)") +
  xlab("Calendar year")# +
  # ylim(0,max(upr))
print(p1)

#Figure 6
mySD <- matrix(fit$parameter_estimates$SD$value[names(fit$parameter_estimates$SD$value)=="Index_cyl"],22,4)
myIndex <- fit$Report$Index_cyl[1,,]
upr <- myIndex+0.674*mySD
lwr <- myIndex-0.674*mySD
DF <- data.frame(strata = rep(strata.limits$STRATA,each=22), 
                 year=rep(1998:2019,4), 
                 index=c(myIndex),
                 upr = c(upr), 
                 lwr = c(lwr),
                 percent=c(t(t(myIndex)/myIndex[1,])),
                 upr_per=c(t(t(myIndex)/upr[1,])), 
                 lwr_per = c(t(t(myIndex)/lwr[1,])))
p1 <- ggplot(DF, aes(x=year, y = index, fill=strata, color=strata)) + 
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.2, colour = NA) +
  ylab("Index (squid/km2)") +
  xlab("Calendar year") +
  ylim(0,max(upr))
p2 <- ggplot(DF, aes(x=year, y = percent, fill=strata, color=strata)) + 
  geom_ribbon(aes(ymin=lwr_per, ymax=upr_per), alpha=0.2) +
  ylab("Ratio of annual index relative to 1998") +
  xlab("Calendar year")
plot_grid(p1, p2, labels = c('A', 'B'), nrow=2, label_size = 12)
print(p)

#Supplemental figure of size distributions
lens <- read.csv("comb_lengths.csv", header=TRUE)
p <- ggplot(lens, aes(x=Length, fill = Survey, color=Survey)) + 
  geom_histogram(alpha=0.1) + 
  facet_wrap(~ Year)
print(p)
