# library(ggplot2)
# library(ggmap)
# ggimage((rep$eps_c_yl), fullpage = TRUE)
theme_set(theme_minimal())

df <- data.frame(y=unique(Data$t_i+1997),
                 c_y=SD$value[names(SD$value)=="eps_c_y"],
                 c_y_sd=SD$sd[names(SD$value)=="eps_c_y"],
                 p_y=SD$value[names(SD$value)=="eps_p_y"],
                 p_y_sd=SD$sd[names(SD$value)=="eps_p_y"])

p_c_y <- ggplot(df, aes(x=y,y=exp(c_y+rep$beta_c[1])/exp(c_y[1]+rep$beta_c[1]))) +
  geom_ribbon(aes(ymin = exp(rep$beta_c[1]+c_y - 0.675*c_y_sd)/exp(rep$beta_c[1]+c_y[1]),
                  ymax = exp(rep$beta_c[1] + c_y + 0.675*c_y_sd)/exp(rep$beta_c[1]+c_y[1])),
              fill = "grey90") +
  ylab("Perecent change CPUE") +
  # scale_y_log10()+
  # ylim(0,600) +
  xlab("Calendar year") +
  geom_line(color = "firebrick", size = 1)
print(p_c_y)

theme_set(theme_minimal())

p_p_y <- ggplot(df, aes(x=y,y=plogis(p_y+rep$beta_p[1]))) +
  geom_ribbon(aes(ymin = plogis(rep$beta_p[1]+p_y - 0.675*p_y_sd), 
                  ymax = plogis(rep$beta_p[1] + p_y + 0.675*p_y_sd)), 
              fill = "grey90") +
  ylab("CPUE = 0") +
  ylim(0,1) +
  xlab("Calendar year") +
  geom_line(color = "firebrick", size = 1) 
print(p_p_y)

df <- data.frame(lat=lati_lu$lat,
                 c_lat=SD$value[names(SD$value)=="eps_c_lat"],
                 c_lat_sd=SD$sd[names(SD$value)=="eps_c_lat"],
                 p_lat=SD$value[names(SD$value)=="eps_p_lat"],
                 p_lat_sd=SD$sd[names(SD$value)=="eps_p_lat"])

p_p_lat <- ggplot(df, aes(y=plogis(p_lat+rep$beta_p[1]),x=lat)) +
  geom_ribbon(aes(ymin = plogis(rep$beta_p[1] + p_lat - 0.675*p_lat_sd),
                  ymax = plogis(rep$beta_p[1] + p_lat + 0.675*p_lat_sd)),
              fill = "grey90") +
  ylab("CPUE = 0") +
  ylim(0,1)+
  xlab("Latitude") +
  geom_line(color = "firebrick", size = 1) +
  coord_flip()
print(p_p_lat)


p_c_lat <- ggplot(df, aes(x=lat,y=exp(c_lat+rep$beta_c[1])/exp(c_lat[1]+rep$beta_c[1])*100-100)) +
  geom_ribbon(aes(ymin = exp(rep$beta_c[1] + c_lat - 0.675*c_lat_sd)/exp(c_lat[1]+rep$beta_c[1])*100-100,
                  ymax = exp(rep$beta_c[1] + c_lat + 0.675*c_lat_sd)/exp(c_lat[1]+rep$beta_c[1])*100-100),
              fill = "grey90") +
  ylab("Perecent change CPUE") +
  xlab("Latitude") +
  ylim(-100,150) +
  geom_line(color = "firebrick", size = 1) +
  coord_flip()
print(p_c_lat)

world_map <- map_data("world")
ggplot(world_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill="white", colour = "black") +
  ylim(40,50)


library(maps)
figure <- ggarrange(p_p_y, p_c_y, 
                    p_p_lat, p_c_lat,
                    ncol = 2, nrow = 2)


ggplot(data = world) +
  geom_sf() +
  # geom_point(data = sites, aes(x = longitude, y = latitude), size = 4, 
  #            shape = 23, fill = "darkred") +
  coord_sf(xlim = c( -126.4235, -122.4033), ylim = c(37.16667, 48.34717), expand = TRUE)

figure
