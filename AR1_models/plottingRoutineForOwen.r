#Load some report object from TMB
load("TMBout.rData")

#Everything below here is for plotting
#Load the libraries you need
library(INLA)
library(INLA); library(sp); library(fields)
library(geoR)
library(viridisLite)
library(TMB)
library(ggthemes)
library(maps)
require(sp)
require(maptools)
require(ggplot2)

n <- 200 #Raster cells
xlim <- c(-125.5,-123.5) #longitude limits
ylim <- c(44,48.5) #latitude limits
#Projected mesh based on raster cells
proj = inla.mesh.projector(mesh, 
                           xlim = xlim,
                           ylim = ylim, 
                           dims=c(n, n))

#Project the deviates for the inla mesh
field.proj = inla.mesh.project(proj, rep$Omega_x+rep$Epsilon_xt[,1])

#Create the data.frame for doing the plotting 
DF <- data.frame(yr=rep(unique(survey$Year)[1],n*n),
                 x=rep(proj$x,n),
                 y=rep(proj$y,each=n),
                 dens = exp(c(as.matrix(field.proj))))
for(i in 2:ncol(rep$Epsilon_xt)){
  field.proj = inla.mesh.project(proj, rep$Omega_x+rep$Epsilon_xt[,i])
  DF_tmp <- data.frame(yr=rep(unique(survey$Year)[i],n*n),
                       x=rep(proj$x,n),
                       y=rep(proj$y,each=n),
                       dens = exp(c(as.matrix(field.proj))))
  DF <- rbind(DF,DF_tmp)
}
DF$yr <- as.factor(DF$yr)


# get the map, don't plot it but `fill` gets us named 'map' polygons
world <- map("world", fill=TRUE, plot=FALSE)

# convert the 'map' to something we can work with via geom_map
IDs <- sapply(strsplit(world$names, ":"), function(x) x[1])
world <- map2SpatialPolygons(world, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

# this does the magic for geom_map
world_map <- fortify(world)

#You can subset the years of the inla output
myYears <- DF$yr%in%c(as.character(1998:2019))

#This is the ggplot
p <- ggplot(DF[myYears,], aes(x, y)) +
  geom_raster(aes(x, y, fill = dens)) +
  scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
  facet_wrap(~yr, ncol=3) +
  scale_x_continuous(limits = c(min(DF$x),max(DF$x)), expand = c(0, 0)) +
  scale_y_continuous(limits = c(min(DF$y),max(DF$y)), expand = c(0, 0)) +
  theme_bw() +
  ylab("Latitude") +
  xlab("Longitude") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  annotation_map(map_data("world")) +
  theme_bw()+theme(plot.background = element_blank(),panel.grid.major = element_blank()
                   ,panel.grid.minor = element_blank()) +labs(title="",y="", x="")+
  theme(axis.title.x = element_text(face="bold", size=8,colour = rgb(0,0,0)),axis.text.x = element_text(size=8,colour = rgb(0,0,0)))+
  theme(axis.title.y = element_text(face="bold", size=8,colour = rgb(0,0,0)),axis.text.y = element_text(size=8,colour = rgb(0,0,0)))+
  theme(plot.title = element_text(lineheight=.8, face="bold",size=10,colour = rgb(0,0,0)))+
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  theme(legend.position="none")
#Add the map as a base layer before the points

# pdf("SpatialAnalysis_INLA.pdf", width = 7, height=10)
print(p)
# ggsave(plot = p, file = "SpatialAnalysis_INLA.png",
# type = "cairo-png",  bg = "transparent",
# width = 8, height = 8, units = "cm", dpi = 800)

# # dev.off()
# 
