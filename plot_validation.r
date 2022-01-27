library(ggplot2)
x <- read.table("predNLlcomp.dat", header=TRUE)
p <- ggplot(x[,],aes(x=Year,y=Val,group=EnvironmentalData)) +
  geom_line(aes(color=EnvironmentalData),size=1.2) +
  facet_grid(vars(Category), scales = "free") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  ylab("")
  
  
print(p)