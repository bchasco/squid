#Plot everything for all years
plot_results(fit=fit, plot_set = c(1,2,3))

#Figure 2 & 3
plot_results(fit=fit, plot_set = c(2,3), years_to_plot = c(1,7,14,21))

#Figure 4
plot_biomass_index(fit$data_list,fit$parameter_estimates$SD)
df <- data.frame(Idx=c(log(fit$Report$Index_cyl)),
                 yr=rep(fit$year_labels,nrow(strata.limits)),
                 strata=rep(strata.limits$STRATA,each=length(fit$year_labels)),
                 sd = fit$parameter_estimates$SD$sd["ln_Index_cyl"==names(fit$parameter_estimates$SD$value)])

library(ggplot2)
p <- ggplot(data=df,aes(y=exp(Idx)/1000,x=yr,colour=strata))
p <- p + geom_line(size=0.5)
# p <- p + geom_point()
p <- p + xlab("Year") + ylab("Index of abundance (1000)")
p <- p + geom_ribbon(aes(ymin=exp(Idx-0.67*sd)/1000,ymax=exp(Idx+0.67*sd)/1000, fill=strata), alpha=0.1)

print(p)

#Figure 5
source("plot_density_dependence.r")

#Figure 6
source('Cog_plot_MLE.r')


#Covariate data
#encounter, density
beta <- fit$parameter_estimates$par[grep("Beta", names(fit$parameter_estimates$par))]
#Temp is the first row
lambda <- matrix(fit$parameter_estimates$par[grep("lambda", names(fit$parameter_estimates$par))],2,2)
mu <- c(plogis(beta[1]),exp(beta[2]))

rawQ_ik <- raw[,c('Top20m_Temp','Top20m_Salinity')] #rep(1,nrow(raw))
apply(rawQ_ik,c(2),function(x){return(quantile(na.omit(x),probs = c(0.025,0.975)))})

meanQ_ik <- apply(na.omit(rawQ_ik),c(2),mean)
sdQ_ik <- apply(na.omit(rawQ_ik),c(2),sd)

plusTef <- c(plogis(beta[1]+lambda[1,1]/sdQ_ik[1])/plogis(beta[1]),exp(beta[2]+lambda[1,2]/sdQ_ik[1])/exp(beta[2]))
plusSef <- c(plogis(beta[1]+lambda[2,1]/sdQ_ik[2])/plogis(beta[1]),exp(beta[2]+lambda[2,2]/sdQ_ik[2])/exp(beta[2]))

#Index results
(fit$Report$Index_cyl[1,dim(fit$Report$Index_cyl)[2],]-fit$Report$Index_cyl[1,1,])/fit$Report$Index_cyl[1,1,]

