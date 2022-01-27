if(!exists("All_fit")){
  load("Allfit.Rdata")
} 

if(!exists("raw")){
  raw <- read.csv("Update_Comb_Catch_wTrawlDist_flat.csv")
  
  #Adjust by wainright paper for JSOES
  raw$catch[raw$Gear=="264NRT+MMED_Down" & raw$Extension=="No"] <- raw$catch[raw$Gear=="264NRT+MMED_Down" & raw$Extension=="No"]/0.48
  raw$catch[raw$Gear=="264NRT+MMED_Up" & raw$Extension=="No"] <- raw$catch[raw$Gear=="264NRT+MMED_Up" & raw$Extension=="No"]/0.88
  
  #Adjustments for SWFSC - REad Cheryl's email from 9/15/2020
  raw$catch[raw$Gear=="264NRT+MMED modified" & raw$Extension=="No"] <- raw$catch[raw$Gear=="264NRT+MMED modified" & raw$Extension=="No"]/0.48
  raw$catch[raw$Gear=="264NRT+MMED" & raw$Extension=="No"] <- raw$catch[raw$Gear=="264NRT+MMED" & raw$Extension=="No"]/0.88
  
  
  #Get rid of any blanks
  raw <- raw[!apply(raw,1,function(x)return(sum(is.na(x)))),]
  
  #9) Catchability associated with surveys
  Q_ik <- raw[,c('Top20m_Temp','Top20m_Salinity')] #rep(1,nrow(raw))

}
parNames <- names(All_fit$parameter_estimates$par)
p_mu <- All_fit$parameter_estimates$par["Beta_mean1_c"]
r_mu <- All_fit$parameter_estimates$par["Beta_mean2_c"]

p_lam <- All_fit$parameter_estimates$par[grep("lambda1", parNames)]
r_lam <- All_fit$parameter_estimates$par[grep("lambda2", parNames)]

p1 <- (plogis(p_mu+p_lam[1]*seq(-2,2,0.1))-plogis(p_mu))/plogis(p_mu)*100
r1 <- (exp(r_mu+r_lam[1]*seq(-2,2,0.1))-exp(r_mu))/exp(r_mu)*100

p2 <- (plogis(p_mu+p_lam[2]*seq(-2,2,0.1))-plogis(p_mu))/plogis(p_mu)*100
r2 <- (exp(r_mu+r_lam[2]*seq(-2,2,0.1))-exp(r_mu))/exp(r_mu)*100

par(mfrow=c(2,1))

matplot(seq(-2,2,0.1)*sd(na.omit(Q_ik[,1]))+mean(na.omit(Q_ik[,1])),
        cbind(p1,r1), type="l", las=1,
        xlab="Change in temperature",
        ylab="Percent change", col=c("orange","black"), lwd=3,lty=1)
legend(8,170,legend=c("Encounter probability","Density (#/km2)"), lwd=3, bty="n", col=c("orange","black"))
matplot(seq(-2,2,0.1)*sd(na.omit(Q_ik[,2]))+mean(na.omit(Q_ik[,2])),
        cbind(p2,r2), type="l", las=1,
        xlab="Change in salinity",
        ylab="Percent change",col=c("orange","black"), lwd=3, lty=1)
legend(30,150,legend=c("Encounter probability","Density (#/km2)"), lwd=3, bty="n", col=c("orange","black"))
