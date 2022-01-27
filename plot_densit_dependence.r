png("plot_density_dependence.png", units="in", width=5, height=8, res=600)
load("AllMLE.Rdata")
fit <- All_fit
myval <- fit$parameter_estimates$SD$value
mysd <- fit$parameter_estimates$SD$sd
ef <- myval[names(myval)=="log_effective_area_ctl"]
ef_sd <- mysd[names(myval)=="log_effective_area_ctl"]

par(mfrow=c(3,1), mai=c(0.8,0.9,0.1,0.9))
plot(1998:2019,
     exp(ef[1:22]),
     ylab="Effective area occupied (km^2)",
     xlab="Year",
     pch=16,
     cex=2,
     col="grey",
     cex.lab=1.5,
     cex.axis=1.5,
     ylim=c(min(exp(ef[1:22]-0.67*ef_sd[1:22])),max(exp(ef[1:22]+0.67*ef_sd[1:22])))
)
text(2018,max(exp(ef[1:22]+0.67*ef_sd[1:22]))*0.9,"(a)")
for(i in 1998:2019){
  segments(i,exp(ef[i-1997]),i,exp(ef[i-1997]+0.67*ef_sd[i-1997]), lty=2, lwd=0.5)
  segments(i,exp(ef[i-1997]),i,exp(ef[i-1997]-0.67*ef_sd[i-1997]), lty=2, lwd=0.5)
}
points(1998:2019,
       exp(ef[1:22]),
       pch=16,
       cex=2,
       col="grey")

idx <- myval[names(myval)=="ln_Index_ctl"]
idx_sd <- mysd[names(myval)=="ln_Index_ctl"]

plot(1998:2019,
     exp(idx[1:22])/1000,
     ylab="Abundance index (1000)",
     xlab="Year",
     pch=16,
     cex=2,
     col="grey",
     cex.lab=1.5,
     cex.axis=1.5,
     ylim=exp(c(min(idx[1:22]-0.67*idx_sd[1:22]),max(idx[1:22]+0.67*idx_sd[1:22])))/1000
)
text(2018,max(exp(idx[1:22]+0.67*idx_sd[1:22]))/1000*0.9,"(b)")
for(i in 1998:2019){
  segments(i,exp(idx[i-1997])/1000,i,exp(idx[i-1997]+0.67*idx_sd[i-1997])/1000, lty=2, lwd=0.5)
  segments(i,exp(idx[i-1997])/1000,i,exp(idx[i-1997]-0.67*idx_sd[i-1997])/1000, lty=2, lwd=0.5)
}
points(1998:2019,
       exp(idx[1:22])/1000,
       pch=16,
       cex=2,
       col="grey")


plot((idx[1:22]),
     (ef[1:22]),
     xlab="log(Abundance index)",
     ylab="log(EAO (km^2))",
     pch=16,
     cex=2,
     col="grey",
     cex.lab=1.5,
     cex.axis=1.5,
     xlim=(c(min(idx[1:22]-0.67*idx_sd[1:22]),max(idx[1:22]+0.67*idx_sd[1:22]))),
     ylim=c(min((ef[1:22]-0.67*ef_sd[1:22])),max((ef[1:22]+0.67*ef_sd[1:22])))
)
text(max(idx[1:22]+0.67*idx_sd[1:22])*0.99,max((ef[1:22]+0.67*ef_sd[1:22]))*0.99,"(c)")
#Error idx
for(i in 1998:2019){
  segments((idx[i-1997]),(ef[i-1997]),(idx[i-1997]+0.67*idx_sd[i-1997]),(ef[i-1997]), lty=2, lwd=0.5)
  segments((idx[i-1997]),(ef[i-1997]),(idx[i-1997]-0.67*idx_sd[i-1997]),(ef[i-1997]), lty=2, lwd=0.5)
}
#Error ef
for(i in 1998:2019){
  segments((idx[i-1997]),(ef[i-1997]),(idx[i-1997]),(ef[i-1997]-0.67*ef_sd[i-1997]), lty=2, lwd=0.5)
  segments((idx[i-1997]),(ef[i-1997]),(idx[i-1997]),(ef[i-1997]+0.67*ef_sd[i-1997]), lty=2, lwd=0.5)
}
points((idx[1:22]),
     (ef[1:22]),
     pch=16,
     cex=2,
     col="grey")

dev.off()

summary(lm(ef~idx
   ))

save(ef,ef_sd,idx,idx_sd, file="dens.rData") #Save density dependent data for estimating slope

