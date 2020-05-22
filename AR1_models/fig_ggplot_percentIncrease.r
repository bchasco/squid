par(mfrow=c(2,2))

plot(1998:2019,
     exp(rep$beta_c[1]+rep$eps_c_y)/exp(rep$beta_c[1]+rep$eps_c_y[1])*100-100,
     type="l",las=1,
     ylab="Percent increase", 
     col="darkred", lwd=3, xlab="")

plot(1998:2019,
     exp(rep$beta_c[1]+rep$eps_c_y),
     type="l",las=1,
     ylab="Percent increase", 
     col="darkred", lwd=3, xlab="")
