####################################################################################
##The data used in this analysis is publicly available on R from the package mdhglm#
####################################################################################
#upload of useful libraries
library(quantreg)
library(nnet)

##################################################################################
##################definition of functions to be used in the analyses##############
##################################################################################

##function to compute the correlation starting from multinomial predictions
Phi <-function(P){
  EstCov <- P[,1]-(P[,1]+P[,3])*(P[,1]+P[,4])
  var1 <- (P[,4]+P[,2])*(1-P[,4]-P[,2])
  var2 <- (P[,3]+P[,2])*(1-P[,3]-P[,2])
  EstCor <- EstCov/(sqrt(var1)*sqrt(var2))
  EstCor
}

##logistic transform
logit <- function(x,m){
  log((x-m+1/100)/(1-x+1/100))
}

##inverse logistic transform
expit <- function(x,m){
  (exp(x)*(1+1/100)+m-1/100)/(exp(x)+1)
}

###################################################################
#################dataset upload and management#####################
###################################################################
library(mdhglm)
data(nmsqt)
help(nmswt)
twin <- nmsqt

#data management to prepare the dataset for the analyses
twin.a <- twin[twin$number %% 2==0,]
twin.b <- twin[twin$number %% 2==1,]

names(twin.a)[1] <- "number.a"
names(twin.b)[1] <- "number.b" 

twin.a$y.a <- twin.a$y1+twin.a$y2+twin.a$y3+twin.a$y4
twin.a <- twin.a[,-c(7:10)]
twin.b$y.b <- twin.b$y1+twin.b$y2+twin.b$y3+twin.b$y4
twin.b <- twin.b[,-c(7:10)]
twin.t <- merge(twin.a,twin.b, by.x="pairnum", by.y="pairnum")
twin.t <- twin.t[,-c(10:14)]

twin.t$z1 <- ifelse(twin.t$x1.x==1,1,0)
twin.t$z2 <- ifelse(twin.t$x4.x>=4,1,0)
twin.t$z3 <- ifelse(twin.t$x2.x>=3 | twin.t$x3.x>=3,1,0)
twin.t$z4 <- 1-twin.t$x5.x

#final dataset
twin <- twin.t

#preliminary analysis of the data
head(twin)
summary(twin)

twin$z1 <- as.factor(twin$z1)
twin$z2 <- as.factor(twin$z2)
twin$z3 <- as.factor(twin$z3)
twin$z4 <- as.factor(twin$z4) 

attach(twin)

cor(y.a,y.b)
plot(y.a,y.b)

#zigosity group. Set grid of values to use in predictions
grid.z4 <- c(0,1)
grid.z4 <- factor(grid.z4)

##############################################
##   Model estimation section             ###
#############################################

##set quantiles 
p <- 100:400/500

##useful initializations
cm1sp <- cm2sp <- NULL
estimates_mm1 <- ci_low_mm1 <- ci_up_mm1 <- NULL
estimates_mm2 <- ci_low_mm2 <- ci_up_mm2 <- NULL
lci_pred <- uci_pred <- NULL
PhiEstlong <- NULL


for(t in p){
  #windows()
  #par(mfrow=c(2,2),las=1,mai=c(1.3,1.4,0.7,0.4))

  ##univariate quantile models
  #outcome 1
  m1 <- rq(y.a ~ z1 + z2 + z3 + z4 + z2*z4, tau = t, data=twin)
  m1s <- summary(m1,se = "boot", R = 1000)
  cm1s <- cbind(m1s$coefficients,t)
  cm1sp <- rbind(cm1sp,cm1s)
  
  #outcome 2
  m2 <- rq(y.b ~ z1 + z2 + z3 + z4 + z2*z4, tau = t, data=twin)
  m2s <- summary(m2,se = "boot", R = 1000)
  cm2s <- cbind(m2s$coefficients,t)
  cm2sp <- rbind(cm2sp,cm2s)
  
  ##Sign concordance of quantile regression residuals
  res1 <- m1$residuals
  res2 <- m2$residuals
  
  sr1 <- sign(res1)
  sr1 <- ifelse(sr1==0,1,sr1)  
  sr2 <- sign(res2)
  sr2 <- ifelse(sr2==0,1,sr2) 
  
  
  ##Z variable of concordance
  Zm <- NA 
  Zm <- ifelse(sr1==1 & sr2==1,0,Zm)
  Zm <- ifelse(sr1==-1 & sr2==-1,1,Zm) 
  Zm <- ifelse(sr1==1 & sr2==-1 | sr1==-1 & sr2==1,2,Zm)
  
  ###Multinomial logistic model
  mod <- multinom(Zm ~   z4 , data=twin)
  mm <- summary(mod)
  coef_mod <- mm$coefficients
  
  mm1 <- coef_mod[1,]
  mm2 <- coef_mod[2,]
  
  ##Set new data values for predictions
  newdata.z4 <- data.frame(z4 = grid.z4)
  
  ##predictions of multinomial model
  p_hat <- predict(mod, type = "probs", newdata=newdata.z4)
  p_hat_n <- cbind(p_hat[,1:2], p_hat[,3]/2, p_hat[,3]/2)
  

  ##computing the correlation for the entire matrix of predictions
  PhiEst <- Phi(p_hat_n)
  ##storing of predictions
  PhiEstlong <- rbind(PhiEstlong,PhiEst)
  
  ###############################################
  ##########bootstrap estimates###################
  ###############################################
  set.seed <- 1234
  n <- nrow(twin)
  Phi.m <- NULL
  boot_coef1 <- boot_coef2 <- boot_coef3 <- NULL
  
  B <- 1000 #bootstrap replicates
  j <-0
  while(j<B){
    
    w <- rexp(n) #exponential weight
    q1 <- rq(y.a ~ z1 + z2 + z3 + z4 + z2*z4, tau = t, data=twin, weights=w)
    q2 <- rq(y.b ~ z1 + z2 + z3 + z4 + z2*z4, tau = t, data=twin, weights=w)
    
    res1b <- q1$residuals
    res2b <- q2$residuals
    
    sr1b <- sign(res1b)
    sr1b <- ifelse(sr1b==0,1,sr1b)  
    sr2b <- sign(res2b)
    sr2b <- ifelse(sr2b==0,1,sr2b) 
    
    ##Variable of concordance for each replicate
    Zb <- NA 
    Zb <- ifelse(sr1b==1 & sr2b==1,0,Zb)
    Zb <- ifelse(sr1b==-1 & sr2b==-1,1,Zb) 
    Zb <- ifelse(sr1b==1 & sr2b==-1 | sr1b==-1 & sr2b==1,2,Zb)
    
    mb <- multinom(Zb ~  z4, data=twin, weights=w)
    mmb <- summary(mb)
    coef_mb <- mmb$coefficients
    mb1 <- coef_mb[1,]
    mb2 <- coef_mb[2,]
    
    
    boot_coef1 <- rbind(boot_coef1, mb1)
    boot_coef2 <- rbind(boot_coef2, mb2)
    
    p_b <- predict(mb, type="probs", newdata=newdata.z4)
    p_b_n <- cbind(p_b[,1:2], p_b[,3]/2, p_b[,3]/2) 
    
    
    phi.m <- Phi(p_b_n)
    Phi.m <- rbind(Phi.m, phi.m)
    j<- j+1
  }
  
  se_1 <- apply(boot_coef1,2,sd)
  se_2 <- apply(boot_coef2,2,sd)
  
  
  ci_low_1 <- mm1 - 1.96*se_1
  ci_low_2 <- mm2 - 1.96*se_2
  
  ci_up_1 <- mm1 + 1.96*se_1
  ci_up_2 <- mm2 + 1.96*se_2
  
  ##computation of phi_min of phi coefficient for each quantile 
  m <- NULL
  if(t <=0.5){
    m <- -t/(1-t)}else{
      m <- (t-1)/t
    }
  #confidence intervals for quantile regression estimates with outcome 1
  estimates_mm1 <- rbind(estimates_mm1, mm1)
  ci_low_mm1 <- rbind(ci_low_mm1, ci_low_1) 
  ci_up_mm1 <- rbind(ci_up_mm1, ci_up_1)
  
  #confidence intervals for quantile regression estimates with outcome 2
  estimates_mm2 <- rbind(estimates_mm2, mm2)
  ci_low_mm2 <- rbind(ci_low_mm2, ci_low_2) 
  ci_up_mm2 <- rbind(ci_up_mm2, ci_up_2)
  
  ##confidence intervals for predictions
  se_pred <- apply(Phi.m,2,sd)
  ci_low_pred <- PhiEst - 1.96*se_pred
  ci_up_pred <- PhiEst + 1.96*se_pred
  lci_pred <- rbind(lci_pred, ci_low_pred)
  uci_pred <- rbind(uci_pred, ci_up_pred)
}

####################################################
####plots in the paper  ############################
####################################################
par(mfrow=c(2,2))
par(mai=c(0.8,1.4,1.2,1))


plot(p,estimates_mm1[,1],type="l", main="Intercept", xlab=expression(tau), ylab = ' ', cex.lab=1.8,yaxt="n", lwd=2, cex.main=1.4)
axis(2, cex.axis=1, las=2)
lines(p,ci_low_mm1[,1],type="l",lty=2)
lines(p,ci_up_mm1[,1], type="l",lty=2)
abline(h=0)

mtext(expression( paste("(a) ", "log ", bgroup("(", frac(P(paste("Z = ","\"11\"")),P(paste("Z = " ,"\"00\""))), ")" ), " = ", 
      gamma[paste(paste("\"11\""):tau,",",0)]+gamma[paste(paste("\"11\""):tau,",", 1)], "Zigosity"  ) ), side = 3, line = -4, outer = TRUE)

plot(p,estimates_mm1[,2],type="l", main="Zigosity",ylim=c(min(ci_low_mm1[,2]),max(ci_up_mm1[,2])), xlab=expression(tau), ylab = ' ', 
     cex.lab=1.8,yaxt="n", lwd=2, cex.main=1.4)
axis(2, cex.axis=1, las=2)
lines(p,ci_low_mm1[,2],type="l",lty=2)
lines(p,ci_up_mm1[,2], type="l",lty=2)  
abline(h=0)

mtext(expression( paste("(b) ", "log ", bgroup("(", frac(P(paste("Z = ","\"01\"+\"10\"")),P(paste("Z = " ,"\"00\""))), ")" ), " = ", 
                        gamma[paste(paste("\"01\"+\"10\""):tau,",",0)]+gamma[paste(paste("\"01\"+\"10\""):tau,",", 1)], "Zigosity"  ) ), side = 3, line = -22, outer = TRUE)


plot(p,estimates_mm2[,1],type="l", main="Intercept",ylim=c(min(ci_low_mm2[,1]),max(ci_up_mm2[,1])), xlab=expression(tau), ylab = ' ', 
     cex.lab=1.8,yaxt="n", lwd=2, cex.main=1.4)
axis(2, cex.axis=1, las=2)
lines(p,ci_low_mm2[,1],type="l",lty=2)
lines(p,ci_up_mm2[,1], type="l",lty=2)
abline(h=0)

plot(p,estimates_mm2[,2],type="l", main="Zigosity",ylim=c(min(ci_low_mm2[,2]),max(ci_up_mm2[,2])), xlab=expression(tau), ylab = ' ', 
     cex.lab=1.8,yaxt="n", lwd=2, cex.main=1.4)
axis(2, cex.axis=1, las=2)
lines(p,ci_low_mm2[,2],type="l",lty=2)
lines(p,ci_up_mm2[,2], type="l",lty=2)
abline(h=0)


windows()
#quartz()
################predictions plot
par(mfrow=c(1,1))
yy1 <- c(lci_pred[,1], tail(uci_pred[,1], 1), rev(uci_pred[,1]), lci_pred[1,1])
yy2 <- c(lci_pred[,2], tail(uci_pred[,2], 1), rev(uci_pred[,2]), lci_pred[1,2])
xx <- c(p, tail(p, 1), rev(p), p[1])


par(mai=c(1,1.2,1,1))
plot(0,0, col = "white",
     ylab = ' ',yaxt="n", xlab=expression(tau), 
     cex.lab=1.5,
     xlim = c(0.2,0.8),
     ylim=c(0,1),
     cex.main = 1.8
)
polygon(xx, yy1, col = "grey", border = NA)
points(p, PhiEstlong[,1], type = "l", lwd = 2, lty=1)
polygon(xx, yy2, col = "grey", border = NA)
points(p, PhiEstlong[,2], type = "l", lwd = 2, lty=2)
mtext(expression(hat(phi)), side=2,line=3, cex=1.5, las=2)
axis(2, cex.axis=1, las=2)
legend("topleft", legend=c("identical", "fraternal"), lty=c(2,1),bty="n")




