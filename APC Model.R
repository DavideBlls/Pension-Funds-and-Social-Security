"setwd _ _ _"
source("Functions.r") 

### APC Model ###
#Forecast
APCfor_BELm <- forecast(APCfit_BELm, h=y.pred, gc.order = c(1, 1, 0))
APCfor_BELf <- forecast(APCfit_BELf, h=y.pred, gc.order = c(1, 1, 0))
plot(APCfor_BELm, only.kt = TRUE)
plot(APCfor_BELf, only.kt = TRUE)
plot(APCfor_BELm, only.gc = TRUE)
plot(APCfor_BELf, only.gc = TRUE)

#####
par(mfrow=c(1, 2), oma=c(2,0,0,0))
plot(c(a.min:a.max),log(APCfor_BELm$rates[,y.pred]), bty="l", xlab="Ages", ylab="Male death rates (log scale)",
  ylim=c(min(log(APCfor_BELm$rates[,y.pred])),max(log((APCfit_BELm$Dxt/APCfit_BELm$Ext)[,length(Y.fit)]))), type="l")
lines(c(a.min:a.max),log((APCfit_BELm$Dxt/APCfit_BELm$Ext)[,length(Y.fit)]), type="p")
plot(c(a.min:a.max),log(APCfor_BELf$rates[,y.pred]), bty="l", xlab="Ages", ylab="Female death rates (log scale)",
  ylim=c(min(log(APCfor_BELf$rates[,y.pred])),max(log((APCfit_BELf$Dxt/APCfit_BELf$Ext)[,length(Y.fit)]))), type="l")
lines(c(a.min:a.max),log((APCfit_BELf$Dxt/APCfit_BELf$Ext)[,length(Y.fit)]), type="p")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("last observed", "last projected"), pch=c(1,-1), lty=c(0,1), lwd=c(1,1), col=c(1:1), cex=0.9)
dev.off()

#
rates_APCfit_BELm <- fitted(APCfit_BELm, type = "rates")
ifelse(is.na(rates_APCfit_BELm),BELmRates,rates_APCfit_BELm)
rates_APCfor_BELm <- APCfor_BELm$rates
#rates_APC_BELm <- cbind(rates_APCfit_BELm,rates_APCfor_BELm)
rates_APC_BELm <- cbind(BELmRates,rates_APCfor_BELm)
q_APC_BELm <- 1- exp(-rates_APC_BELm)
q_APC_BELm.ext  <- extrapolation.fit(q_APC_BELm)

rates_APCfit_BELf <- fitted(APCfit_BELf, type = "rates")
ifelse(is.na(rates_APCfit_BELf),BELfRates,rates_APCfit_BELf)
rates_APCfor_BELf <- APCfor_BELf$rates
#rates_APC_BELf <- cbind(rates_APCfit_BELf,rates_APCfor_BELf)
rates_APC_BELf <- cbind(BELfRates,rates_APCfor_BELf)
q_APC_BELf <- 1- exp(-rates_APC_BELf)
q_APC_BELf.ext  <- extrapolation.fit(q_APC_BELf)

q_APC.BELm<-as.data.frame(q_APC_BELm.ext)
q_APC.BELf<-as.data.frame(q_APC_BELf.ext)

plot(extractCohort(q_APC_BELm.ext, cohort = 1980), type = "p", log = "y", xlab = "age", ylab = "q(x)",
 main = "Mortality rates for a cohort (log scale)")
lines(extractCohort(q_APC_BELm, cohort = 1980))

plot(extractCohort(log(q_APC_BELm.ext/(1-q_APC_BELm.ext)), cohort = 1980), type = "p", xlab = "age", ylab = "logit(q(x))",
 main = "Logit transform of mortality rates for a cohort")
lines(extractCohort(log(q_APC_BELm/(1-q_APC_BELm)), cohort = 1980))


### 1 year survival probability, n year survival probability and life expectancy
p_APC_BELm.ext <- 1-q_APC_BELm.ext # 1 year survival probability
p0n_APC_BELm.ext <- apply(p_APC_BELm.ext, 2, cumprod) # n year survival probability
ex_APC_BELm.ext <- life.exp(q_APC_BELm.ext)

p_APC_BELf.ext <- 1-q_APC_BELf.ext # 1 year survival probability
p0n_APC_BELf.ext <- apply(p_APC_BELf.ext, 2, cumprod) # n year survival probability
ex_APC_BELf.ext <- life.exp(q_APC_BELf.ext)

p_APC.BELm<-as.data.frame(p_APC_BELm.ext)
p0n_APC.BELm<-as.data.frame(p0n_APC_BELm.ext)
ex_APC.BELm<-as.data.frame(ex_APC_BELm.ext)
p_APC.BELf<-as.data.frame(p_APC_BELf.ext)
p0n_APC.BELf<-as.data.frame(p0n_APC_BELf.ext)
ex_APC.BELf<-as.data.frame(ex_APC_BELf.ext)

####
par(mfrow=c(1, 2), oma=c(2,0,0,0))
plot(c(a.min:a.max),log(APCfor_BELm$rates[,y.pred]), bty="l", xlab="Ages", ylab="Male death rates (log scale)",
     ylim=c(min(log(APCfor_BELm$rates[,y.pred])),max(log((APCfit_BELm$Dxt/APCfit_BELm$Ext)[,length(Y.fit)]))), type="l")
lines(c(a.min:a.max),log((APCfit_BELm$Dxt/APCfit_BELm$Ext)[,length(Y.fit)]), type="p")
plot(c(a.min:a.max),log(APCfor_BELf$rates[,y.pred]), bty="l", xlab="Ages", ylab="Female death rates (log scale)",
     ylim=c(min(log(APCfor_BELf$rates[,y.pred])),max(log((APCfit_BELf$Dxt/APCfit_BELf$Ext)[,length(Y.fit)]))), type="l")
lines(c(a.min:a.max),log((APCfit_BELf$Dxt/APCfit_BELf$Ext)[,length(Y.fit)]), type="p")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("last observed", "last projected"), pch=c(1,-1), lty=c(0,1), lwd=c(1,1), col=c(1:1), cex=0.9)

###
par(mfrow=c(1, 2))
plot(APCfit_BELm$years, APCfit_BELm$kt[1, ], xlim = range(APCfit_BELm$years, APCsim_BELm.mrwd$kt.s$years),
     ylim = range(APCfit_BELm$kt, APCsim_BELm.mrwd$kt.s$sim), type = "l", xlab = "year", ylab = "kt",
     main = "Lee-Carter: Period index (mrwd)")
matlines(APCsim_BELm.mrwd$kt.s$years, APCsim_BELm.mrwd$kt.s$sim[1, , ], type = "l", lty = 1)

plot(APCfit_BELm$years, (APCfit_BELm$Dxt / APCfit_BELm$Ext)["65", ], xlim = range(APCfit_BELm$years, APCsim_BELm.mrwd$years),
     ylim = range((APCfit_BELm$Dxt / APCfit_BELm$Ext)["65", ], APCsim_BELm.mrwd$rates["65", , ]), type = "l", xlab = "year", ylab = "rate",
     main = "Lee-Carter: Simulated mortality rates at age 65")
matlines(APCsim_BELm.mrwd$years, APCsim_BELm.mrwd$rates["65", , ], type = "l", lty = 1)


### SIMULATION WITH RANDOM WALK WITH DRIFT
# Male
APCsim_BELm.mrwd <- simulate(APCfit_BELm, nsim = n.sim, h=y.pred, gc.order = c(2, 2, 0))
rates_APC_BELm.st <- APCsim_BELm.mrwd$rates
q_APC_BELm.st <- 1- exp(-rates_APC_BELm.st)
q_APC_BELm.st.ext <-  extrapolation.sim(q_APC_BELm.st)

par(mfrow=c(1, 2))
plot(APCfit_BELm$years, APCfit_BELm$kt[1, ], xlim = range(APCfit_BELm$years, APCsim_BELm.mrwd$kt.s$years),
  ylim = range(APCfit_BELm$kt, APCsim_BELm.mrwd$kt.s$sim), type = "l", xlab = "year", ylab = "kt",
  main = "Lee-Carter: Period index (mrwd)")
matlines(APCsim_BELm.mrwd$kt.s$years, APCsim_BELm.mrwd$kt.s$sim[1, , ], type = "l", lty = 1)
plot(APCfit_BELm$years, (APCfit_BELm$Dxt / APCfit_BELm$Ext)["65", ], xlim = range(APCfit_BELm$years, APCsim_BELm.mrwd$years),
  ylim = range((APCfit_BELm$Dxt / APCfit_BELm$Ext)["65", ], APCsim_BELm.mrwd$rates["65", , ]), type = "l", xlab = "year", ylab = "rate",
  main = "Lee-Carter: Simulated mortality rates at age 65")
matlines(APCsim_BELm.mrwd$years, APCsim_BELm.mrwd$rates["65", , ], type = "l", lty = 1)
dev.off()

#female
APCsim_BELf.mrwd <- simulate(APCfit_BELf, nsim = n.sim, h=y.pred, gc.order = c(2, 2, 0))
rates_APC_BELf.st <- APCsim_BELf.mrwd$rates
q_APC_BELf.st <- 1- exp(-rates_APC_BELf.st)
q_APC_BELf.st.ext <-  extrapolation.sim(q_APC_BELf.st)

par(mfrow=c(1, 2))
plot(APCfit_BELf$years, APCfit_BELf$kt[1, ], xlim = range(APCfit_BELf$years, APCsim_BELf.mrwd$kt.s$years),
     ylim = range(APCfit_BELf$kt, APCsim_BELf.mrwd$kt.s$sim), type = "l", xlab = "year", ylab = "kt",
     main = "Lee-Carter: Period index (mrwd)")
matlines(APCsim_BELf.mrwd$kt.s$years, APCsim_BELf.mrwd$kt.s$sim[1, , ], type = "l", lty = 1)
plot(APCfit_BELf$years, (APCfit_BELf$Dxt / APCfit_BELf$Ext)["65", ], xlim = range(APCfit_BELf$years, APCsim_BELf.mrwd$years),
     ylim = range((APCfit_BELf$Dxt / APCfit_BELf$Ext)["65", ], APCsim_BELf.mrwd$rates["65", , ]), type = "l", xlab = "year", ylab = "rate",
     main = "Lee-Carter: Simulated mortality rates at age 65")
matlines(APCsim_BELf.mrwd$years, APCsim_BELf.mrwd$rates["65", , ], type = "l", lty = 1)

########

library("fanplot")
probs = c(2.5, 10, 25, 50, 75, 90, 97.5)
par(mfrow=c(1, 2))
matplot(APCfit_BELm$years, t(q_APC_BELm[c("65", "75", "85"),c(1:length(Y.fit)) ]),
 xlim = c(y.fit.min, (y.fit.min+y.pred)), ylim = c(0.001, 0.2), pch = 20, col = "black",
 log = "y", xlab = "year", ylab = "male mortality rate (log scale)")
fan(t(APCsim_BELm.mrwd$rates["65", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
 fan.col = colorRampPalette(c("black", "white")), ln = NULL)
fan(t(APCsim_BELm.mrwd$rates["75", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
 fan.col = colorRampPalette(c("red", "white")), ln = NULL)
fan(t(APCsim_BELm.mrwd$rates["85", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
 fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
text(y.fit.min+4, q_APC_BELm[c("65", "75", "85"), 30],
 labels = c("x = 65", "x = 75", "x = 85"))
matplot(APCfit_BELf$years, t(q_APC_BELf[c("65", "75", "85"),c(1:length(Y.fit)) ]),
 xlim = c(y.fit.min, (y.fit.min+y.pred)), ylim = c(0.001, 0.2), pch = 20, col = "black",
 log = "y", xlab = "year", ylab = "female mortality rate (log scale)")
fan(t(APCsim_BELf.mrwd$rates["65", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
 fan.col = colorRampPalette(c("black", "white")), ln = NULL)
fan(t(APCsim_BELf.mrwd$rates["75", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
 fan.col = colorRampPalette(c("red", "white")), ln = NULL)
fan(t(APCsim_BELf.mrwd$rates["85", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
 fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
text(y.fit.min+4, q_APC_BELf[c("65", "75", "85"), 30],
 labels = c("x = 65", "x = 75", "x = 85"))
dev.off()

### Bootstrap
APCboot_BELm <- bootstrap(APCfit_BELm, nBoot = n.boot, type = "semiparametric")
plot(APCboot_BELm, nCol = 3)
APCsim_BELm.boot <- simulate(APCboot_BELm, nsim = n.sim/n.boot, h = y.pred, gc.order = c(1, 1, 0))
rates_APC_BELm.boot.st <- APCsim_BELm.boot$rates
q_APC_BELm.boot.st <- 1- exp(-rates_APC_BELm.boot.st)
q_APC_BELm.boot.st.ext <-  extrapolation.sim(q_APC_BELm.boot.st)

APCboot_BELf <- bootstrap(APCfit_BELf, nBoot = n.boot, type = "semiparametric")
plot(APCboot_BELf, nCol = 3)
APCsim_BELf.boot <- simulate(APCboot_BELf, nsim = n.sim/n.boot, h = y.pred, gc.order = c(1, 1, 0))
rates_APC_BELf.boot.st <- APCsim_BELf.boot$rates
q_APC_BELf.boot.st <- 1- exp(-rates_APC_BELf.boot.st)
q_APC_BELf.boot.st.ext <-  extrapolation.sim(q_APC_BELf.boot.st)

### Confidence intervals at 90%
# 1 year death probability
conf.lev <- 0.9
q_APC_BELm.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_APC_BELm.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_APC_BELm.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_APC_BELm.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
      q_APC_BELm.q95[j,k] <- quantile(q_APC_BELm.st.ext[j,k,], probs=0.5+conf.lev/2)
      q_APC_BELm.q05[j,k] <- quantile(q_APC_BELm.st.ext[j,k,], probs=(1-conf.lev)/2)
      q_APC_BELm.boot.q95[j,k] <- quantile(q_APC_BELm.boot.st.ext[j,k,], probs=0.5+conf.lev/2)
      q_APC_BELm.boot.q05[j,k] <- quantile(q_APC_BELm.boot.st.ext[j,k,], probs=(1-conf.lev)/2)
  }
}

q_APC_BELf.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_APC_BELf.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_APC_BELf.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_APC_BELf.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
      q_APC_BELf.q95[j,k] <- quantile(q_APC_BELf.st.ext[j,k,], probs=0.5+conf.lev/2)
      q_APC_BELf.q05[j,k] <- quantile(q_APC_BELf.st.ext[j,k,], probs=(1-conf.lev)/2)
      q_APC_BELf.boot.q95[j,k] <- quantile(q_APC_BELf.boot.st.ext[j,k,], probs=0.5+conf.lev/2)
      q_APC_BELf.boot.q05[j,k] <- quantile(q_APC_BELf.boot.st.ext[j,k,], probs=(1-conf.lev)/2)
  }
}

q_APC.BELm.q95<-as.data.frame(q_APC_BELm.q95)
q_APC.BELm.q05<-as.data.frame(q_APC_BELm.q05)
q_APC.BELm.boot.q95<-as.data.frame(q_APC_BELm.boot.q95)
q_APC.BELm.boot.q05<-as.data.frame(q_APC_BELm.boot.q05)
q_APC.BELf.q95<-as.data.frame(q_APC_BELf.q95)
q_APC.BELf.q05<-as.data.frame(q_APC_BELf.q05)
q_APC.BELf.boot.q95<-as.data.frame(q_APC_BELf.boot.q95)
q_APC.BELf.boot.q05<-as.data.frame(q_APC_BELf.boot.q05)

#male
plot(q_APC_BELm.ext["65",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_APC_BELm.q95["65",], x=Y.pred, type="l", lty=2)
lines(q_APC_BELm.q05["65",], x=Y.pred, type="l", lty=2)
lines(q_APC_BELm.boot.q95["65",], x=Y.pred, type="l", lty=2, col=2)
lines(q_APC_BELm.boot.q05["65",], x=Y.pred, type="l", lty=2, col=2)

plot(q_APC_BELm.ext["20",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_APC_BELm.q95["20",], x=Y.pred, type="l", lty=2)
lines(q_APC_BELm.q05["20",], x=Y.pred, type="l", lty=2)
lines(q_APC_BELm.boot.q95["20",], x=Y.pred, type="l", lty=2, col=2)
lines(q_APC_BELm.boot.q05["20",], x=Y.pred, type="l", lty=2, col=2)

# female
plot(q_APC_BELf.ext["65",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_APC_BELf.q95["65",], x=Y.pred, type="l", lty=2)
lines(q_APC_BELf.q05["65",], x=Y.pred, type="l", lty=2)
lines(q_APC_BELf.boot.q95["65",], x=Y.pred, type="l", lty=2, col=2)
lines(q_APC_BELf.boot.q05["65",], x=Y.pred, type="l", lty=2, col=2)

plot(q_APC_BELf.ext["20",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_APC_BELf.q95["20",], x=Y.pred, type="l", lty=2)
lines(q_APC_BELf.q05["20",], x=Y.pred, type="l", lty=2)
lines(q_APC_BELf.boot.q95["20",], x=Y.pred, type="l", lty=2, col=2)
lines(q_APC_BELf.boot.q05["20",], x=Y.pred, type="l", lty=2, col=2)


##### 1 year survival probability, n year survival probability and life expectancy confidence interval
p_APC_BELm.q95 <- 1-q_APC_BELm.q95# 1 year survival probability
p0n_APC_BELm.q95 <- apply(p_APC_BELm.q95, 2, cumprod) # n year survival probability
ex_APC_BELm.q95 <- life.exp(q_APC_BELm.q95) # life expectancy

p_APC_BELm.q05 <- 1-q_APC_BELm.q05# 1 year survival probability
p0n_APC_BELm.q05 <- apply(p_APC_BELm.q05, 2, cumprod) # n year survival probability
ex_APC_BELm.q05 <- life.exp(q_APC_BELm.q05) # life expectancy

p_APC_BELm.boot.q95 <- 1-q_APC_BELm.boot.q95# 1 year survival probability
p0n_APC_BELm.boot.q95 <- apply(p_APC_BELm.boot.q95, 2, cumprod) # n year survival probability
ex_APC_BELm.boot.q95 <- life.exp(q_APC_BELm.boot.q95) # life expectancy

p_APC_BELm.boot.q05 <- 1-q_APC_BELm.boot.q05# 1 year survival probability
p0n_APC_BELm.boot.q05 <- apply(p_APC_BELm.boot.q05, 2, cumprod) # n year survival probability
ex_APC_BELm.boot.q05 <- life.exp(q_APC_BELm.boot.q05) # life expectancy

p_APC_BELf.q95 <- 1-q_APC_BELf.q95# 1 year survival probability
p0n_APC_BELf.q95 <- apply(p_APC_BELf.q95, 2, cumprod) # n year survival probability
ex_APC_BELf.q95 <- life.exp(q_APC_BELf.q95) # life expectancy

p_APC_BELf.q05 <- 1-q_APC_BELf.q05# 1 year survival probability
p0n_APC_BELf.q05 <- apply(p_APC_BELf.q05, 2, cumprod) # n year survival probability
ex_APC_BELf.q05 <- life.exp(q_APC_BELf.q05) # life expectancy

p_APC_BELf.boot.q95 <- 1-q_APC_BELf.boot.q95# 1 year survival probability
p0n_APC_BELf.boot.q95 <- apply(p_APC_BELf.boot.q95, 2, cumprod) # n year survival probability
ex_APC_BELf.boot.q95 <- life.exp(q_APC_BELf.boot.q95) # life expectancy

p_APC_BELf.boot.q05 <- 1-q_APC_BELf.boot.q05# 1 year survival probability
p0n_APC_BELf.boot.q05 <- apply(p_APC_BELf.boot.q05, 2, cumprod) # n year survival probability
ex_APC_BELf.boot.q05 <- life.exp(q_APC_BELf.boot.q05) # life expectancy

ex_APC.BELm.q05<-as.data.frame(ex_APC_BELm.q05)
ex_APC.BELm.q95<-as.data.frame(ex_APC_BELm.q95)
ex_APC.BELm.boot.q05<-as.data.frame(ex_APC_BELm.boot.q05)
ex_APC.BELm.boot.q95<-as.data.frame(ex_APC_BELm.boot.q95)
ex_APC.BELf.q05<-as.data.frame(ex_APC_BELf.q05)
ex_APC.BELf.q95<-as.data.frame(ex_APC_BELf.q95)
ex_APC.BELf.boot.q05<-as.data.frame(ex_APC_BELf.boot.q05)
ex_APC.BELf.boot.q95<-as.data.frame(ex_APC_BELf.boot.q95)

#### Present value of stocastic annuities
p_APC_BELm.st.ext <- 1-q_APC_BELm.st.ext # 1 year survival probability
ann_APC_BELm.st <- annuity.st(q_APC_BELm.st.ext,v_vect,0.9)
ann_APC_BELm.mean <- ann_APC_BELm.st[[1]]
ann_APC_BELm.q95 <- ann_APC_BELm.st[[2]]
ann_APC_BELm.q05 <- ann_APC_BELm.st[[3]]

# Male
par(mfrow=c(1, 2))
plot(ann_APC_BELm.mean["65",], type="l", lty=1) # evolution with t (same age)
lines(ann_APC_BELm.q95["65",], type="l", lty=2)
lines(ann_APC_BELm.q05["65",], type="l", lty=2)
plot(ann_APC_BELm.mean[,"2050"], type="l", lty=1) # evolution with x (same year)
lines(ann_APC_BELm.q95[,"2050"], type="l", lty=2)
lines(ann_APC_BELm.q05[,"2050"], type="l", lty=2)
dev.off()

p_APC_BELf.st.ext <- 1-q_APC_BELf.st.ext # 1 year survival probability
ann_APC_BELf.st <- annuity.st(q_APC_BELf.st.ext,v_vect,0.9)
ann_APC_BELf.mean <- ann_APC_BELf.st[[1]]
ann_APC_BELf.q95 <- ann_APC_BELf.st[[2]]
ann_APC_BELf.q05 <- ann_APC_BELf.st[[3]]

# Female
par(mfrow=c(1, 2))
plot(ann_APC_BELf.mean["65",], type="l", lty=1) # evolution with t (same age)
lines(ann_APC_BELf.q95["65",], type="l", lty=2)
lines(ann_APC_BELf.q05["65",], type="l", lty=2)
plot(ann_APC_BELf.mean[,"2050"], type="l", lty=1) # evolution with x (same year)
lines(ann_APC_BELf.q95[,"2050"], type="l", lty=2)
lines(ann_APC_BELf.q05[,"2050"], type="l", lty=2)
dev.off()

ann_APC.BELm.m <- as.data.frame(ann_APC_BELm.mean)
ann_APC.BELm.q95 <- as.data.frame(ann_APC_BELm.q95)
ann_APC.BELm.q05 <- as.data.frame(ann_APC_BELm.q05)
ann_APC.BELf.m <- as.data.frame(ann_APC_BELf.mean)
ann_APC.BELf.q95 <-as.data.frame(ann_APC_BELf.q95)
ann_APC.BELf.q05<- as.data.frame(ann_APC_BELf.q05)

### Present value of stocastic annuities : II method
ann_APC.BELm.c <- as.data.frame(annuity(q_APC_BELm.ext, v_vect))
ann_APC.BELm.q05q <- as.data.frame(annuity(q_APC_BELm.q95, v_vect))
ann_APC.BELm.q95q <- as.data.frame(annuity(q_APC_BELm.q05, v_vect))
ann_APC.BELf.c <- as.data.frame(annuity(q_APC_BELf.ext, v_vect))
ann_APC.BELf.q05q <- as.data.frame(annuity(q_APC_BELf.q95, v_vect))
ann_APC.BELf.q95q <- as.data.frame(annuity(q_APC_BELf.q05, v_vect))
