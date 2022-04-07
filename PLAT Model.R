"setwd _ _ _"
source("Functions.r") 

##############################
######### Plat model #########
##############################
# Forecast
PLATfor_BELm <- forecast(PLATfit_BELm, h=y.pred, gc.order = c(2, 0, 0))
# default ARIMA order for gc. Alternative order can be choosen via auto.arima
PLATfor_BELf <- forecast(PLATfit_BELf, h=y.pred, gc.order = c(2, 0, 0))
plot(PLATfor_BELm, only.kt = TRUE)
plot(PLATfor_BELf, only.kt = TRUE)
plot(PLATfor_BELm, only.gc = TRUE)
plot(PLATfor_BELf, only.gc = TRUE)
#####
par(mfrow=c(1, 2), oma=c(2,0,0,0))
plot(c(a.min:a.max),log(PLATfor_BELm$rates[,y.pred]), bty="l", xlab="Ages", ylab="Male death rates (log scale)",
  ylim=c(min(log(PLATfor_BELm$rates[,y.pred])),max(log((PLATfit_BELm$Dxt/PLATfit_BELm$Ext)[,length(Y.fit)]))), type="l")
lines(c(a.min:a.max),log((PLATfit_BELm$Dxt/PLATfit_BELm$Ext)[,length(Y.fit)]), type="p")
plot(c(a.min:a.max),log(PLATfor_BELf$rates[,y.pred]), bty="l", xlab="Ages", ylab="Female death rates (log scale)",
  ylim=c(min(log(PLATfor_BELf$rates[,y.pred])),max(log((PLATfit_BELf$Dxt/PLATfit_BELf$Ext)[,length(Y.fit)]))), type="l")
lines(c(a.min:a.max),log((PLATfit_BELf$Dxt/PLATfit_BELf$Ext)[,length(Y.fit)]), type="p")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("last observed", "last projected"), pch=c(1,-1), lty=c(0,1), lwd=c(1,1), col=c(1:1), cex=0.9)
dev.off()
########

rates_PLATfit_BELm <- fitted(PLATfit_BELm, type = "rates")
ifelse(is.na(rates_PLATfit_BELm),BELmRates,rates_PLATfit_BELm)
rates_PLATfor_BELm <- PLATfor_BELm$rates
#rates_PLAT_BELm <- cbind(rates_PLATfit_BELm,rates_PLATfor_BELm)
rates_PLAT_BELm <- cbind(BELmRates,rates_PLATfor_BELm)
q_PLAT_BELm <- 1- exp(-rates_PLAT_BELm)
q_PLAT_BELm.ext  <- extrapolation.fit(q_PLAT_BELm)

rates_PLATfit_BELf <- fitted(PLATfit_BELf, type = "rates")
ifelse(is.na(rates_PLATfit_BELf),BELfRates,rates_PLATfit_BELf)
rates_PLATfor_BELf <- PLATfor_BELf$rates
#rates_PLAT_BELf <- cbind(rates_PLATfit_BELf,rates_PLATfor_BELf)
rates_PLAT_BELf <- cbind(BELfRates,rates_PLATfor_BELf)
q_PLAT_BELf <- 1- exp(-rates_PLAT_BELf)
q_PLAT_BELf.ext  <- extrapolation.fit(q_PLAT_BELf)

#write.table(q_PLAT_BELm.ext,file="q_PLAT.BELm.txt",sep=",")
#write.table(q_PLAT_BELf.ext,file="q_PLAT.BELf.txt",sep=",")
q_PLAT.BELm<-as.data.frame(q_PLAT_BELm.ext)
q_PLAT.BELf<-as.data.frame(q_PLAT_BELf.ext)

plot(extractCohort(q_PLAT_BELm.ext, cohort = 1980), type = "p", log = "y", xlab = "age", ylab = "q(x)",
 main = "Mortality rates for a cohort (log scale)")
lines(extractCohort(q_PLAT_BELm, cohort = 1980))

plot(extractCohort(log(q_PLAT_BELm.ext/(1-q_PLAT_BELm.ext)), cohort = 1980), type = "p", xlab = "age", ylab = "logit(q(x))",
 main = "Logit transform of mortality rates for a cohort")
lines(extractCohort(log(q_PLAT_BELm/(1-q_PLAT_BELm)), cohort = 1980))

### 1 year survival probability, n year survival probability and life expectancy
p_PLAT_BELm.ext <- 1-q_PLAT_BELm.ext # 1 year survival probability
p0n_PLAT_BELm.ext <- apply(p_PLAT_BELm.ext, 2, cumprod) # n year survival probability
ex_PLAT_BELm.ext <- life.exp(q_PLAT_BELm.ext)

p_PLAT_BELf.ext <- 1-q_PLAT_BELf.ext # 1 year survival probability
p0n_PLAT_BELf.ext <- apply(p_PLAT_BELf.ext, 2, cumprod) # n year survival probability
ex_PLAT_BELf.ext <- life.exp(q_PLAT_BELf.ext)

#write.table(ex_PLAT_BELm.ext,file="ex_PLAT.BELm.txt",sep=",")
#write.table(ex_PLAT_BELf.ext,file="ex_PLAT.BELf.txt",sep=",")
p_PLAT.BELm<-as.data.frame(p_PLAT_BELm.ext)
p0n_PLAT.BELm<-as.data.frame(p0n_PLAT_BELm.ext)
ex_PLAT.BELm<-as.data.frame(ex_PLAT_BELm.ext)
p_PLAT.BELf<-as.data.frame(p_PLAT_BELf.ext)
p0n_PLAT.BELf<-as.data.frame(p0n_PLAT_BELf.ext)
ex_PLAT.BELf<-as.data.frame(ex_PLAT_BELf.ext)

### SIMULATION WITH RANDOM WALK WITH DRIFT
PLATsim_BELm.mrwd <- simulate(PLATfit_BELm, nsim = n.sim, h=y.pred, gc.order = c(2, 0, 0))
### Alternative istructions
#PLATsim_BELmArima <- forecast(PLATfit_BELm, nsim = n.sim, h=y.pred, kt.method = "iarima",  kt.order = c(1, 1, 2))  # ARIMA order specified
#PLATsim_BELmBArima <- forecast(PLATfit_BELm, nsim = n.sim, h=y.pred, kt.method = "iarima",  kt.order = NULL)  # auto.arima
rates_PLAT_BELm.st <- PLATsim_BELm.mrwd$rates
q_PLAT_BELm.st <- 1- exp(-rates_PLAT_BELm.st)
q_PLAT_BELm.st.ext <-  extrapolation.sim(q_PLAT_BELm.st)

PLATsim_BELf.mrwd <- simulate(PLATfit_BELf, nsim = n.sim, h=y.pred, gc.order = c(2, 0, 0))
rates_PLAT_BELf.st <- PLATsim_BELf.mrwd$rates
q_PLAT_BELf.st <- 1- exp(-rates_PLAT_BELf.st)
q_PLAT_BELf.st.ext <-  extrapolation.sim(q_PLAT_BELf.st)

par(mfrow=c(1, 2))
plot(PLATfit_BELm$years, PLATfit_BELm$kt[1, ], xlim = range(PLATfit_BELm$years, PLATsim_BELm.mrwd$kt.s$years),
  ylim = range(PLATfit_BELm$kt, PLATsim_BELm.mrwd$kt.s$sim), type = "l", xlab = "year", ylab = "kt",
  main = "PLAT: Period index (mrwd)")
matlines(PLATsim_BELm.mrwd$kt.s$years, PLATsim_BELm.mrwd$kt.s$sim[1, , ], type = "l", lty = 1)
plot(PLATfit_BELm$years, (PLATfit_BELm$Dxt / PLATfit_BELm$Ext)["65", ], xlim = range(PLATfit_BELm$years, PLATsim_BELm.mrwd$years),
  ylim = range((PLATfit_BELm$Dxt / PLATfit_BELm$Ext)["65", ], PLATsim_BELm.mrwd$rates["65", , ]), type = "l", xlab = "year", ylab = "rate",
  main = "PLAT: Simulated mortality rates at age 65")
matlines(PLATsim_BELm.mrwd$years, PLATsim_BELm.mrwd$rates["65", , ], type = "l", lty = 1)
dev.off()
########

library("fanplot")
probs = c(2.5, 10, 25, 50, 75, 90, 97.5)
par(mfrow=c(1, 2))
matplot(PLATfit_BELm$years, t(q_PLAT_BELm[c("65", "75", "85"),c(1:length(Y.fit)) ]),
 xlim = c(y.fit.min, (y.fit.min+y.pred)), ylim = c(0.001, 0.2), pch = 20, col = "black",
 log = "y", xlab = "year", ylab = "male mortality rate (log scale)")
fan(t(PLATsim_BELm.mrwd$rates["65", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
 fan.col = colorRampPalette(c("black", "white")), ln = NULL)
fan(t(PLATsim_BELm.mrwd$rates["75", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
 fan.col = colorRampPalette(c("red", "white")), ln = NULL)
fan(t(PLATsim_BELm.mrwd$rates["85", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
 fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
text(y.fit.min+4, q_PLAT_BELm[c("65", "75", "85"), 30],
 labels = c("x = 65", "x = 75", "x = 85"))
matplot(PLATfit_BELf$years, t(q_PLAT_BELf[c("65", "75", "85"),c(1:length(Y.fit)) ]),
 xlim = c(y.fit.min, (y.fit.min+y.pred)), ylim = c(0.001, 0.2), pch = 20, col = "black",
 log = "y", xlab = "year", ylab = "female mortality rate (log scale)")
fan(t(PLATsim_BELf.mrwd$rates["65", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
 fan.col = colorRampPalette(c("black", "white")), ln = NULL)
fan(t(PLATsim_BELf.mrwd$rates["75", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
 fan.col = colorRampPalette(c("red", "white")), ln = NULL)
fan(t(PLATsim_BELf.mrwd$rates["85", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
 fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
text(y.fit.min+4, q_PLAT_BELf[c("65", "75", "85"), 30],
 labels = c("x = 65", "x = 75", "x = 85"))
dev.off()

### Bootstrap
PLATboot_BELm <- bootstrap(PLATfit_BELm, nBoot = n.boot, type = "semiparametric")
plot(PLATboot_BELm, nCol = 3)
PLATsim_BELm.boot <- simulate(PLATboot_BELm, nsim = n.sim/n.boot, h = y.pred, gc.order = c(2, 0, 0))
rates_PLAT_BELm.boot.st <- PLATsim_BELm.boot$rates
q_PLAT_BELm.boot.st <- 1- exp(-rates_PLAT_BELm.boot.st)
q_PLAT_BELm.boot.st.ext <-  extrapolation.sim(q_PLAT_BELm.boot.st)

PLATboot_BELf <- bootstrap(PLATfit_BELf, nBoot = n.boot, type = "semiparametric")
plot(PLATboot_BELf, nCol = 3)
PLATsim_BELf.boot <- simulate(PLATboot_BELf, nsim = n.sim/n.boot, h = y.pred, gc.order = c(2, 0, 0))
rates_PLAT_BELf.boot.st <- PLATsim_BELf.boot$rates
q_PLAT_BELf.boot.st <- 1- exp(-rates_PLAT_BELf.boot.st)
q_PLAT_BELf.boot.st.ext <-  extrapolation.sim(q_PLAT_BELf.boot.st)

### Confidence intervals at 90%
# 1 year death probability
conf.lev <- 0.9
q_PLAT_BELm.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_BELm.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_BELm.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_BELm.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
      q_PLAT_BELm.q95[j,k] <- quantile(q_PLAT_BELm.st.ext[j,k,], probs=0.5+conf.lev/2)
      q_PLAT_BELm.q05[j,k] <- quantile(q_PLAT_BELm.st.ext[j,k,], probs=(1-conf.lev)/2)
      q_PLAT_BELm.boot.q95[j,k] <- quantile(q_PLAT_BELm.boot.st.ext[j,k,], probs=0.5+conf.lev/2)
      q_PLAT_BELm.boot.q05[j,k] <- quantile(q_PLAT_BELm.boot.st.ext[j,k,], probs=(1-conf.lev)/2)
  }
}

q_PLAT_BELf.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_BELf.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_BELf.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_BELf.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
      q_PLAT_BELf.q95[j,k] <- quantile(q_PLAT_BELf.st.ext[j,k,], probs=0.5+conf.lev/2)
      q_PLAT_BELf.q05[j,k] <- quantile(q_PLAT_BELf.st.ext[j,k,], probs=(1-conf.lev)/2)
      q_PLAT_BELf.boot.q95[j,k] <- quantile(q_PLAT_BELf.boot.st.ext[j,k,], probs=0.5+conf.lev/2)
      q_PLAT_BELf.boot.q05[j,k] <- quantile(q_PLAT_BELf.boot.st.ext[j,k,], probs=(1-conf.lev)/2)
  }
}

q_PLAT.BELm.q95<-as.data.frame(q_PLAT_BELm.q95)
q_PLAT.BELm.q05<-as.data.frame(q_PLAT_BELm.q05)
q_PLAT.BELm.boot.q95<-as.data.frame(q_PLAT_BELm.boot.q95)
q_PLAT.BELm.boot.q05<-as.data.frame(q_PLAT_BELm.boot.q05)
q_PLAT.BELf.q95<-as.data.frame(q_PLAT_BELf.q95)
q_PLAT.BELf.q05<-as.data.frame(q_PLAT_BELf.q05)
q_PLAT.BELf.boot.q95<-as.data.frame(q_PLAT_BELf.boot.q95)
q_PLAT.BELf.boot.q05<-as.data.frame(q_PLAT_BELf.boot.q05)

plot(q_PLAT_BELm.ext["65",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_PLAT_BELm.q95["65",], x=Y.pred, type="l", lty=2)
lines(q_PLAT_BELm.q05["65",], x=Y.pred, type="l", lty=2)
lines(q_PLAT_BELm.boot.q95["65",], x=Y.pred, type="l", lty=2, col=2)
lines(q_PLAT_BELm.boot.q05["65",], x=Y.pred, type="l", lty=2, col=2)

plot(q_PLAT_BELm.ext["20",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_PLAT_BELm.q95["20",], x=Y.pred, type="l", lty=2)
lines(q_PLAT_BELm.q05["20",], x=Y.pred, type="l", lty=2)
lines(q_PLAT_BELm.boot.q95["20",], x=Y.pred, type="l", lty=2, col=2)
lines(q_PLAT_BELm.boot.q05["20",], x=Y.pred, type="l", lty=2, col=2)


#scritto io
plot(q_PLAT_BELf.ext["65",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_PLAT_BELf.q95["65",], x=Y.pred, type="l", lty=2)
lines(q_PLAT_BELf.q05["65",], x=Y.pred, type="l", lty=2)
lines(q_PLAT_BELf.boot.q95["65",], x=Y.pred, type="l", lty=2, col=2)
lines(q_PLAT_BELf.boot.q05["65",], x=Y.pred, type="l", lty=2, col=2)

plot(q_PLAT_BELf.ext["20",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_PLAT_BELf.q95["20",], x=Y.pred, type="l", lty=2)
lines(q_PLAT_BELf.q05["20",], x=Y.pred, type="l", lty=2)
lines(q_PLAT_BELf.boot.q95["20",], x=Y.pred, type="l", lty=2, col=2)
lines(q_PLAT_BELf.boot.q05["20",], x=Y.pred, type="l", lty=2, col=2)

####################################

##### 1 year survival probability, n year survival probability and life expectancy confidence interval
p_PLAT_BELm.q95 <- 1-q_PLAT_BELm.q95# 1 year survival probability
p0n_PLAT_BELm.q95 <- apply(p_PLAT_BELm.q95, 2, cumprod) # n year survival probability
ex_PLAT_BELm.q95 <- life.exp(q_PLAT_BELm.q95) # life expectancy

p_PLAT_BELm.q05 <- 1-q_PLAT_BELm.q05# 1 year survival probability
p0n_PLAT_BELm.q05 <- apply(p_PLAT_BELm.q05, 2, cumprod) # n year survival probability
ex_PLAT_BELm.q05 <- life.exp(q_PLAT_BELm.q05) # life expectancy

p_PLAT_BELm.boot.q95 <- 1-q_PLAT_BELm.boot.q95# 1 year survival probability
p0n_PLAT_BELm.boot.q95 <- apply(p_PLAT_BELm.boot.q95, 2, cumprod) # n year survival probability
ex_PLAT_BELm.boot.q95 <- life.exp(q_PLAT_BELm.boot.q95) # life expectancy

p_PLAT_BELm.boot.q05 <- 1-q_PLAT_BELm.boot.q05# 1 year survival probability
p0n_PLAT_BELm.boot.q05 <- apply(p_PLAT_BELm.boot.q05, 2, cumprod) # n year survival probability
ex_PLAT_BELm.boot.q05 <- life.exp(q_PLAT_BELm.boot.q05) # life expectancy

p_PLAT_BELf.q95 <- 1-q_PLAT_BELf.q95# 1 year survival probability
p0n_PLAT_BELf.q95 <- apply(p_PLAT_BELf.q95, 2, cumprod) # n year survival probability
ex_PLAT_BELf.q95 <- life.exp(q_PLAT_BELf.q95) # life expectancy

p_PLAT_BELf.q05 <- 1-q_PLAT_BELf.q05# 1 year survival probability
p0n_PLAT_BELf.q05 <- apply(p_PLAT_BELf.q05, 2, cumprod) # n year survival probability
ex_PLAT_BELf.q05 <- life.exp(q_PLAT_BELf.q05) # life expectancy

p_PLAT_BELf.boot.q95 <- 1-q_PLAT_BELf.boot.q95# 1 year survival probability
p0n_PLAT_BELf.boot.q95 <- apply(p_PLAT_BELf.boot.q95, 2, cumprod) # n year survival probability
ex_PLAT_BELf.boot.q95 <- life.exp(q_PLAT_BELf.boot.q95) # life expectancy

p_PLAT_BELf.boot.q05 <- 1-q_PLAT_BELf.boot.q05# 1 year survival probability
p0n_PLAT_BELf.boot.q05 <- apply(p_PLAT_BELf.boot.q05, 2, cumprod) # n year survival probability
ex_PLAT_BELf.boot.q05 <- life.exp(q_PLAT_BELf.boot.q05) # life expectancy

ex_PLAT.BELm.q05<-as.data.frame(ex_PLAT_BELm.q05)
ex_PLAT.BELm.q95<-as.data.frame(ex_PLAT_BELm.q95)
ex_PLAT.BELm.boot.q05<-as.data.frame(ex_PLAT_BELm.boot.q05)
ex_PLAT.BELm.boot.q95<-as.data.frame(ex_PLAT_BELm.boot.q95)
ex_PLAT.BELf.q05<-as.data.frame(ex_PLAT_BELf.q05)
ex_PLAT.BELf.q95<-as.data.frame(ex_PLAT_BELf.q95)
ex_PLAT.BELf.boot.q05<-as.data.frame(ex_PLAT_BELf.boot.q05)
ex_PLAT.BELf.boot.q95<-as.data.frame(ex_PLAT_BELf.boot.q95)

### Present value of stocastic annuities
p_PLAT_BELm.st.ext <- 1-q_PLAT_BELm.st.ext # 1 year survival probability
ann_PLAT_BELm.st <- annuity.st(q_PLAT_BELm.st.ext,v_vect,0.9)
ann_PLAT_BELm.mean <- ann_PLAT_BELm.st[[1]]
ann_PLAT_BELm.q95 <- ann_PLAT_BELm.st[[2]]
ann_PLAT_BELm.q05 <- ann_PLAT_BELm.st[[3]]

########
par(mfrow=c(1, 2))
plot(ann_PLAT_BELm.mean["65",], type="l", lty=1) # evolution with t (same age)
lines(ann_PLAT_BELm.q95["65",], type="l", lty=2)
lines(ann_PLAT_BELm.q05["65",], type="l", lty=2)
plot(ann_PLAT_BELm.mean[,"2050"], type="l", lty=1) # evolution with x (same year)
lines(ann_PLAT_BELm.q95[,"2050"], type="l", lty=2)
lines(ann_PLAT_BELm.q05[,"2050"], type="l", lty=2)
dev.off()
########

p_PLAT_BELf.st.ext <- 1-q_PLAT_BELf.st.ext # 1 year survival probability
ann_PLAT_BELf.st <- annuity.st(q_PLAT_BELf.st.ext,v_vect,0.9)
ann_PLAT_BELf.mean <- ann_PLAT_BELf.st[[1]]
ann_PLAT_BELf.q95 <- ann_PLAT_BELf.st[[2]]
ann_PLAT_BELf.q05 <- ann_PLAT_BELf.st[[3]]

########
par(mfrow=c(1, 2))
plot(ann_PLAT_BELf.mean["65",], type="l", lty=1) # evolution with t (same age)
lines(ann_PLAT_BELf.q95["65",], type="l", lty=2)
lines(ann_PLAT_BELf.q05["65",], type="l", lty=2)
plot(ann_PLAT_BELf.mean[,"2050"], type="l", lty=1) # evolution with x (same year)
lines(ann_PLAT_BELf.q95[,"2050"], type="l", lty=2)
lines(ann_PLAT_BELf.q05[,"2050"], type="l", lty=2)
dev.off()
########

ann_PLAT.BELm.m <- as.data.frame(ann_PLAT_BELm.mean)
ann_PLAT.BELm.q95 <- as.data.frame(ann_PLAT_BELm.q95)
ann_PLAT.BELm.q05 <- as.data.frame(ann_PLAT_BELm.q05)
ann_PLAT.BELf.m <- as.data.frame(ann_PLAT_BELf.mean)
ann_PLAT.BELf.q95 <-as.data.frame(ann_PLAT_BELf.q95)
ann_PLAT.BELf.q05<- as.data.frame(ann_PLAT_BELf.q05)

##### Present value of stocastic annuities : II method
ann_PLAT.BELm.c <- as.data.frame(annuity(q_PLAT_BELm.ext, v_vect))
ann_PLAT.BELm.q05q <- as.data.frame(annuity(q_PLAT_BELm.q95, v_vect))
ann_PLAT.BELm.q95q <- as.data.frame(annuity(q_PLAT_BELm.q05, v_vect))
ann_PLAT.BELf.c <- as.data.frame(annuity(q_PLAT_BELf.ext, v_vect))
ann_PLAT.BELf.q05q <- as.data.frame(annuity(q_PLAT_BELf.q95, v_vect))
ann_PLAT.BELf.q95q <- as.data.frame(annuity(q_PLAT_BELf.q05, v_vect))
#################################

