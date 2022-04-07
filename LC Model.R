"setwd _ _ _"
source("Functions.r") #uploading functions from an other script in the directory

####################################
#######Lee-Carter model############ 
####################################
################### Forecast
LCfor_BELm <- forecast(LCfit_BELm, h=y.pred)
### Alternative istructions
#LCfor_BELmArima <- forecast(LCfit_BELm, h=y.pred, kt.method = "iarima",  kt.order = c(1, 1, 2))  # ARIMA order specified
#LCfor_BELmBArima <- forecast(LCfit_BELm, h=y.pred, kt.method = "iarima",  kt.order = NULL)  # auto.arima
LCfor_BELf <- forecast(LCfit_BELf, h=y.pred)
plot(LCfor_BELm, only.kt = TRUE)
plot(LCfor_BELf, only.kt = TRUE)
#####
par(mfrow=c(1, 2), oma=c(2,0,0,0))
plot(c(a.min:a.max),log(LCfor_BELm$rates[,y.pred]), bty="l", xlab="Ages", ylab="Male death rates (log scale)",
     ylim=c(min(log(LCfor_BELm$rates[,y.pred])),max(log((LCfit_BELm$Dxt/LCfit_BELm$Ext)[,length(Y.fit)]))), type="l")
lines(c(a.min:a.max),log((LCfit_BELm$Dxt/LCfit_BELm$Ext)[,length(Y.fit)]), type="p")
plot(c(a.min:a.max),log(LCfor_BELf$rates[,y.pred]), bty="l", xlab="Ages", ylab="Female death rates (log scale)",
     ylim=c(min(log(LCfor_BELf$rates[,y.pred])),max(log((LCfit_BELf$Dxt/LCfit_BELf$Ext)[,length(Y.fit)]))), type="l")
lines(c(a.min:a.max),log((LCfit_BELf$Dxt/LCfit_BELf$Ext)[,length(Y.fit)]), type="p")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("last observed", "last projected"), pch=c(1,-1), lty=c(0,1), lwd=c(1,1), col=c(1:1), cex=0.9)
dev.off()

rates_LCfit_BELm <- fitted(LCfit_BELm, type = "rates")
rates_LCfor_BELm <- LCfor_BELm$rates
#rates_LC_BELm <- cbind(rates_LCfit_BELm,rates_LCfor_BELm)
rates_LC_BELm <- cbind(BELmRates,rates_LCfor_BELm)
q_LC_BELm <- 1- exp(-rates_LC_BELm)
q_LC_BELm.ext  <- extrapolation.fit(q_LC_BELm)

rates_LCfit_BELf <- fitted(LCfit_BELf, type = "rates")
rates_LCfor_BELf <- LCfor_BELf$rates
#rates_LC_BELf <- cbind(rates_LCfit_BELf,rates_LCfor_BELf)
rates_LC_BELf <- cbind(BELfRates,rates_LCfor_BELf)
q_LC_BELf <- 1- exp(-rates_LC_BELf)
q_LC_BELf.ext  <- extrapolation.fit(q_LC_BELf)

#write.table(q_LC_BELm.ext,file="q_LC.BELm.txt",sep=",")
#write.table(q_LC_BELf.ext,file="q_LC.BELf.txt",sep=",")
q_LC.BELm<-as.data.frame(q_LC_BELm.ext)
q_LC.BELf<-as.data.frame(q_LC_BELf.ext)

plot(extractCohort(q_LC_BELm.ext, cohort = 1980), type = "p", log = "y", xlab = "age", ylab = "q(x)",
     main = "Mortality rates for a cohort (log scale)")
lines(extractCohort(q_LC_BELm, cohort = 1980))
plot(extractCohort(log(q_LC_BELm.ext/(1-q_LC_BELm.ext)), cohort = 1980), type = "p", xlab = "age", ylab = "logit(q(x))",
     main = "Logit transform of mortality rates for a cohort")
lines(extractCohort(log(q_LC_BELm/(1-q_LC_BELm)), cohort = 1980))

### 1 year survival probability, n year survival probability and life expectancy
p_LC_BELm.ext <- 1-q_LC_BELm.ext # 1 year survival probability
p0n_LC_BELm.ext <- apply(p_LC_BELm.ext, 2, cumprod) # n year survival probability
ex_LC_BELm.ext <- life.exp(q_LC_BELm.ext)

p_LC_BELf.ext <- 1-q_LC_BELf.ext # 1 year survival probability
p0n_LC_BELf.ext <- apply(p_LC_BELf.ext, 2, cumprod) # n year survival probability
ex_LC_BELf.ext <- life.exp(q_LC_BELf.ext)

#write.table(ex_LC_BELm.ext,file="ex_LC.BELm.txt",sep=",")
#write.table(ex_LC_BELf.ext,file="ex_LC.BELf.txt",sep=",")
p_LC.BELm<-as.data.frame(p_LC_BELm.ext)
p0n_LC.BELm<-as.data.frame(p0n_LC_BELm.ext)
ex_LC.BELm<-as.data.frame(ex_LC_BELm.ext)
p_LC.BELf<-as.data.frame(p_LC_BELf.ext)
p0n_LC.BELf<-as.data.frame(p0n_LC_BELf.ext)
ex_LC.BELf<-as.data.frame(ex_LC_BELf.ext)


##### SIMULATION WITH RANDOM WALK WITH DRIFT #####
#MALE
LCsim_BELm.mrwd <- simulate(LCfit_BELm, nsim = n.sim, h=y.pred)
rates_LC_BELm.st <- LCsim_BELm.mrwd$rates
q_LC_BELm.st <- 1- exp(-rates_LC_BELm.st)
q_LC_BELm.st.ext <-  extrapolation.sim(q_LC_BELm.st)

par(mfrow=c(1, 2))
plot(LCfit_BELm$years, LCfit_BELm$kt[1, ], xlim = range(LCfit_BELm$years, LCsim_BELm.mrwd$kt.s$years),
     ylim = range(LCfit_BELm$kt, LCsim_BELm.mrwd$kt.s$sim), type = "l", xlab = "year", ylab = "kt",
     main = "Lee-Carter: Period index (mrwd)")
matlines(LCsim_BELm.mrwd$kt.s$years, LCsim_BELm.mrwd$kt.s$sim[1, , ], type = "l", lty = 1)

plot(LCfit_BELm$years, (LCfit_BELm$Dxt / LCfit_BELm$Ext)["65", ], xlim = range(LCfit_BELm$years, LCsim_BELm.mrwd$years),
     ylim = range((LCfit_BELm$Dxt / LCfit_BELm$Ext)["65", ], LCsim_BELm.mrwd$rates["65", , ]), type = "l", xlab = "year", ylab = "rate",
     main = "Lee-Carter: Simulated mortality rates at age 65")
matlines(LCsim_BELm.mrwd$years, LCsim_BELm.mrwd$rates["65", , ], type = "l", lty = 1)
dev.off()


#FEMALE
LCsim_BELf.mrwd <- simulate(LCfit_BELf, nsim = n.sim, h=y.pred)
rates_LC_BELf.st <- LCsim_BELf.mrwd$rates
q_LC_BELf.st <- 1- exp(-rates_LC_BELf.st)
q_LC_BELf.st.ext <-  extrapolation.sim(q_LC_BELf.st)

par(mfrow=c(1, 2))
plot(LCfit_BELf$years, LCfit_BELf$kt[1, ], xlim = range(LCfit_BELf$years, LCsim_BELf.mrwd$kt.s$years),
     ylim = range(LCfit_BELf$kt, LCsim_BELf.mrwd$kt.s$sim), type = "l", xlab = "year", ylab = "kt",
     main = "Lee-Carter: Period index (mrwd)")
matlines(LCsim_BELf.mrwd$kt.s$years, LCsim_BELf.mrwd$kt.s$sim[1, , ], type = "l", lty = 1)

plot(LCfit_BELf$years, (LCfit_BELf$Dxt / LCfit_BELf$Ext)["65", ], xlim = range(LCfit_BELf$years, LCsim_BELf.mrwd$years),
     ylim = range((LCfit_BELf$Dxt / LCfit_BELf$Ext)["65", ], LCsim_BELf.mrwd$rates["65", , ]), type = "l", xlab = "year", ylab = "rate",
     main = "Lee-Carter: Simulated mortality rates at age 65")
matlines(LCsim_BELf.mrwd$years, LCsim_BELf.mrwd$rates["65", , ], type = "l", lty = 1)



library("fanplot")
probs = c(2.5, 10, 25, 50, 75, 90, 97.5)
par(mfrow=c(1, 2))
matplot(LCfit_BELm$years, t(q_LC_BELm[c("65", "75", "85"),c(1:length(Y.fit)) ]),
        xlim = c(y.fit.min, (y.fit.min+y.pred)), ylim = c(0.001, 0.2), pch = 20, col = "black",
        log = "y", xlab = "year", ylab = "male mortality rate (log scale)")
fan(t(LCsim_BELm.mrwd$rates["65", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("black", "white")), ln = NULL)
fan(t(LCsim_BELm.mrwd$rates["75", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("red", "white")), ln = NULL)
fan(t(LCsim_BELm.mrwd$rates["85", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
text(y.fit.min+4, q_LC_BELm[c("65", "75", "85"), 30],
     labels = c("x = 65", "x = 75", "x = 85"))
matplot(LCfit_BELf$years, t(q_LC_BELf[c("65", "75", "85"),c(1:length(Y.fit)) ]),
        xlim = c(y.fit.min, (y.fit.min+y.pred)), ylim = c(0.001, 0.2), pch = 20, col = "black",
        log = "y", xlab = "year", ylab = "female mortality rate (log scale)")
fan(t(LCsim_BELf.mrwd$rates["65", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("black", "white")), ln = NULL)
fan(t(LCsim_BELf.mrwd$rates["75", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("red", "white")), ln = NULL)
fan(t(LCsim_BELf.mrwd$rates["85", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
text(y.fit.min+4, q_LC_BELf[c("65", "75", "85"), 30],
     labels = c("x = 65", "x = 75", "x = 85"))
dev.off()


####################################
####### Bootstrap #################
####################################
LCboot_BELm <- bootstrap(LCfit_BELm, nBoot = n.boot, type = "semiparametric")
plot(LCboot_BELm, nCol = 3)
LCsim_BELm.boot <- simulate(LCboot_BELm, nsim = n.sim/n.boot, h = y.pred)
rates_LC_BELm.boot.st <- LCsim_BELm.boot$rates
q_LC_BELm.boot.st <- 1- exp(-rates_LC_BELm.boot.st)
q_LC_BELm.boot.st.ext <-  extrapolation.sim(q_LC_BELm.boot.st)

LCboot_BELf <- bootstrap(LCfit_BELf, nBoot = n.boot, type = "semiparametric")
plot(LCboot_BELf, nCol = 3)
LCsim_BELf.boot <- simulate(LCboot_BELf, nsim = n.sim/n.boot, h = y.pred)
rates_LC_BELf.boot.st <- LCsim_BELf.boot$rates
q_LC_BELf.boot.st <- 1- exp(-rates_LC_BELf.boot.st)
q_LC_BELf.boot.st.ext <-  extrapolation.sim(q_LC_BELf.boot.st)


####################################
#### Confidence intervals at 90% ####
####################################
##### 1 year death probability
conf.lev <- 0.9
q_LC_BELm.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_BELm.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_BELm.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_BELm.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
    q_LC_BELm.q95[j,k] <- quantile(q_LC_BELm.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_LC_BELm.q05[j,k] <- quantile(q_LC_BELm.st.ext[j,k,], probs=(1-conf.lev)/2)
    q_LC_BELm.boot.q95[j,k] <- quantile(q_LC_BELm.boot.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_LC_BELm.boot.q05[j,k] <- quantile(q_LC_BELm.boot.st.ext[j,k,], probs=(1-conf.lev)/2)
  }
}

q_LC_BELf.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_BELf.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_BELf.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_BELf.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
    q_LC_BELf.q95[j,k] <- quantile(q_LC_BELf.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_LC_BELf.q05[j,k] <- quantile(q_LC_BELf.st.ext[j,k,], probs=(1-conf.lev)/2)
    q_LC_BELf.boot.q95[j,k] <- quantile(q_LC_BELf.boot.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_LC_BELf.boot.q05[j,k] <- quantile(q_LC_BELf.boot.st.ext[j,k,], probs=(1-conf.lev)/2)
  }
}

q_LC.BELm.q95<-as.data.frame(q_LC_BELm.q95)
q_LC.BELm.q05<-as.data.frame(q_LC_BELm.q05)
q_LC.BELm.boot.q95<-as.data.frame(q_LC_BELm.boot.q95)
q_LC.BELm.boot.q05<-as.data.frame(q_LC_BELm.boot.q05)
q_LC.BELf.q95<-as.data.frame(q_LC_BELf.q95)
q_LC.BELf.q05<-as.data.frame(q_LC_BELf.q05)
q_LC.BELf.boot.q95<-as.data.frame(q_LC_BELf.boot.q95)
q_LC.BELf.boot.q05<-as.data.frame(q_LC_BELf.boot.q05)

plot(q_LC_BELm.ext["65",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_LC_BELm.q95["65",], x=Y.pred, type="l", lty=2)
lines(q_LC_BELm.q05["65",], x=Y.pred, type="l", lty=2)
lines(q_LC_BELm.boot.q95["65",], x=Y.pred, type="l", lty=2, col=2)
lines(q_LC_BELm.boot.q05["65",], x=Y.pred, type="l", lty=2, col=2)

plot(q_LC_BELm.ext["20",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_LC_BELm.q95["20",], x=Y.pred, type="l", lty=2)
lines(q_LC_BELm.q05["20",], x=Y.pred, type="l", lty=2)
lines(q_LC_BELm.boot.q95["20",], x=Y.pred, type="l", lty=2, col=2)
lines(q_LC_BELm.boot.q05["20",], x=Y.pred, type="l", lty=2, col=2)

#aggiunto io
plot(q_LC_BELf.ext["20",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_LC_BELf.q95["20",], x=Y.pred, type="l", lty=2)
lines(q_LC_BELf.q05["20",], x=Y.pred, type="l", lty=2)
lines(q_LC_BELf.boot.q95["20",], x=Y.pred, type="l", lty=2, col=2)
lines(q_LC_BELf.boot.q05["20",], x=Y.pred, type="l", lty=2, col=2)

plot(q_LC_BELf.ext["65",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_LC_BELf.q95["65",], x=Y.pred, type="l", lty=2)
lines(q_LC_BELf.q05["65",], x=Y.pred, type="l", lty=2)
lines(q_LC_BELf.boot.q95["65",], x=Y.pred, type="l", lty=2, col=2)
lines(q_LC_BELf.boot.q05["65",], x=Y.pred, type="l", lty=2, col=2)



##### 1 year survival probability, n year survival probability and life expectancy confidence interval
p_LC_BELm.q95 <- 1-q_LC_BELm.q95 # 1 year survival probability
p0n_LC_BELm.q95 <- apply(p_LC_BELm.q95, 2, cumprod) # n year survival probability
ex_LC_BELm.q95 <- life.exp(q_LC_BELm.q95) # life expectancy

p_LC_BELm.q05 <- 1-q_LC_BELm.q05 # 1 year survival probability
p0n_LC_BELm.q05 <- apply(p_LC_BELm.q05, 2, cumprod) # n year survival probability
ex_LC_BELm.q05 <- life.exp(q_LC_BELm.q05) # life expectancy

p_LC_BELm.boot.q95 <- 1-q_LC_BELm.boot.q95 # 1 year survival probability
p0n_LC_BELm.boot.q95 <- apply(p_LC_BELm.boot.q95, 2, cumprod) # n year survival probability
ex_LC_BELm.boot.q95 <- life.exp(q_LC_BELm.boot.q95) # life expectancy

p_LC_BELm.boot.q05 <- 1-q_LC_BELm.boot.q05 # 1 year survival probability
p0n_LC_BELm.boot.q05 <- apply(p_LC_BELm.boot.q05, 2, cumprod) # n year survival probability
ex_LC_BELm.boot.q05 <- life.exp(q_LC_BELm.boot.q05) # life expectancy

p_LC_BELf.q95 <- 1-q_LC_BELf.q95 # 1 year survival probability
p0n_LC_BELf.q95 <- apply(p_LC_BELf.q95, 2, cumprod) # n year survival probability
ex_LC_BELf.q95 <- life.exp(q_LC_BELf.q95) # life expectancy

p_LC_BELf.q05 <- 1-q_LC_BELf.q05 # 1 year survival probability
p0n_LC_BELf.q05 <- apply(p_LC_BELf.q05, 2, cumprod) # n year survival probability
ex_LC_BELf.q05 <- life.exp(q_LC_BELf.q05) # life expectancy

p_LC_BELf.boot.q95 <- 1-q_LC_BELf.boot.q95 # 1 year survival probability
p0n_LC_BELf.boot.q95 <- apply(p_LC_BELf.boot.q95, 2, cumprod) # n year survival probability
ex_LC_BELf.boot.q95 <- life.exp(q_LC_BELf.boot.q95) # life expectancy

p_LC_BELf.boot.q05 <- 1-q_LC_BELf.boot.q05 # 1 year survival probability
p0n_LC_BELf.boot.q05 <- apply(p_LC_BELf.boot.q05, 2, cumprod) # n year survival probability
ex_LC_BELf.boot.q05 <- life.exp(q_LC_BELf.boot.q05) # life expectancy

ex_LC.BELm.q05<-as.data.frame(ex_LC_BELm.q05)
ex_LC.BELm.q95<-as.data.frame(ex_LC_BELm.q95)
ex_LC.BELm.boot.q05<-as.data.frame(ex_LC_BELm.boot.q05)
ex_LC.BELm.boot.q95<-as.data.frame(ex_LC_BELm.boot.q95)
ex_LC.BELf.q05<-as.data.frame(ex_LC_BELf.q05)
ex_LC.BELf.q95<-as.data.frame(ex_LC_BELf.q95)
ex_LC.BELf.boot.q05<-as.data.frame(ex_LC_BELf.boot.q05)
ex_LC.BELf.boot.q95<-as.data.frame(ex_LC_BELf.boot.q95)


#################################################
##### Present value of stochastic annuities #####
#################################################
p_LC_BELm.st.ext <- 1-q_LC_BELm.st.ext # 1 year survival probability
ann_LC_BELm.st <- annuity.st(q_LC_BELm.st.ext,v_vect,0.9)
ann_LC_BELm.mean <- ann_LC_BELm.st[[1]]
ann_LC_BELm.q95 <- ann_LC_BELm.st[[2]]
ann_LC_BELm.q05 <- ann_LC_BELm.st[[3]]

########
par(mfrow=c(1, 2))
plot(ann_LC_BELm.mean["65",], type="l", lty=1) # evolution with t (same age)
lines(ann_LC_BELm.q95["65",], type="l", lty=2)
lines(ann_LC_BELm.q05["65",], type="l", lty=2)
plot(ann_LC_BELm.mean[,"2050"], type="l", lty=1) # evolution with x (same year)
lines(ann_LC_BELm.q95[,"2050"], type="l", lty=2)
lines(ann_LC_BELm.q05[,"2050"], type="l", lty=2)
dev.off()
########

p_LC_BELf.st.ext <- 1-q_LC_BELf.st.ext # 1 year survival probability
ann_LC_BELf.st <- annuity.st(q_LC_BELf.st.ext,v_vect,0.9)
ann_LC_BELf.mean <- ann_LC_BELf.st[[1]]
ann_LC_BELf.q95 <- ann_LC_BELf.st[[2]]
ann_LC_BELf.q05 <- ann_LC_BELf.st[[3]]

########
par(mfrow=c(1, 2))
plot(ann_LC_BELf.mean["65",], type="l", lty=1) # evolution with t (same age)
lines(ann_LC_BELf.q95["65",], type="l", lty=2)
lines(ann_LC_BELf.q05["65",], type="l", lty=2)
plot(ann_LC_BELf.mean[,"2050"], type="l", lty=1) # evolution with x (same year)
lines(ann_LC_BELf.q95[,"2050"], type="l", lty=2)
lines(ann_LC_BELf.q05[,"2050"], type="l", lty=2)
dev.off()
########

ann_LC.BELm.m <- as.data.frame(ann_LC_BELm.mean)
ann_LC.BELm.q95 <- as.data.frame(ann_LC_BELm.q95)
ann_LC.BELm.q05 <- as.data.frame(ann_LC_BELm.q05)
ann_LC.BELf.m <- as.data.frame(ann_LC_BELf.mean)
ann_LC.BELf.q95 <-as.data.frame(ann_LC_BELf.q95)
ann_LC.BELf.q05<- as.data.frame(ann_LC_BELf.q05)

##### Present value of stocastic annuities : II method
ann_LC.BELm.c <- as.data.frame(annuity(q_LC_BELm.ext, v_vect))
ann_LC.BELm.q05q <- as.data.frame(annuity(q_LC_BELm.q95, v_vect))
ann_LC.BELm.q95q <- as.data.frame(annuity(q_LC_BELm.q05, v_vect))
ann_LC.BELf.c <- as.data.frame(annuity(q_LC_BELf.ext, v_vect))
ann_LC.BELf.q05q <- as.data.frame(annuity(q_LC_BELf.q95, v_vect))
ann_LC.BELf.q95q <- as.data.frame(annuity(q_LC_BELf.q05, v_vect))
#################################



q_LC_BELf.1955 <- extractCohort(q_LC_BELf.ext, cohort = 1955)
q_LC_BELf.st.1955 <- extractCohort(q_LC_BELf.st.ext, cohort = 1955)
q_LC_BELf.1980 <- extractCohort(q_LC_BELf.ext, cohort = 1980)
q_LC_BELf.st.1980 <- extractCohort(q_LC_BELf.st.ext, cohort = 1980)

delta.i <- 0.0001
q_LC_BELf.1955.up <- q_LC_BELf.1955.up+delta.i
q_LC_BELf.1980.up <- q_LC_BELf.1980.up+delta.i

q_LC_BELf.1955.dw <- q_LC_BELf.1955.up-delta.i
q_LC_BELf.1980.dw <- q_LC_BELf.1980.up-delta.i

V.ann.1955 <- annuity.coh(q_LC_BELf.1955, v_vect)
V.ann.1955.up <- annuity.coh(q_LC_BELf.1955.up, v_vect)
V.ann.1955.dw <- annuity.coh(q_LC_BELf.1955.dw, v_vect)
dur.ann.1955 <- (V.ann.1955.up -V.ann.1955.dw)/(2*V.ann.1955*delta.i)

V.LI.1980 <- LI.coh(q_LC_BELf.1980, v_vect)
V.LI.1980.up <- LI.coh(q_LC_BELf.1980.up, v_vect)
V.LI.1980.dw <- LI.coh(q_LC_BELf.1980.dw, v_vect)
dur.LI.1980 <- (V.LI.1980.up -V.LI.1980.dw)/(2*V.LI.1980*delta.i)

w.LI <- -dur.ann.1955["65"]/(dur.LI.1980["40"]-dur.ann.1955["65"])
w.ann <- 1- w.LI
V.ptf <- 1000
q.LI <- V.ptf*w.LI/V.LI.1980["40"]
q.ann <- V.ptf*w.ann/V.ann.1955["65"]

V.ann.st.1955 <-q.ann*annuity.st.coh(q_LC_BELf.st.1955, v_vect)["65",]
V.LI.st.1980 <-q.LI*LI.st.coh(q_LC_BELf.st.1980, v_vect)["40",]
V.ptf.st <- V.ann.st.1955 + V.LI.st.1980

cv.ann <- sd(V.ann.st.1955)/mean(V.ann.st.1955)
cv.LI <- sd(V.LI.st.1980)/mean(V.LI.st.1980)
cv.ptf <- sd(V.ptf.st)/mean(V.ptf.st)

