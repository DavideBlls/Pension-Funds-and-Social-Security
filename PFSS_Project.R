library(demography)
library("StMoMo")
library("MTS")
library(devtools)
library(xlsx)#it allows to read directly an xlsx file
library(writexl)#it allows to write directly an xlsx file

"setwd _ _ _"
source("Functions.r") #uploading functions from an other script in the directory


##################################
######### Data setting ###########
##################################
a.min <- 0              # Start Age
a.max <- 100            # Max Age
omega <- 121            # Extreme Age
A.fit <-c(a.min:a.max)  # Age interval

# Fitting time interval
y.fit.min <- 1968
y.fit.max <- 2018
Y.fit <- c(y.fit.min:y.fit.max)

# Projectation interval
y.pred <- 100
Y.pred <- c((y.fit.max+1):(y.fit.max+y.pred))
n.sim <- 1000
n.boot <- 20
int_vect <- c(rep(0.02,omega-a.min)) # vector of interest rates (flat)
v_vect <- cumprod((1+int_vect)^-1)   # vector of discount factors


#####################################
########       BELGIUM      #########
### Demographic Data - Source HMD ###
#####################################
BELdata <- hmd.mx(country = "BEL", username = ("_ _ _ @ _.it"),
                  password = "_ _ _", label = "Belgium")

# I divide the population into males and females #
BELmStMoMoC <- StMoMoData(BELdata, series = "male")
BELfStMoMoC <- StMoMoData(BELdata, series = "female")

# from central exposed to risk to initial exposed to risk #
BELmStMoMoI <- central2initial(BELmStMoMoC)
BELfStMoMoI <- central2initial(BELfStMoMoC)

# specify the weight matrix #
wxt <- genWeightMat(ages = A.fit, years = Y.fit, clip = 4)

# death rates for male population #
BELmRates <- BELmStMoMoC$Dxt/BELmStMoMoC$Ext

# matrix of the males death rate #
BELmRates <- BELmRates[A.fit+1,tail(BELmStMoMoC$years+1,length(Y.fit))-BELmStMoMoC$years[1]]

# death rate for female population #
BELfRates <- BELfStMoMoC$Dxt/BELfStMoMoC$Ext

# matrix of the females death rate #
BELfRates <- BELfRates[A.fit+1,tail(BELfStMoMoC$years+1,length(Y.fit))-BELfStMoMoC$years[1]]


###################################
######## Models fitting ###########
###################################
LCfit_BELm <- fit(lc(link = "logit"), data = BELmStMoMoI, ages.fit = a.min:a.max,years.fit=y.fit.min:y.fit.max, wxt = wxt)
LCfit_BELf <- fit(lc(link = "logit"), data = BELfStMoMoI, ages.fit = a.min:a.max,years.fit=y.fit.min:y.fit.max, wxt = wxt)

# in RH model you need to specify an initial value for alpha and beta
# although this feature may make it better, it also makes it much slower in execution than the others
# we use the initial value of LC to initialize the algorithm
RHfit_BELm <- fit(rh(link = "logit", cohortAgeFun="1"), data = BELmStMoMoI, ages.fit = a.min:a.max,years.fit=y.fit.min:y.fit.max,
                  wxt = wxt, start.ax = LCfit_BELm$ax, start.bx = LCfit_BELm$bx, start.kt = LCfit_BELm$kt)
RHfit_BELf <- fit(rh(link = "logit", cohortAgeFun="1"), data = BELfStMoMoI, ages.fit = a.min:a.max,years.fit=y.fit.min:y.fit.max,
                  wxt = wxt, start.ax = LCfit_BELf$ax, start.bx = LCfit_BELf$bx, start.kt = LCfit_BELf$kt)
APCfit_BELm <- fit(apc(link = "logit"), data = BELmStMoMoI, ages.fit = a.min:a.max,years.fit=y.fit.min:y.fit.max, wxt = wxt)
APCfit_BELf <- fit(apc(link = "logit"), data = BELfStMoMoI, ages.fit = a.min:a.max,years.fit=y.fit.min:y.fit.max, wxt = wxt)
CBDfit_BELm <- fit(cbd(), data = BELmStMoMoI, ages.fit = a.min:a.max,years.fit=y.fit.min:y.fit.max, wxt = wxt)
CBDfit_BELf <- fit(cbd(), data = BELfStMoMoI, ages.fit = a.min:a.max,years.fit=y.fit.min:y.fit.max, wxt = wxt)
M7fit_BELm <- fit(m7(link = "logit"), data = BELmStMoMoI, ages.fit = a.min:a.max,years.fit=y.fit.min:y.fit.max, wxt = wxt)
M7fit_BELf <- fit(m7(link = "logit"), data = BELfStMoMoI, ages.fit = a.min:a.max,years.fit=y.fit.min:y.fit.max, wxt = wxt)
# the function PLAT is defined in the script functions.R
PLATfit_BELm <- fit(PLAT, data = BELmStMoMoI, ages.fit = a.min:a.max,years.fit=y.fit.min:y.fit.max, wxt = wxt)
PLATfit_BELf <- fit(PLAT, data = BELfStMoMoI, ages.fit = a.min:a.max,years.fit=y.fit.min:y.fit.max, wxt = wxt)


##################################
### Parameters representation #### 
##################################
plot(LCfit_BELm, nCol = 3)
plot(LCfit_BELf, nCol = 3)
plot(RHfit_BELm, nCol = 3)
plot(RHfit_BELf, nCol = 3)
plot(APCfit_BELm, parametricbx = F)
plot(APCfit_BELf, parametricbx = F)
plot(CBDfit_BELm, parametricbx = F)
plot(CBDfit_BELf, parametricbx = F)
plot(M7fit_BELm, parametricbx = F)
plot(M7fit_BELf, parametricbx = F)
plot(PLATfit_BELm, parametricbx = F)
plot(PLATfit_BELf, parametricbx = F)

#### lines below refer to Lee-Carter model - similar istructions can be written for the others
par(mfrow=c(2, 3))
plot(LCfit_BELm$ax, main=c(expression(paste("Male population - ",alpha[x]))))
plot(LCfit_BELm$bx, main=c(expression(paste("Male population - ",beta[x]))))
plot(t(LCfit_BELm$kt), main=c(expression(paste("Male population - ",kappa[t]))))
plot(LCfit_BELf$ax, main=c(expression(paste("Female population - ",alpha[x]))))
plot(LCfit_BELf$bx, main=c(expression(paste("Female population - ",beta[x]))))
plot(t(LCfit_BELf$kt), main=c(expression(paste("Female population - ",kappa[t]))))
dev.off()

######################################
######### Models comparison ##########
######### Log-likelihood   ###########
######################################

#Males
LCfit_BELm$loglik
CBDfit_BELm$loglik
RHfit_BELm$loglik
APCfit_BELm$loglik
M7fit_BELm$loglik
PLATfit_BELm$loglik
#Females
LCfit_BELf$loglik
CBDfit_BELf$loglik
RHfit_BELf$loglik
APCfit_BELf$loglik
M7fit_BELf$loglik
PLATfit_BELf$loglik

###### BIC ####
#Males
BIC(LCfit_BELm)
BIC(CBDfit_BELm)
BIC(RHfit_BELm)
BIC(APCfit_BELm)
BIC(M7fit_BELm)
BIC(PLATfit_BELm)
#Females
BIC(LCfit_BELf)
BIC(CBDfit_BELf)
BIC(RHfit_BELf)
BIC(APCfit_BELf)
BIC(M7fit_BELf)
BIC(PLATfit_BELf)

#####AIC#####
#MALE
AIC(LCfit_BELm)
AIC(CBDfit_BELm)
AIC(RHfit_BELm)
AIC(APCfit_BELm)
AIC(M7fit_BELm)
AIC(PLATfit_BELm)
#Females
AIC(LCfit_BELf)
AIC(CBDfit_BELf)
AIC(RHfit_BELf)
AIC(APCfit_BELf)
AIC(M7fit_BELf)
AIC(PLATfit_BELf)

###### Residuals
LCres_BELm <- residuals(LCfit_BELm)
LCres_BELf <- residuals(LCfit_BELf)
CBDres_BELm <- residuals(CBDfit_BELm)
CBDres_BELf <- residuals(CBDfit_BELf)
RHres_BELm <- residuals(RHfit_BELm)
RHres_BELf <- residuals(RHfit_BELf)
APCres_BELm <- residuals(APCfit_BELm)
APCres_BELf <- residuals(APCfit_BELf)
M7res_BELm <- residuals(M7fit_BELm)
M7res_BELf <- residuals(M7fit_BELf)
PLATres_BELm <- residuals(PLATfit_BELm)
PLATres_BELf <- residuals(PLATfit_BELf)

### Alternative residuals representations
TYPE="colourmap"#"scatter","signplot"
plot(LCres_BELm, type=TYPE)
plot(LCres_BELf, type=TYPE)
plot(CBDres_BELm, type=TYPE)
plot(CBDres_BELf, type=TYPE)
plot(RHres_BELm, type=TYPE)
plot(RHres_BELf, type=TYPE)
plot(APCres_BELm, type=TYPE)
plot(APCres_BELf, type=TYPE)
plot(M7res_BELm, type=TYPE)
plot(M7res_BELf, type=TYPE)
plot(PLATres_BELm, type=TYPE)
plot(PLATres_BELf, type=TYPE)
