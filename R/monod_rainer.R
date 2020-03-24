####################
##  MONOD'S MODEL ##
####################

INPATH <- "~/work/hhu_talks/biq980/uebung_201912/biq980_2019"

## plot parameters
fig.type <- "png"
W <- 3.5
H <- 1.5


## loading all experiment result files from BioLector
DATPATH <- file.path(INPATH,"data")
RESPATH <- file.path(INPATH,"results") # results
setwd(RESPATH)


## TEST MONOD Model

library(mugro)

## initial conditions and parameters
yinit <- c(s=270, y=1)
poinit <- c(phi=0, sin=0, mumax=.3, K=.5, Y=0.35, d=0)
time <- seq(0,24,.1)
## simulate ode from starting parameters
odem <- ode(yinit, time, ode_monod, poinit)
plot(odem[,"time"], odem[,"y"],col=2, log="y")


### COMPARE TO DATA
library(platexpress)
library(growthrates)

load(file.path(RESPATH, "23_BIQ980_2019_REDUCTION-1.rda"))

## get glucose data
dat <- dat_glc

dfg <- data2growthrates(dat, "biomass", plate=dat$layout)


## initial conditions and parameters
## NOTE: yield was calculated as 0.36 and maximal growth rate ca. 0.25
##       in the preprocessing script!
yinit <- c(s=270, y=1)
poinit <- c(phi=0, sin=0, mumax=.25, K=5, Y=0.35, d=0)

## NOTE: i ran the code below a couple of times and gradually fixed
##       more and more parameters until only initial conditions
##       were allowed to vary
fixed <- c("sin", "phi", "d", "Y", "mumax", "K")
fitted <- c(names(yinit),names(poinit))
fitted <- fitted[!fitted%in%fixed]

## simulate ode from starting parameters
well <- "C1"
odem <- ode(yinit, dfg$time[dfg$well==well], ode_monod, poinit)
plot(dfg$time[dfg$well==well], dfg$value[dfg$well==well])
lines(odem[,"time"], odem[,"y"],col=2)
plot(odem[,"time"], odem[,"y"],col=2)
lines(dfg$time[dfg$well==well], dfg$value[dfg$well==well])

## inspect modle results
plot(odem[,"time"], odem[,"s"])
lines(odem[,"time"], odem[,"y"],col=2)
mu <- poinit["mumax"] * odem[,"s"] / (odem[,"s"] + poinit["K"])
par(new=TRUE)
plot(odem[,"time"], mu, col=4, type="l", axes=FALSE)
axis(1)
axis(4)

## lower/upper boundaries for parameters
parms <- c(yinit, poinit)
lowp <- uppp <- parms
lowp[] <- 0
uppp[] <- 60
uppp["mumax"] <- 2
lowp["s"] <- 0
uppp["s"] <- 300
lowp["y"] <- 0.1
uppp["y"] <- 2.0
lowp["Y"] <- .01
uppp["Y"] <- 1

fits7 <- all_growthmodels(value ~ time | well, data = dfg,
                          FUN = grow_monod,
                          p = parms, lower = lowp, upper = uppp,
                          which = fitted, ncores=1) 
apply(coef(fits7),2,sd)

## TODO: check results in fits7@fits[[i]]@fit$message and
## fits7@fits[[i]]@fit$convergence
## "Relative error between `par' and the solution is at most `ptol'."

## ADD FITTED MODEL TO PLATEXPRESS
dat <- addModel(fits7, dat, ID="ODEMonod", replace=TRUE, col="#0000AA")

viewPlate(dat, yids=c("biomass","ODEMonod"),axes=c(2,4))

# plotdev("biolector_glucose_monod", width=W, height=W,
#         type=fig.type, res=300)
par(mfcol=c(2,1), mai=c(.3,.5,.1,.1), mgp=c(1.3,.4,0), tcl=-.25)
viewGroups(dat, yids="biomass", embed=TRUE, no.par=TRUE,log="",
           g1.legend=FALSE, g2.legend=FALSE,
           ylab="biomass, C-mMol", xlab="time, h")
viewGroups(dat, yids="ODEMonod", embed=TRUE, no.par=TRUE,log="",
           g1.legend=FALSE, g2.legend=FALSE,
           ylab=expression("ODE model:"~y(t)), xlab="time, h")
#dev.off()


## DOSE RESPONSE CURVES for FITTED PARAMETERS

## get parameters
coefs <- coef(fits7)

## match amounts and parameters
amnt <- dat$layout$amount
names(amnt) <-dat$layout$well
amnt <- amnt[names(amnt)%in%rownames(coefs)]

# plotdev("biolector_glucose_monod_parameters", width=W, height=W,
#         type=fig.type, res=300)
par(mfcol=c(1,1), mgp=c(1.3,.4,0),tcl=-.25,mai=c(.5,.5,.1,.1))
plot(amnt, coefs[names(amnt),"s"],
     ylab=expression("ODE model, fitted"~S(0)*", C-mMol"),
     xlab="glucose, C-mMol")
abline(a=0, b=1, lwd=.5)
##plot(amnt, coefs[names(amnt),"mumax"], ylab=expression(mu[max]*","~h^-1),
##     xlab="glucose, C-mMol")
##plot(amnt, coefs[names(amnt),"K"], ylab="K, C-mMol",xlab="glucose, C-mMol")
legend("bottomright",title="fixed:",
       paste(names(poinit[fixed]), poinit[fixed],sep="="))
# dev.off()


## compare invidual results

biomass <- getData(dat, "biomass")
model <- getData(dat, "ODEMonod")

## only experiment wells
gwells <- as.character(dat$layout$well[dat$layout$strain!="blank"])

## compare all values

## calculate global residuals over all experiments
Ntime <- nrow(model)
Nexp <- length(gwells)
N <- Ntime*Nexp
residuals <- model[,gwells] - biomass[,gwells]
RMSE <- sqrt(sum(c(residuals^2))/N)

## RMSE per experiment
resm <- apply(residuals,2,function(x) sqrt(sum(x^2)/Ntime))

plot(amnt, resm[names(amnt)]) # note: higher RSME at high carbon

## Pearson correlation
## global, with test for p-values
crg <- cor.test(c( model[,gwells]), c(biomass[,gwells]))
## and for each experiment
crm <- sapply(1:length(gwells), function(i) cor(model[,gwells[i]],
                                                biomass[,gwells[i]]))
names(crm) <- gwells
plot(amnt, crm[names(amnt)]) # note: higher cor at higher C

## NOTE: choosing a transparent color to get an idea
## about local density (overllaying points)
# plotdev("biolector_glucose_monod_cor", width=W, height=W,
#         type=fig.type, res=300)
par(mfcol=c(1,1), mgp=c(1.3,.4,0),tcl=-.25,mai=c(.5,.5,.1,.1))
plot(c(biomass[,gwells]), c(model[,gwells]), cex=.5, pch=19, col="#00000011",
     ylab="modelled biomass, C-mMol", xlab="measured biomass, C-mMol")
abline(a=0, b=1, col=2, lwd=3, lty=2)
legend("bottomright", title="Pearson r:", legend=round(crg$estimate,4))
# dev.off()

## compare individual curves to show where it deviates

well <- "A5"
for ( well in gwells ) { 
    # plotdev(paste0("biolector_glucose_monod_residuals_", well),
    #         width=W, height=W, type=fig.type, res=300)
    par(mfcol=c(1,1), mgp=c(1.3,.4,0),tcl=-.25,mai=c(.5,.5,.1,.5))
    plot(dat$Time, biomass[,well], #ylim=c(0,60),
         xlab="time, h", ylab="biomass, C-mMol", pch=19, cex=.5)
    lines(dat$Time, model[,well], col=2, lwd=2)
    legend("topleft", legend=c("measured", "modelled"), col=1:2,
           pch=c(19,NA), lty=c(NA,1), pt.cex=.5)
    par(new=TRUE)
    plot(dat$Time, model[,well]-biomass[,well], col=4,
         axes=FALSE, xlab=NA, ylab=NA)
    axis(4, col=4, col.axis=4)
    mtext("residuals", 4, par("mgp")[1], col=4)
    abline(h=0, col=4)
    # dev.off()
}

