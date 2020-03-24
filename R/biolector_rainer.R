#######################################################
## MAKE SURE YOU INSTALLED THE LATEST PLATEXPRESS!!! ##
#######################################################
library(platexpress)

##########################################################
## ADAPT THIS PATH TO THE LOCATION OF THE GIT DIRECTORY ##
##########################################################
INPATH <- "C://Users/tpeng/OneDrive/OD__Universität/H__SynB/part_Axmann/git_repos - Kopie"


## plot parameters

PLOT2PDF <- FALSE # TRUE # re-direct all plots to PDF

fig.type <- "png"
W <- 3.5
H <- 1.5


### PARSING and INITIAL ANALYSIS OF BIOLECTOR DATA

## loading all experiment result files from BioLector

DATPATH <- file.path(INPATH,"data")
RESPATH <- file.path(INPATH,"results") # results

setwd(RESPATH)

expnum <- 23 #20# 
expid <- paste0(expnum, "_BIQ980_2019_REDUCTION-1")
data.file <- file.path(DATPATH, paste0(expid,".csv"))
layout.file <- file.path(DATPATH, paste0(expid,"_layout.csv"))

## plot all to PDF
if ( PLOT2PDF )
  pdf(file.path(RESPATH, paste0(expid,".pdf")), width=12, height=8)


### BioLector Pro E.coli experiment with either glc or ace growth

## inspect full experiment
dat <- readExperiment(data.file, type="BioLectorPro", time.conversion=1/3600,
                      layout=layout.file, 
                      fields=c("strain","glc","ace","aTc"), 
                      afields=c("glc","ace","aTc"),blank.id="blank",
                      group1="glc", group2 = c("glc","amount"),
                      group2.color="color")
viewPlate(dat)

## check plate-layout
viewMap(map=dat$layout, nrow=6, ncol=8, color="color", text="amount", title="plate layout: glucose")


## load glucose data - skip all other wells
dat <- readExperiment(data.file, type="BioLectorPro", time.conversion=1/3600,
                      layout=layout.file, 
                      fields=c("strain","glc","ace","aTc"), 
                      afields=c("glc","ace","aTc"),
                      blank.id="blank",blank.data="Biomass",
                      skip.wells = c(paste0(toupper(letters[4:6]), rep(1:7,3))),
                      group1="glc", group2 = c("glc","amount"),
                      group2.color="color")


## parse raw data
dat <- prettyData(dat, yids=c(ribof="Riboflavine",O2="DO(Pst3)",
                              scatter="Biomass", pH="pH(HP8)",
                              NADH="NADH - NADPH"))


## FIX AMOUNTS to actual values, in C-mMol
dat$layout$amount <- dat$layout$amount/1.2 * 270
dat$layout$amount[dat$layout$strain=="blank"] <- NA
dat$layout$color[dat$layout$strain=="blank"] <- "#999999"

## inspect the data
viewPlate(dat)
## check plate-layout
viewMap(map=dat$layout, nrow=6, ncol=8, color="color", text="amount", title="plate layout: glucose")

## plot to file for script!
plotdev("biolector_layout_glucose", width=W, height=1.5*H,
        type=fig.type, res=300)
par(mai=c(.3,.3,.2,.05), mgp=c(1.3,.4,0),tcl=-.25)
viewMap(map=dat$layout, nrow=6, ncol=8, color="color",
        text="amount", title="plate layout: glucose")
dev.off()

## cut to exclude long stationary phase
dat <- cutData(dat, xrng=c(0,23))


viewGroups(dat, yids="scatter", show.ci95=TRUE, lwd.orig=0)
viewGroups(dat, yids="ribof", show.ci95=TRUE, lwd.orig=0)
viewGroups(dat, yids="O2", show.ci95=TRUE, lwd.orig=0)
viewGroups(dat, yids="pH", show.ci95=TRUE, lwd.orig=0)
viewGroups(dat, yids="NADH", show.ci95=TRUE, lwd.orig=0)

## CONVERT SCATTER TO BIOMASS: to C-mMol
## using measurements by Sergej

scatter <- getData(dat, "scatter") # a.u., arbitrary unit
gps <- 0.603 # mg/ml/scatter, biomass per scatter
cpg <- 0.474 # g/g, gram carbon per gram biomass, Folsom&Carlson 2015
mpc <- 1/12 # mol/g, mol carbon per gram carbon
biomass <- scatter * 1e3 * gps * cpg * mpc # mmol/L
dat <- addData(dat, "biomass", dat=biomass)

## YIELD on GLUCOSE

cmol <- dat$layout$amount
names(cmol) <- dat$layout$well

## without blanks
cmol <- cmol[dat$layout$strain!="blank"]

## get mean biomass between 15 and 20h
bm <- apply(getData(dat, "biomass", xrng=c(15,20)),2,mean)[names(cmol)]

plot(cmol, bm)

## fit for substrate range before saturation
fit <- lm(bm[cmol<150] ~  cmol[cmol<150])
pars <- coef(fit)

## plot yield calculation and fit
plotdev(paste0("biolector_yield_", "glc"),
        type=fig.type, width=W, height=W)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.4,0), tcl=-.25)
plot(cmol, bm, xlab=paste0("glucose", ", C-mMol"), ylab="biomass, C-mMol")
abline(a=pars[1], b=pars[2])
legend("right",title="yield:", legend=paste(round(pars[2],2)))
dev.off()

## CONVERT O2% to mmol O2
## using Henry's law

## CALCULATE O2 RATE
## using kLa of BioLector

## TODO: what is the maximal O2 consumption rate?
## compare to Andersen et al. 1981: 20 mmol/h/g

## save glucose data
dat_glc <- dat

## load acetate data - skip all other wells
dat <- readExperiment(data.file, type="BioLectorPro", time.conversion=1/3600,
                      layout=layout.file, 
                      fields=c("strain","glc","ace","aTc"), 
                      afields=c("glc","ace","aTc"),
                      blank.id="blank", blank.data=c("Biomass"),
                      skip.wells = c(paste0(toupper(letters[1:3]), rep(1:7,3))),
                      group1="ace", group2 = c("ace","ace.amount"),
                      group2.color="ace.color")


## parse raw data
dat <- prettyData(dat, yids=c(ribof="Riboflavine",O2="DO(Pst3)",
                              scatter="Biomass", pH="pH(HP8)",
                              NADH="NADH - NADPH"))

## FIX AMOUNTS to actual values
dat$layout$ace.amount <- dat$layout$ace.amount/1.2 * 270
dat$layout$ace.amount[dat$layout$strain=="blank"] <- NA
dat$layout$ace.color[dat$layout$strain=="blank"] <- "#999999"

scatter <- getData(dat, "scatter") # a.u., arbitrary unit
gps <- 0.603 # mg/ml/scatter, biomass per scatter
cpg <- 0.474 # g/g, gram carbon per gram biomass
mpc <- 1/12 # mol/g, mol carbon per gram carbon
biomass <- scatter * 1e3 * gps * cpg * mpc # mmol/L
dat <- addData(dat, "biomass", dat=biomass)

## inspect the data
viewPlate(dat)
## check plate-layout
viewMap(map=dat$layout, nrow=6, ncol=8, color="ace.color", text="ace.amount")

## plot to file for script!
plotdev("biolector_layout_acetate", width=W, height=1.5*H,
        type=fig.type, res=300)
par(mai=c(.3,.3,.2,.05), mgp=c(1.3,.4,0),tcl=-.25)
viewMap(map=dat$layout, nrow=6, ncol=8, color="ace.color", text="ace.amount", title="plate layout: acetate")
dev.off()

## view by replicate grouping
viewGroups(dat, yids="scatter",show.ci95=TRUE, lwd.orig=0)
viewGroups(dat, yids="biomass",show.ci95=TRUE, lwd.orig=0)
viewGroups(dat, yids="ribof", show.ci95=TRUE, lwd.orig=0)
viewGroups(dat, yids="O2", show.ci95=TRUE, lwd.orig=0)
viewGroups(dat, yids="pH", show.ci95=TRUE, lwd.orig=0)
viewGroups(dat, yids="NADH", show.ci95=TRUE, lwd.orig=0)

dat_ace <- dat

save(dat_ace, dat_glc, file=file.path(RESPATH, paste0(expid,".rda")))

if ( PLOT2PDF )
    dev.off()

## nicely design a plot

for ( nutrient in c("ace", "glc") ) {

    ## get data
    sdat <- get(paste0("dat_", nutrient))
    
    plotdev(paste0("biolector_growth_", nutrient),
            type=fig.type, width=W, height=6*H)
    par(mfcol=c(6,1), mai=c(.3,.35,.1,.1), mgp=c(1.3,.4,0), tcl=-.25)
    did <- ifelse(nutrient=="glc","",paste0(nutrient,"."))
    viewMap(map=sdat$layout, nrow=6, ncol=8,
            color=paste0(did,"color"),
            text=paste0(did,"amount"),
            title=paste("plate layout:", nutrient))
    viewGroups(sdat, yids="scatter",  show.ci95=FALSE, lwd.orig=0, log="",
               embed=TRUE, no.par=TRUE, g1.legend=FALSE,g2.legend=FALSE,
               ylab="scatter, a.u.", xlab="time, h", ylim=c(0,6))
    viewGroups(sdat, yids="ribof", show.ci95=FALSE, lwd.orig=0,
               embed=TRUE, no.par=TRUE, g1.legend=FALSE, g2.legend=FALSE,
               ylab="riboflavin fluor., a.u.", xlab="time, h", ylim=c(0,2.5))
    viewGroups(sdat, yids="O2", show.ci95=FALSE,lwd.orig=0,
           embed=TRUE,no.par=TRUE, g1.legend=FALSE, g2.legend=FALSE,
           ylab="diss. O2, %", xlab="time, h", ylim=c(0,120))
    viewGroups(sdat, yids="pH", show.ci95=FALSE,lwd.orig=0,
               embed=TRUE,no.par=TRUE, g1.legend=FALSE, g2.legend=FALSE,
               ylab="pH", xlab="time, h", ylim=c(6,7.2))
    viewGroups(sdat, yids="NADH", show.ci95=FALSE,lwd.orig=0,
               embed=TRUE,no.par=TRUE,g1.legend=FALSE, g2.legend=FALSE,
               ylab="NADH fluor., a.u.", xlab="time, h", ylim=c(0,.8))
    dev.off()
    
}

### GROWTH RATES WITH DPSEG

for ( nutrient in c("glc", "ace") ) {

    ## get data
    dat <- get(paste0("dat_", nutrient))

    ## smooth scatter data
    scat <- getData(dat, "biomass")
    smscat <- apply(scat, 2, ma, 5)

    ## inspect individual well
    id <- "B5"
    if ( nutrient=="ace" ) id <- "E5"
    plot(scat[,id])
    lines(smscat[,id],col=2, type="b", cex=.5)

    dat <- addData(dat, ID="smscat", smscat)
    
    dpseg <- dpseg_plate(dat, "smscat",P=.0001)
    
    ## add to plate express object to inspect
    mdat <- addModel(dpseg, dat, ID="smscat.dpseg")
    viewPlate(mdat, yids=c("smscat","smscat.dpseg"), log="y")
    viewGroups(mdat, yids=c("smscat.dpseg"),log="y")
    
    ## add growth rates
    mdat <- addModel(dpseg, mdat, ID="mu.dpseg", add.slopes=TRUE, col="#FF0000")
    viewPlate(mdat, yids=c("smscat","mu.dpseg"), log="")

    ## inspect a single well:
    
### NOTE: requires knowledge of the data structures
### use RStudio data browser to inspect
    id <- "A1"
    if ( nutrient=="ace" ) id <- "E5"

    par(mfcol=c(1,1), mai=c(.5,.5,.1,.5),mgp=c(1.3,.4,0),tcl=-.25)
    plot(mdat$smscat$data[,id],log="y",cex=.5)
    lines(mdat$smscat.dpseg$data[,id], lwd=2)
    par(new=TRUE)
    plot(mdat$mu.dpseg$data[,id], axes=FALSE, type="l",ylab=NA)
    abline(v=dpseg[[id]]$segments[,"start"])
    abline(v=dpseg[[id]]$segments[,"end"], col=2)
    axis(4)
    mtext("growth rate, h-1", 4, par("mgp")[1])

    ## inspect to get growth rates
    
    viewPlate(mdat, yids=c("smscat","mu.dpseg"))

    ## REPORT FIGURE
    plotdev(paste0("biolector_dpseg_",nutrient), width=W, height=3*H,
            type=fig.type, res=300)
    par(mai=c(.3,.5,.2,.05), mgp=c(1.3,.4,0),tcl=-.25)
    par(mfcol=c(2,1))
    viewGroups(mdat, yids=c("biomass"), show.ci95=FALSE,lwd.orig=0,
               xlab="time, h", ylab="biomass, c-mMol",
               log="y", ylim=c(1,100),
               embed=TRUE,no.par=TRUE, g1.legend=FALSE, g2.legend=FALSE)
    viewGroups(mdat, yids=c("mu.dpseg"), show.ci95=FALSE, lwd.orig=0,
               ylim=c(-.01,.32),
               xlab="time, h", ylab=expression("growth rate"~mu*","~h^-1),
               embed=TRUE,no.par=TRUE, g1.legend=FALSE, g2.legend=FALSE)
    abline(h=0, lwd=2)
    dev.off()
    
    ## TODO - platexpress: fix orig colors!
    colors.fixed <- FALSE
    if ( colors.fixed ) {
        par(mfcol=c(2,1))
        viewGroups(mdat, yids=c("smscat"), log="y",embed=TRUE,
                   show.ci95=FALSE, lwd.orig=1)
        viewGroups(mdat, yids=c("mu.dpseg"), embed=TRUE, ylim=c(-.01,.22),
                   g2.legend=FALSE, show.ci95=FALSE, lwd.orig=1)
        abline(h=0, lwd=2)
    }
}


for ( nutrient in c("glc", "ace") ) {

    ## get data
    dat <- get(paste0("dat_", nutrient))
 
    dpseg <- dpseg_plate(dat, "ribof",P=.0001)
    
    ## add to plate express object to inspect
    mdat <- addModel(dpseg, dat, ID="ribof.dpseg")
    viewPlate(mdat, yids=c("ribof","ribof.dpseg"), log="y")
    viewGroups(mdat, yids=c("ribof.dpseg"),log="y")
    
    ## add growth rates
    mdat <- addModel(dpseg, mdat, ID="mu.dpseg", add.slopes=TRUE, col="#FF0000")
    viewPlate(mdat, yids=c("ribof","mu.dpseg"), log="")

    ## inspect a single well:
    
### NOTE: requires knowledge of the data structures
### use RStudio data browser to inspect
    id <- "B5"
    if ( nutrient=="ace" ) id <- "E5"

    par(mfcol=c(1,1), mai=c(.5,.5,.1,.5),mgp=c(1.3,.4,0),tcl=-.25)
    plot(mdat$ribof$data[,id],log="y",cex=.5)
    lines(mdat$ribof.dpseg$data[,id], lwd=2)
    par(new=TRUE)
    plot(mdat$mu.dpseg$data[,id], axes=FALSE, type="l",ylab=NA)
    abline(v=dpseg[[id]]$segments[,"start"])
    abline(v=dpseg[[id]]$segments[,"end"], col=2)
    axis(4)
    mtext("growth rate, h-1", 4, par("mgp")[1])

    ## inspect to get growth rates
    viewPlate(mdat, yids=c("ribof","mu.dpseg"))
    par(mfcol=c(2,1))
    viewGroups(mdat, yids=c("ribof"), log="y",embed=TRUE, show.ci95=FALSE)
    viewGroups(mdat, yids=c("mu.dpseg"), embed=TRUE, ylim=c(-.01,.32),
               g2.legend=FALSE, show.ci95=FALSE)
    abline(h=0, lwd=2)

    ## TODO - platexpress: fix orig colors!
    colors.fixed <- FALSE
    if ( colors.fixed ) {
        par(mfcol=c(2,1))
        viewGroups(mdat, yids=c("ribof"), log="y",embed=TRUE,
                   show.ci95=FALSE, lwd.orig=1)
        viewGroups(mdat, yids=c("mu.dpseg"), embed=TRUE, ylim=c(-.01,.22),
                   g2.legend=FALSE, show.ci95=FALSE, lwd.orig=1)
        abline(h=0, lwd=2)
    }
}

## TODO:
## dose response curves: growth rate and yield vs. glucose vs. acetate
## respiration rate in mmol O2 per g DCW per hour

