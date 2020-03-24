library(platexpress)
library(reshape2)
library(pracma)
library(ggplot2)

INPATH <- "C://Users/tpeng/OneDrive/OD__Universität/H__SynB/part_Axmann/git_repos - Kopie"
expnum <- 23 #20#

glcMaxConc = 22.2*6*0.9
aceMaxConc = 66.6*2*0.9

#======================================= functions =======================================

findSat = function(mat, npoints, timeOnCol = T){
  
  if(!timeOnCol){mat = t(mat)}
  
  endVal = ncol(mat)
  
  mat = melt(mat)
  mat$Var2 = as.numeric(mat$Var2)
  
  startVar = var(mat$value)
  minVar = startVar
  
  minVar_x = 1
  
  currRange = c(1,endVal)
  
  for(x in 1:(as.integer(endVal/npoints)+1)*2){
    
    evalsAll = round(linspace(x1 = currRange[1], x2 = currRange[2], n = npoints + 2)) #evaluated points
    evalsAll = unique(evalsAll)
    evals = evalsAll[c(-1, -length(evalsAll))]
    
    evals_var=
    sapply(evals, function(newStart){
      var(mat[mat$Var2>=newStart,"value"])
    })
    
    evals_var = evals_var/(endVal - evals)
    
    evals_best = which.min(evals_var)
    evals_min = evals_var[evals_best]
    
    if(evals_min<minVar){
      minVar_x = evals[evals_best]
      minVar = evals_min
      
      currRange = evalsAll[c(evals_best, evals_best+2)]
    }else{
      return(minVar_x)
    }
  }
  return(minVar_x)
}

bestFrontLm = function(y, x){
  
  df = data.frame(x,y)
  
  groups = sort(unique(x))
  
  lms = lapply(groups[-1], function(i){
    lm(y~x, data=df[x<=i,])
  })
  
  rsqs = sapply(lms, function(i){
    summary(i)$adj.r.squared
  })
  
  # print(rsqs)
  # return(lms)
  
  bestLm = which.max(rsqs)
  
  res = list(includedX = groups[1:(bestLm+1)],
             lm = lms[[bestLm]])
  return(res)
}

O2PercToConc = function(perc, p = 21278.25){#perc = percent of O2 saturation, p = partial pressure of O2 in Pa (21% * 101,325 Pa)
  #Henrys law : c_max = H_cp * p , max gas concentration in mol/m^3
  H_cp = 1.3 * 10^(-5) #mol/(m^3*Pa)
  
  c_max = H_cp * p * 10^(-3) #max O2 concentration in mol/l
  
  c = c_max*perc*10 #*1000 (mol -> mmol), *0.01 (percent -> fraction=
  
  return(c) #concentration in mmol/l
}

# O2concToRate = function(conc, c_max, kLa = 260){ #conc = O2 concentration in mmol/l, c_max = maximal O2 concentration after Henry-law, kLa = Oxygen transfer-coefficient, given by aeration setup
#   # for this BioLector the kLa is approximately 260 h^(-1) for 1ml and 500 rpm
#   rate = kLa * (c_max - conc) #mmol/(l h)
#   
#   return(rate)
# }

ggplotDf = function(mat, groups, times){
  mat = t(mat)
  mat = as.data.frame(mat)
  mat = cbind(groups, labels = row.names(mat), mat)
  mat = melt(mat, id.vars = c("groups", "labels"),   variable.name = "times")
  
  mat$labels = levels(mat$labels)[mat$labels]
  mat$times = times[mat$times]
  
  if(is.factor(mat$groups)){mat$groups = levels(mat$groups)[mat$groups]}
  
  return(mat)
}

#=========================================================================================

DATPATH <- file.path(INPATH,"data")
RESPATH <- file.path(INPATH,"results") # results

setwd(RESPATH)


expid <- paste0(expnum, "_BIQ980_2019_REDUCTION-1")
data.file <- file.path(DATPATH, paste0(expid,".csv"))
layout.file <- file.path(DATPATH, paste0(expid,"_layout.csv"))


#view full dataset
dat <- readExperiment(data.file, type="BioLectorPro", time.conversion=1/3600,
                      layout=layout.file, 
                      fields=c("strain","glc","ace","aTc"), 
                      afields=c("glc","ace","aTc"),blank.id="blank",
                      group1="glc", group2 = c("glc","amount"),
                      group2.color="color")

viewPlate(dat)

viewMap(map=dat$layout, nrow=6, ncol=8, color="color", text="amount", title="plate layout: glucose")
viewMap(map=dat$layout, nrow=6, ncol=8, color="ace.color", text="ace.amount", title="plate layout: acetate")


#======================================= glucose =======================================

for(i in list(c("glucose", "glc"), c("acetate")))
#read data
dat <- readExperiment(data.file, type="BioLectorPro", time.conversion=1/3600,
                      layout=layout.file, 
                      fields=c("strain","glc","ace","aTc"), 
                      afields=c("glc","ace","aTc"),
                      blank.id="blank",blank.data="Biomass",
                      skip.wells = c(paste0(toupper(letters[4:6]), rep(1:7,3))),
                      group1="glc", group2 = c("glc","amount"),
                      group2.color="color")

dat <- prettyData(dat, yids=c(ribof="Riboflavine",O2="DO(Pst3)",
                              scatter="Biomass", pH="pH(HP8)",
                              NADH="NADH - NADPH"))

#format data
dat$layout$amount <- (dat$layout$amount/max(dat$layout$amount)) * glcMaxConc
dat$layout$amount[dat$layout$strain=="blank"] <- NA
dat$layout$color[dat$layout$strain=="blank"] <- "#999999"


#cut off first lag phase
dat <- cutData(dat, xrng=c(0,23))


viewGroups(dat, yids="scatter", show.ci95=TRUE, lwd.orig=0)
viewGroups(dat, yids="ribof", show.ci95=TRUE, lwd.orig=0)
viewGroups(dat, yids="O2", show.ci95=TRUE, lwd.orig=0)
viewGroups(dat, yids="pH", show.ci95=TRUE, lwd.orig=0)
viewGroups(dat, yids="NADH", show.ci95=TRUE, lwd.orig=0)

#convert scatter to biomass
scatter <- getData(dat, "scatter") # a.u., arbitrary unit
gps <- 0.603 # mg/ml/scatter, biomass per scatter
cpg <- 0.474 # g/g, gram carbon per gram biomass, Folsom&Carlson 2015
mpc <- 1/12 # mol/g, mol carbon per gram carbon
biomass <- scatter * 1e3 * gps * cpg * mpc # mmol/L
dat <- addData(dat, "biomass", dat=biomass)


#get yield
cmol <- dat$layout$amount
names(cmol) <- dat$layout$well

cmol <- cmol[dat$layout$strain!="blank"] #exclude blanks

sapply(unique(dat$layout$amount)[!is.na(unique(dat$layout$amount))], function(i){
  
  
  subCols = dat$layout$amount == i
  subCols = subCols[order(dat$layout$well)]
  subCols[is.na(subCols)] = F
  
  subDat = dat$biomass$data[,subCols]
  
  findSat(subDat,10,F)
}) #find a good beginning for the saturation

bm <- apply(getData(dat, "biomass", xrng=c(15,20)),2,mean)[names(cmol)] #biomass between 15h and 20h

fitFull <- bestForcedLm(x = cmol, y = bm)
fit  = fitFull$lm
pars <- coef(fit)

par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.4,0), tcl=-.25)
plot(cmol, bm, xlab=paste0("glucose", ", C-mMol"), ylab="biomass, C-mMol")
abline(a=pars[1], b=pars[2])
legend("right",title="yield:", legend=paste(round(pars[2],2)))


#Oxygen
O2 = getData(dat, "O2")
O2 = O2/max(O2)*100

O2Conc = apply(O2, 2, O2PercToConc) #mmol/l

dat <- addData(dat, "O2conc", dat=O2Conc)


# O2AddRate = apply(O2Conc, 2, O2concToRate, c_max = O2PercToConc(100), kLa = 260) #mmol/(l h)
O2AddRate = (O2PercToConc(100) - O2Conc) * 260 #(c_max - c) * kLa (kLa = Oxygen transfer-coefficient, given by aeration setup) (mmol/(l h))

O2Change = (rbind(O2Conc[-1,], NA) - O2Conc) #O2 difference of two consecutive datapoints (mmol/l)
row.names(O2Change) = row.names(O2Conc)
O2Change = O2Change/ c(dat$Time[-1] - dat$Time[-168], NA) #devide by time difference (mmol/(l h))

O2MinusRate = O2Change - O2AddRate # difference in real and theoretical rate is consumption rate (mmol (l h))
O2MinusRate = O2MinusRate *10^(-3) # mmol/(ml h) (*1 ml) mmol/h

biomassG = scatter*gps*10^(-3) #g/ml (*1ml) -> g

O2rate = O2MinusRate/biomassG # O2 consumpion rate of the bacteria (mmol/(h g))

# max(O2rate[!is.infinite(O2rate)])
# boxplot(O2rate[!is.infinite(O2rate)])
# 
# plot(y = O2rate, x=matrix(dat$Time, nrow = nrow(O2rateNorm), ncol = ncol(O2rateNorm)))

dat <- addData(dat, "O2rate", dat=O2rate)

viewGroups(dat, yids="O2rate", show.ci95=TRUE, lwd.orig=0)

wells = dat$layout$amount[order(dat$layout$well)]
times = dat$Time

O2rate_df = ggplotDf(O2rate, wells, times)
O2rate_df = O2rate_df[complete.cases(O2rate_df),]

O2rate_plt = ggplot(O2rate_df[!O2rate_df$groups==0,], aes(x=times, y=value, group=as.factor(labels)))+
  geom_line()+
  facet_grid(groups~.)+
  geom_hline(yintercept = 0, linetype ="dashed")+
  scale_x_continuous(expand=c(0,0))+
  labs(x="time (in h)", y=expression("O_2 consumprion rate"))+
  theme_classic()

O2rate_plt + coord_cartesian(xlim = c(3,23), ylim= c(-50,0)) #Only after 3h, since the data looks suspicious before

#save data for glucose
dat_glc <- dat