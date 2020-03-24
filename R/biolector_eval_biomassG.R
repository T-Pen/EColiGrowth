INPATHBM <- "tpeng/OneDrive/OD__Universität/BA__Axmann"
expidBM <- "27_Ecoli_2020_REDUCTION-1"

if(PC){INPATHBMHeader = "D://"}else{INPATHBMHeader = "C://Users/"}

INPATHBM = paste0(INPATHBMHeader,INPATHBM)

DATPATHBM <- file.path(INPATHBM,"data")

data.fileBM <- file.path(DATPATHBM, paste0(expidBM,".csv"))
layout.fileBM <- file.path(DATPATHBM, paste0(expidBM,"_layout.csv"))

cpg <- 0.474 # g/g, gram carbon per gram biomass, Folsom&Carlson 2015
mpc <- 1/12 # mol/g, mol carbon per gram carbon
gps <- 0.603 # mg/ml/scatter, biomass per scatter

skip.wells=sapply(1:8,function(x){paste0(LETTERS[1:6],x)})
skip.wells=as.vector(skip.wells)[-44:-48]


#only reference wells
BMdat <- readExperiment(data.fileBM, type="BioLectorPro", time.conversion=1/3600,
                        layout=layout.fileBM, 
                        fields=c("strain","glc","ace","aTc"), 
                        afields=c("glc","ace","aTc"),
                        blank.id="blank",#blank.data=c("Biomass","Riboflavine"), blank.type = "mean",
                        skip.wells = skip.wells,
                        group1="glc", group2 = c("glc","amount"),
                        group2.color="color")

BMdat = correctBlanks(BMdat, BMdat$layout, yids = c("Biomass","Riboflavine"), type = "mean")


BMdat =  prettyData(BMdat, yids=c(ribof="Riboflavine",O2="DO(Pst3)",
                                scatter="Biomass", pH="pH(HP8)",
                                NADH="NADH - NADPH"))

#all wells
BMdatFull <- readExperiment(data.fileBM, type="BioLectorPro", time.conversion=1/3600,
                        layout=layout.fileBM, 
                        fields=c("strain","glc","ace","aTc"), 
                        afields=c("glc","ace","aTc"),
                        blank.id="blank",#blank.data=c("Biomass","Riboflavine"), blank.type = "mean",
                        group1="glc", group2 = c("glc","amount"),
                        group2.color="color")

BMdatFull = correctBlanks(BMdatFull, BMdatFull$layout, yids = c("Biomass","Riboflavine"), type = "mean")


BMdatFull =  prettyData(BMdatFull, yids=c(ribof="Riboflavine",O2="DO(Pst3)",
                                  scatter="Biomass", pH="pH(HP8)",
                                  NADH="NADH - NADPH"))


viewGroups(BMdat, "ribof")

viewGroups(BMdat, "scatter")

BMdat <- cutData(BMdat, xrng=c(0,25))
BMdatFull <- cutData(BMdatFull, xrng=c(0,25))
nTime = length(BMdat$Time)
concs = BMdat$layout[["amount"]][order(BMdat$layout$well)]
times = BMdat$Time


#read biomass measurement
BMass = read.csv(paste0(INPATHBM,"/biomass030320.csv"))

BMass = list(data = BMass[,-1], time = BMass[,1], description = "Biomass of usually 5ml of Ecoli Z1 in M9 + 0.36% Glucose. Timepoints are Scater 0.5 + X hours")

timeBool = getData(BMdat,"Events")=="PAUSE MEASUREMENT"
BMass$time = BMdat$Events$time[timeBool]

BMass$vol = c(5,5.09,5,5,5,5,5)
BMass$groups = sub("(.*)_[0-9]+", "\\1", colnames(BMass$data))

BMass$normdat = BMass$data/BMass$vol #filterweight [mg/ml] (no meaning)

BMassNormDf = data.frame(BMass$normdat)
BMassNormDf$time = BMass$time

BMassNormDf = melt(BMassNormDf, id.vars = "time")
BMassNormDf$variable = sub("(.*)_[0-9]+", "\\1", BMassNormDf$variable)

#overview plot
ggplot(BMassNormDf)+
  geom_point(aes(x=time, y=value, group=variable, color=variable))+
  labs(x="time [h]", y="filterweight [mg]", color="measurement")+
  #scale_x_continuous(expand=c(0,0))+
  #geom_hline(yintercept = 0, linetype ="dashed")+
  theme_classic()



#compare scatter and riboflavine of sample wells with reference wells
extractTimes = rep(c(BMass$time, Inf), 6) #when were the columns sampled, last column wasn't (excenpt for the firt well)
extractTimes[8] = BMass$time[7]

for(compareVar in c("scatter", "ribof")){
  varDat = getData(BMdatFull, compareVar)
  
  for(i in 1:ncol(varDat)) varDat[c(BMdatFull$Time[-1],Inf) > extractTimes[i],i] = NA #exclude all timepoints after sampeling + one prior, to exclude datapoints, where sampling took place during measurment
  
  sampleIDs = as.numeric(factor(extractTimes))
  sampleIDs[sampleIDs==8] = "not sampled"
  
  pltDf = ggplotDf(varDat, sampleIDs, times)
  pltDf = pltDf[complete.cases(pltDf),]
  
  pltPlt = ggplot(pltDf)+
    geom_line(aes(x=times, y=value, group=as.factor(labels)))+
    geom_hline(yintercept = 0, linetype ="dashed")+
    geom_vline(xintercept=BMass$time, alpha=0.2)+
    geom_rect(data = data.frame(1), fill= "royalblue3", xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf, alpha = 0.2)+
    scale_x_continuous(expand=c(0,0))+
    labs(x="time [h]", y=sprintf("%s [AU]", compareVar), color="sample number")+
    theme_classic()
  
  plts[["glcBM"]][[paste0(compareVar,"CompareGG")]] = pltPlt
}



#subtract pre-weight
BMass$diffdat = apply(BMass$normdat[BMass$groups == "pre"], 1, mean)
BMass$diffdat = BMass$normdat - BMass$diffdat

BMass$Cdat = BMass$diffdat*cpg*mpc * 1000 #mmol/ml


BMassCDf = data.frame(BMass$Cdat[BMass$groups != "pre"])
BMassCDf$time = BMass$time

BMassCDf = melt(BMassCDf, id.vars = "time")
BMassCDf$variable = sub("(.*)_[0-9]+", "\\1", BMassCDf$variable)

measureConf = c("after 1d", "after 2d")
names(measureConf) = c("post1", "post2")
BMassCDf$variable = measureConf[BMassCDf$variable]

ggplot(BMassCDf)+
  geom_point(aes(x=time, y=value, group=variable, color=variable))+
  labs(x="time [h]", y="Biomass [C-mmol/ml]", color="measurement")+
  #scale_x_continuous(expand=c(0,0))+
  geom_hline(yintercept = 0, linetype ="dashed")+
  theme_classic()

BMassGG = geom_point(data=BMassCDf, aes(x=time, y=value, group=variable, color=variable))




#compare biomass measurements

#scatter
biomassScat <- getData(BMdat, "scatter") * 1e3 * gps * cpg * mpc # mmol/L
BMdat <- addData(BMdat, "biomassScat", dat=biomassScat)


#riboflavine
ribof = getData(BMdat, "ribof")
scatter = getData(BMdat, "scatter")

#find stable and linear ("saturated") areas near the end of measurement for scatter and riboflavin for each group of C-molarity
SatBeginnings = data.frame(
  ribof = findSatGrouped(BMdat, "ribof", "amount"),
  scatter = findSatGrouped(BMdat, "scatter", "amount")
)

SatBeginnings = apply(SatBeginnings, 1, max) #take the later saturation for each group

#Each group of Concentrations evaluated at smaller saturation range (ribof or scatter)
#SatBeginnings = rep(SatBeginnings,each = 3)[order(BMdat$layout$well)]

#Alternatively: All groups evaluated in smallest saturation range
SatBeginnings = rep(max(SatBeginnings), 3*length(SatBeginnings))[order(BMdat$layout$well)]


ribof_SatVec = lapply(1:length(SatBeginnings), function(x){
  if(is.na(SatBeginnings[x])){return(NULL)}
  
  ribof[SatBeginnings[x]:nTime,x]
})
ribof_SatVec=unlist(ribof_SatVec)

scatter_SatVec = lapply(1:length(SatBeginnings), function(x){
  if(is.na(SatBeginnings[x])){return(NULL)}
  
  scatter[SatBeginnings[x]:nTime,x]
})
scatter_SatVec=unlist(scatter_SatVec)

ribofScatterLm = lm(ribof_SatVec~0+scatter_SatVec)

ribofScatterScale = coef(ribofScatterLm)
ribofScatterRsq = summary(ribofScatterLm)$r.squared

ribofScatterDf = data.frame(ribof=ribof_SatVec, scatter =scatter_SatVec)


ribofScatterPlt_label = c(sprintf("slope: %.3f", ribofScatterScale), sprintf(": %.3f",ribofScatterRsq))
ribofScatterPlt_label = list(bquote(atop(.(ribofScatterPlt_label[1]),~~~~R^2*.(ribofScatterPlt_label[2]))))


ribofScatterPlt = ggplot()+
  geom_point(data=ribofScatterDf, aes(x=scatter, y=ribof))+
  geom_point(aes(x=0,y=0), color="red", size=2)+
  geom_abline(slope=ribofScatterScale, color="red", linetype="dashed")+
  labs(x="scatter [AU]", y="riboflavine [AU]")+
  geom_label(label=ribofScatterPlt_label, aes(x=max(scatter_SatVec),y=0.5*max(ribof_SatVec)), hjust="right", color="red", parse=T)+
  theme_classic()

plts[["glcBM"]][["ribof_scatterGG"]]= ribofScatterPlt

biomassRibof <- ribof/ ribofScatterScale * 1e3 * gps * cpg * mpc   # C-mmol/L, scale riboflavin to scater and convert to biomass with scatter scaling
BMdat <- addData(BMdat, "biomassRibof", dat=biomassRibof)


#plot biomasses
biomassScatDf = ggplotDf(biomassScat, concs, times)
biomassScatDf = biomassScatDf[complete.cases(biomassScatDf),]
biomassScatPlt = ggplot(biomassScatDf[!biomassScatDf$groups==0 &!is.na(biomassScatDf$groups),])+
  geom_line(aes(x=times, y=value, group=as.factor(labels)))+
  geom_rect(data = data.frame(1), fill= "royalblue3", xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf, alpha = 0.2)+
  geom_hline(yintercept = 0, linetype ="dashed")+
  BMassGG+
  scale_x_continuous(expand=c(0,0))+
  labs(x="time [h]", y="biomass from scatter [C-mmol/l]", group="measurement", color="measurement")+
  theme_classic()
plts[["glcBM"]][["biomassScatGG"]] = biomassScatPlt


biomassRibofDf = ggplotDf(biomassRibof, concs, times)
biomassRibofDf = biomassRibofDf[complete.cases(biomassRibofDf),]
biomassRibofPlt = ggplot(biomassRibofDf[!biomassRibofDf$groups==0 &!is.na(biomassRibofDf$groups),])+
  geom_line(aes(x=times, y=value, group=as.factor(labels)))+
  geom_rect(data = data.frame(1), fill= "royalblue3", xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf, alpha = 0.2)+
  geom_hline(yintercept = 0, linetype ="dashed")+
  scale_x_continuous(expand=c(0,0))+
  labs(x="time [h]", y="biomass from riboflavine [C-mmol/l]", group="measurement", color="measurement")+
  BMassGG+
  theme_classic()
plts[["glcBM"]][["biomassRibofGG"]] = biomassRibofPlt