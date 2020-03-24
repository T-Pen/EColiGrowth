PC = T

#======================================= functions =======================================
ode_monodMu = function (time, init, parms, ...) {
  with(as.list(c(parms, init)), {
    mu <- mumax * s/(s + K)
    ds <- phi * (sin - s) - (mu/Y) * y
    dy <- (mu - phi) * y - d * y
    list(c(ds, dy, mu))
  })
}

ode_monodT = function (time, init, parms, ...){
  with(as.list(c(parms, init)), {
    mu_t = mumax_t * s/(s + Kst)
    mu <- mumax * t/(t + Kty) #*Regulation !!! e.g. (1 - s/(s + Ksty))
    
    ds <- phi * (sin - s) - mu_t * y
    dt = phi * (tin - t) + (n_st * mu_t - (mu/(Y * n_st))) * y
    dy <- (mu - phi) * y - d * y
    list(c(ds, dy))
  })
}

ode_anacatT = function (time, init, parms, ...){
  with(as.list(c(parms, init)), {
    adp <- axp - atp
    mu_t <- mumax_t * s/(s + Kst) * adp/(adp + Kat)
    mu_ab <- mumax_ab * t/(t + Ktab) * atp/(atp + Kaab)
    mu_cd <- mumax_cd * t/(t + Ktcd) * adp/(adp + Kacd) #s Regulation?
    mu_m <- mumax_m
    args <- list(...)
    if (!is.null(args[["funcd"]])) {
      for (i in 1:length(funcd)) assign(names(funcd)[i], 
                                        do.call(funcd[[i]], as.list(c(init, parms))))
    }
    
    ds <- phi * (sin - s) - mu_t * y
    dt <- phi * (tin - t) + n_st * mu_t * y - (mu_ab + mu_cd) * y
    dy <- (mu_ab - phi) * y
    datp = (n_cd * mu_cd - n_ab * mu_ab - mu_m) * Cc/Vc + 
      n_t * mu_t - mu_ab * atp
    fnctdefs <- NULL
    if (!is.null(args[["funcd"]])) {
      fnctdefs <- rep(NA, length(funcd))
      names(fnctdefs) <- names(funcd)[i]
      for (i in 1:length(funcd)) fnctdefs[i] <- get(names(funcd)[i])
    }
    list(c(ds, dt, dy, datp), adp = adp, mu = mu_ab, mu_cd = mu_cd, mu_t = mu_t, 
         fnctdefs)
  })
}

summarizeModel = function(datMat, model, mainVar, datGroups, datTimes, modelGroups, modelTimes, datType, modelType, units=list(), mods, relPlts = list(names(model)), replaceVars = c()){
  res = list()
  
  relVarGG = list()
  
  modelMat = model[[mainVar]]
  
  groupMatch = match(datGroups, modelGroups)
  timesMatch = match(datTimes, modelTimes)
  
  datDf = ggplotDf(datMat, datGroups, datTimes)
  datDf = datDf[complete.cases(datDf),]
  datPlt = ggplot(datDf)+
    geom_line(aes(x=times, y=value, group=as.factor(labels)))+
    facet_grid(groups~.)+
    geom_hline(yintercept = 0, linetype ="dashed")+
    scale_x_continuous(expand=c(0,0))+
    labs(x="time [h]", y=sprintf("%s [%s]", datType, units[[mainVar]]))+
    mods[["bgColsGG"]]+
    mods[["timeStartGG"]]+
    mods[["labelFacetGG"]]+
    mods[["scatterPeaksGG"]]+
    theme_classic()
  
  #fit overlay
  modelDf = ggplotDf(modelMat, modelGroups, modelTimes)
  modelDf = modelDf[complete.cases(modelDf),]
  
  modelGG = geom_line(data = modelDf, aes(x=times, y=value, group=as.factor(labels), color=modelType), size=1, )
  
  modelFitPlt = datPlt+
    modelGG+
    scale_color_manual(name = "fitted curves", values="purple")
  res[["fit"]] = modelFitPlt
  
  
  #residuals
  modelResids = modelMat[timesMatch,groupMatch] - datMat
  
  modelResidsDf = ggplotDf(modelResids, datGroups, datTimes)
  modelResidsDf = modelResidsDf[complete.cases(modelResidsDf),]
  
  modelResidsPlt = ggplot(modelResidsDf)+
    geom_point(aes(x=times, y=value))+
    facet_grid(groups~.)+
    geom_hline(yintercept = 0, linetype ="dashed")+
    scale_x_continuous(expand=c(0,0))+
    labs(x="time [h]", y=sprintf("model redsuduals [%s]", units[[mainVar]]))+
    mods[["bgColsGG"]]+
    mods[["timeStartGG"]]+
    mods[["labelFacetGG"]]+
    mods[["scatterPeaksGG"]]+
    theme_classic()
  res[["residuals"]] = modelResidsPlt
  
  #model MSE
  modelMSE = colMeans(modelResids^2)
  
  modelMSEDf = data.frame(MSE = modelMSE, groups = datGroups)
  modelMSEDf = modelMSEDf[complete.cases(modelMSEDf),]
  
  bgColsBox = (modelGroups[-1] - modelGroups[-length(modelGroups)])/2
  bgColsBox = bgColsBox[c(1:length(bgColsBox), length(bgColsBox))]
  bgColsBoxDf = data.frame(xmin=modelGroups - bgColsBox, xmax=modelGroups+bgColsBox, col= mods[["bgColsGG"]]$data$col)
  
  modelMSEPlt = ggplot(modelMSEDf)+
    geom_boxplot(aes(x=groups, y=MSE, group=groups))+
    geom_hline(yintercept = 0, linetype ="dashed")+
    labs(x="glucose concentration [C-mmol/l]", y=bquote(.(sprintf("model MSE [(%s)", units[mainVar]))^2~"]"))+
    geom_rect(data = bgColsBoxDf, aes(xmin =xmin, xmax=xmax, fill = col), fill= bgColsBoxDf$col ,ymin = -Inf,ymax = Inf, alpha = 0.25)+
    scale_x_continuous(expand=c(0,0), breaks = modelGroups)+
    theme_classic()
  res[["mse"]] = modelMSEPlt
  
  #model variables
  for(nam in names(model)){
    varMat = model[[nam]]
    
    namVis = if(nam %in% names(replaceVars)){replaceVars[nam]}else{nam}
    
    varDf = ggplotDf(varMat, modelGroups, modelTimes)
    varDf = varDf[complete.cases(varDf),]
    varPlt = ggplot(varDf)+
      geom_line(aes(x=times, y=value, group=as.factor(labels)))+
      facet_grid(groups~.)+
      geom_hline(yintercept = 0, linetype ="dashed")+
      scale_x_continuous(expand=c(0,0))+
      labs(x="time [h]", y=sprintf("%s [%s]", namVis, units[[nam]]))+
      mods[["bgColsGG"]]+
      mods[["labelFacetGG"]]+
      theme_classic()
    
    res[[nam]] = varPlt
    
    #data for relative plots
    varMax = suppressWarnings(apply(varMat, 2, max, na.rm=T))
    
    varMax[!complete.cases(t(varMat))] = NA
    
    varMat_rel = t(t(varMat)/varMax) * 100
    
    varDf_rel = ggplotDf(varMat_rel, modelGroups, modelTimes)
    varDf_rel = varDf_rel[complete.cases(varDf_rel),]
    varDf_rel$col = namVis
    
    relVarGG[[nam]] = geom_line(data = varDf_rel, aes(x=times, y=value, group=as.factor(labels), color = col))
  }
  
  #relative plots
  for(setNum in 1:length(relPlts)){
    varSet = relPlts[[setNum]]
    
    varSetPlt = ggplot()+
      facet_grid(groups~.)+
      geom_hline(yintercept = 0, linetype ="dashed")+
      scale_x_continuous(expand=c(0,0))+
      labs(x="time [h]", y="percent max value of model variable [%]")+
      mods[["bgColsGG"]]+
      mods[["labelFacetGG"]]+
      scale_color_discrete(name = "model variable")+
      theme_classic()
      
      for(varNam in varSet) varSetPlt = varSetPlt + relVarGG[[varNam]]
    
    res[[paste0("relVar", setNum)]] = varSetPlt
  }
  
  
  return(res)
}

mat2Latex = function(mat, fstring, prnt = T){
  fstring = gsub(" ", " & ", fstring)
  fstring = sprintf("%s \\", fstring)
  
  res = apply(mat, 1, function(x){do.call(sprintf, c(list(fstring), x))})
  res = matrix(res)
  
  headr = paste(colnames(mat) ,collapse = " & ")
  headr = sprintf("%s \\", headr)
  
  res = rbind(headr, res)
  dimnames(res) = NULL
  
  if(prnt){
    apply(res, 1, function(x){
      cat(sprintf("%s\\ \n", x))
    })
  }
  return(res)
}

#==========================================================================================


if(PC){INPATHHeader = "D://"}else{INPATHHeader = "C://Users/"}

load(paste0(INPATHHeader,"/tpeng/OneDrive/OD__Universität/BA__Axmann/biolector_evalDAT.RData"))

dat = dats$glc

dfg <- data2growthrates(dat, "biomass", plate=dat$layout)

#concentrations
concs = dat$layout[["amount"]][order(dat$layout$well)]
concsUniq = sort(unique(dat$layout$amount))
times=dat$Time
yield = coefs$glc$yield
scatterPeakThresh = -0.15
scatterPeakMinL = 5

bgColsGG = plts$glc$mods$bgColsGG
timeStartGG = plts$glc$mods$timeStartGG
labelFacetGG = plts$glc$mods$labelFacetGG
scatterPeaksGG = plts$glc$mods$scatterPeaksGG



# #find Peaks and their times
# scatterPeaksAllDf = 
#   findDrops(mat = getData(dat, "scatterSM"),
#             times = dat$Time,
#             groups = concs,
#             slopeThresh = scatterPeakThresh,
#             minl = scatterPeakMinL,
#             timeThresh = 1)
# 
# scatterPeaksAllDf$well = dat$layout$well[order(dat$layout$well)][scatterPeaksAllDf$column]
# scatterPeaksAllDf$row = sub("([A-Z]+)([1-9])", "\\1", scatterPeaksAllDf$well)
# scatterPeaksAllDf$col = sub("([A-Z]+)([1-9])", "\\2", scatterPeaksAllDf$well)
# 
# scatterPeaksAllDf_appear = table(scatterPeaksAllDf$column)
# 
# #show peaks in wells
# scatter = getData(dat,"scatter")
# 
# scatterPlateDf = ggplotDf(scatter, concs, times)
# scatterPlateDf = scatterPlateDf[complete.cases(scatterPlateDf),]
# scatterPlateDf$row = sub("([A-Z]+)([1-9])", "\\1", scatterPlateDf$labels)
# scatterPlateDf$col = sub("([A-Z]+)([1-9])", "\\2", scatterPlateDf$labels)
# 
# bgColsPlate = data.frame(row=c("A","B","C"),
#                          col=rep(1:7, each=3),
#                          cols=rep(bgColsGG$data$col, each=3))
# 
# scatterPlatePlt = ggplot(scatterPlateDf)+
#   geom_line(aes(x=times, y=value, group=as.factor(labels)))+
#   facet_grid(row~col)+
#   scale_x_continuous(expand=c(0,0))+
#   labs(x="time [h]", y="scatter [AU]")+
#   geom_rect(data = bgColsPlate, aes(fill = cols), fill= bgColsPlate$cols, xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf, alpha = 0.25)+
#   geom_rect(data = scatterPeaksAllDf, aes(xmin = t1, xmax = t2), fill= "grey", ymin = -Inf,ymax = Inf, alpha = 0.4)+
#   theme_classic()
# plts$glc[["scatterPeaksPlateGG"]] = scatterPlatePlt
# 
# 
# #only take the peaks beginning times
# scatterPeaksPos = matrix(NA,length(scatterPeaksAllDf_appear),max(scatterPeaksAllDf_appear))
# row.names(scatterPeaksPos) = names(scatterPeaksAllDf_appear)
# 
# for(k in names(scatterPeaksAllDf_appear)){
#   columnPos = which(scatterPeaksAllDf$column == k)
#   scatterPeaksPos[k,1:length(columnPos)] = scatterPeaksAllDf[columnPos, "t1"]
# }
# 
# scatterPeaksPos = data.frame(scatterPeaksPos)
# colnames(scatterPeaksPos) = sub("X", "Peak",colnames(scatterPeaksPos))
# 
# #also mark differences in peak times per well
# scatterPeaksPos_difs = sapply(2:ncol(scatterPeaksPos), function(colNum){
#   scatterPeaksPos[,colNum] - scatterPeaksPos[,colNum-1]
# })
# 
# #colnames(scatterPeaksPos_difs) = sapply(2:ncol(scatterPeaksPos), function(colNum){paste(colnames(scatterPeaksPos)[c(colNum-1,colNum)], collapse="_")})
# colnames(scatterPeaksPos_difs) = paste0("Delta*", colnames(scatterPeaksPos)[-1])
# 
# scatterPeaksPos = cbind(scatterPeaksPos, scatterPeaksPos_difs)
# 
# #add the respectibe glucose level
# scatterPeaksPos$groups = sapply(names(scatterPeaksAllDf_appear), function(x){
#   scatterPeaksAllDf[scatterPeaksAllDf$column == x, "groups"][1]
# })
# 
# scatterPeaksPosDf = melt(scatterPeaksPos, id.vars = "groups")
# scatterPeaksPosDf = scatterPeaksPosDf[complete.cases(scatterPeaksPosDf),]
# 
# 
# bgColsBox = (concsUniq[-1] - concsUniq[-length(concsUniq)])/2
# bgColsBox = bgColsBox[c(1:length(bgColsBox), length(bgColsBox))]
# bgColsBoxDf = data.frame(xmin=concsUniq - bgColsBox, xmax=concsUniq+bgColsBox, cols= bgColsGG$data$col)
# 
# scatterPeaksPosPlt = ggplot()+
#   geom_point(data=scatterPeaksPosDf, aes(x=groups, y=value, group=groups))+
#   facet_grid(variable~., labeller = label_parsed)+
#   geom_rect(data = bgColsBoxDf, aes(xmin =xmin, xmax=xmax, fill = cols), fill= rep(bgColsBoxDf$cols, 3),ymin = -Inf,ymax = Inf, alpha = 0.25)+
#   geom_hline(yintercept = 0, linetype ="dashed")+
#   labs(y="time of respective peak [h]", x="glucose concentration [C-mmol/l]")+
#   scale_x_continuous(breaks=concsUniq, expand=c(0,0))+
#   theme_classic()
# plts$glc[["scatterPeaksPosGG"]] = scatterPeaksPosPlt


#======================================= monod model =======================================

# mu <- mumax * s/(s + K)
# ds <- phi * (sin - s) - (mu/Y) * y
# dy <- (mu - phi) * y - d * y


#average the biomass per group for finding the parameters
biomassMean = sapply(concsUniq, function(x){
  rowMeans(getData(dat, "biomass")[,concs == x], na.rm = T)
})
colnames(biomassMean) = concsUniq

#average the growthrates per group for finding K
growthrateMean = sapply(concsUniq, function(x){
  rowMeans(getData(dat, "biomassSM.mu.dpseg")[,concs == x], na.rm = T)
})
colnames(growthrateMean) = concsUniq

#find Saturations for determining the final Biomass
sats = findSatGrouped(dat, "biomass", "amount")

ntimes=length(dat$Time)


#######################
odeTime <- inflateTimes(times, 3)
odeCoefsMintime = 3 #h


modelBM = matrix(nrow= length(odeTime), ncol = length(concsUniq))
modelBM = list(s=modelBM, y=modelBM, mu=modelBM)

odeKs = c()
odeCoefs = matrix(nrow=length(concsUniq), ncol=8, dimnames = list(concsUniq, c("s", "y", "phi", "sin", "mumax", "K", "Y", "d")))

odeCoefsMintime = which(times > odeCoefsMintime)[1]



for(concID in 2:length(concsUniq)){
startTime = 0

conc = concsUniq[concID]
initBM = biomassMean[1,concID]

maxBM = mean(biomassMean[sats[concID]:ntimes,concID])

maxRate = max(growthrateMean[odeCoefsMintime : length(times),concID], na.rm = T)


# K is the substrate concentration, where mu(max) is half-maximal
KTime = tail(which(growthrateMean[, concID]> maxRate*0.5), 1)
#KTime = times[KTime]

KBM = biomassMean[KTime, concID]

K = KBM/ yield #used up Substrate
K = conc - K

odeKs[concID] = K

names(odeKs)[concID] = as.character(KTime)


#initial values
yinit <- c(s=conc, y=initBM, mu=0) #initial values for s (substrate) y (biomass) (mu, only for output)
poinit <- c(phi=0, sin=0, mumax= maxRate, K = K, Y= yield, d=0) #phi: dilution rate, sin: incoming substrate concentration if diluted, mumax: maximal growthrate, K: Saturation constant treating the growth as an enzymatic reaction 

names(poinit) = c("phi", "sin", "mumax", "K", "Y", "d")

odeCoefs[concID,] = c(yinit[-3], poinit)

odem <- ode(yinit, odeTime, ode_monodMu, poinit)

modelBM[["y"]][,concID] = odem[,"y"]
modelBM[["s"]][,concID] = odem[,"s"]
modelBM[["mu"]][,concID] = odem[,"mu"]
}

BMMonod=
summarizeModel(datMat = getData(dat, "biomass"),
               model = modelBM,
               mainVar = "y",
               datGroups = concs,
               datTimes = times,
               modelGroups = concsUniq,
               modelTimes = odeTime,
               datType = "biomass",
               modelType = "monod model",
               units = c("y"= "C-mmol/l",
                         "s" = "C-mmol/l",
                         "mu" = "C-mmol/(l h)"),
               mods = plts$glc$mods,
               replaceVars = c("mu" = "\u03BC"))
plts$glcMdl[["monod"]] = BMMonod

#mat2Latex(odeCoefs, "%.1f %.3f %d %d %.3f %.2f %.3f %d")

# #======================================= modified monod model =======================================
# 
# # mu_t = mumax_t * s/(s + Kst)
# # mu <- mumax * t/(t + Kty) * (s/(s + Rs))
# # ds <- phi * (sin - s) - mu_t * y
# # dt = phi * (tin - t) + (n_st * mu_t - (mu/(Y * n_st))) * y
# # dy <- (mu - phi) * y - d * y
# 
# #for useful model (n_st * mu_t) > (mu/(Y * n_st))
# 
# yinit <- c(s=conc, t=0, y=initBM) #initial values for s (substrate), t (transfer metabolite), y (biomass)
# poinit <- c(phi=0,
#             sin=0,
#             tin = 0,
#             mumax_t = NA, #max glucose uptake rate
#             mumax = maxRate, #max growth rate
#             Kst = NA, #(Glucose-K_M of glucose transporter?)
#             Kty = NA, #(if glycolysis if fast then Kty = K, else: ???)
#             n_st = NA, #number of transfer molecules per glucose
#             Y= yield, #yield
#             d=0)
# time <- seq(0,24,.1)
# 
# 
# #======================================= modified anacat model =======================================
# 
# # adp <- axp - atp
# # mu_t <- mumax_t * s/(s + Kst) * adp/(adp + Kat)
# # mu_ab <- mumax_ab * t/(t + Ktab) * atp/(atp + Kaab)
# # mu_cd <- mumax_cd * t/(t + Ktcd) * adp/(adp + Kacd)
# # mu_m <- mumax_m
# # ds <- phi * (sin - s) - mu_t * y
# # dt <- phi * (tin - t) + n_st * mu_t * y - (mu_ab + mu_cd) * y
# # dy <- (mu_ab - phi) * y
# # datp = (n_cd * mu_cd - n_ab * mu_ab - mu_m) * Cc/Vc + 
# #   n_t * mu_t - mu_ab * atp
# 
# yinit <- c(s=conc, t=0, y=initBM) #initial values for s (substrate), t (transfer metabolite), y (biomass)
# poinit <- c(phi=0,
#             sin=0,
#             tin=0,
#             mumax_t=NA, #max glucose uptake rate
#             mumax_ab=NA, #max growth rate
#             mumax_cd=NA, #(max acetate uptake rate?)
#             mumax_m=NA, #(max glucose uptake rate in steady state)
#             Kst=NA, #(Glucose-K_M of glucose transporter?)
#             Kat=NA, #(Smallest ATP-K_M in glycolysis: pyruvatkinase / phosphoglyceratkinase)
#             Ktab=NA, #
#             Kaab=NA,
#             Ktcd=NA,
#             Kacd=NA,
#             n_st=NA,
#             n_cd=NA,
#             n_ab=NA,
#             n_t=NA,
#             Cc=NA,
#             Vc=NA)
# time <- seq(0,24,.1)
# 
# ## simulate ode from starting parameters
# odem <- ode(yinit, time, ode_monod, poinit)
# plot(odem[,"time"], odem[,"y"],col=2, log="y")



savePlts(plts)

savePltsPDF(pltList = plts,
            pdfName = "summaryPlts",
            pltVec = findPlts("../documents/summaryBioLector.txt"),
            filePath = "../documents/")
