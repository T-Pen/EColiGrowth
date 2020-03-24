library(platexpress)
library(reshape2)
library(pracma)
library(ggplot2)
library(mugro)
library(growthrates)
library(pspline)
library(tidyr)

PC = T

INPATH <- "tpeng/OneDrive/OD__Universität/BA__Axmann"
expid <- "25_Ecoli_2020_REDUCTION-1"

MaxConc = c("glc" = 22.2*6,
            "ace" = 66.6*2)

biomassEqivalent = "ribof"

statVar = "scatter"

settleTime = 0.5



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
    
    evals_var = evals_var/(endVal - evals) #reduce accounted variance by number of included points to prefer larger ranges
    
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

findSatGrouped = function(dat, ID, group, npoints = 10, timeOnCol = T){
  satDat = getData(dat, ID)
  
  grps = dat$layout[,group]
  
  blnk = dat$layout$strain == "blank"
  
  grps_uniq = unique(grps)
  grps_uniq = grps_uniq[!is.na(grps_uniq)]
  
  
  sapply(grps_uniq, function(i){
    
    
    subCols = (grps == i) & !blnk
    subCols = subCols[order(dat$layout$well)]
    subCols[is.na(subCols)] = F
    
    subDat = satDat[,subCols]
    
    findSat(subDat,10,F)
  })
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

savePlts = function(pltList, filePath = ""){
  cats = names(pltList)
  
  for(category in cats){
    nams = names(pltList[[category]])
    nams = nams[nams != "mods"]
    
    mdlBool = grepl("Mdl$",category)
    
    for(nam in nams){
      
      if(!mdlBool){
        pltPath = paste0(filePath, category, "_", nam)
        
        svg(paste0(pltPath, ".svg"))
        print(pltList[[category]][[nam]])
        dev.off()
        
        png(paste0(pltPath, ".png"))
        print(pltList[[category]][[nam]])
        dev.off()
      
      }else{
        subnams = names(pltList[[category]][[nam]])
        
        for(subnam in subnams){
          pltPath = paste0(filePath, category, "_", nam, toupper(sub("^(.).*", "\\1", subnam)), sub("^.(.*)", "\\1", subnam))
        
          svg(paste0(pltPath, ".svg"))
          print(pltList[[category]][[nam]][[subnam]])
          dev.off()
          
          png(paste0(pltPath, ".png"))
          print(pltList[[category]][[nam]][[subnam]])
          dev.off()
        }
      }
    }
  }
}

savePltsSingle = function(pltList, pltVec, filePath = ""){
  for(pltNam in pltVec){
    
    pltNamSplt = strsplit(pltNam, "_")[[1]]
    category = pltNamSplt[1]
    nam = pltNamSplt[2]
    
    mdlBool = grepl("Mdl$",category)
    pltPath = paste0(filePath, pltNam)
      
    if(!mdlBool){
      svg(paste0(pltPath, ".svg"))
      print(pltList[[category]][[nam]])
      dev.off()
      
      png(paste0(pltPath, ".png"))
      print(pltList[[category]][[nam]])
      dev.off()
      
    }else{
      subnam = regexpr("[A-Z]", nam)
      subnam = substring(nam,subnam)
      subnam = paste0(tolower(sub("^(.).*", "\\1", subnam)), sub("^.(.*)", "\\1", subnam))
      
      nam = substring(nam,1,subnam-1)
      
      svg(paste0(pltPath, ".svg"))
      print(pltList[[category]][[nam]][[subnam]])
      dev.off()
      
      png(paste0(pltPath, ".png"))
      print(pltList[[category]][[nam]][[subnam]])
      dev.off()
    }
  }
}

savePltsPDF = function(pltList, pltVec, pdfName, filePath = ""){
  pdf(paste0(filePath, pdfName, ".pdf"))
  
  for(pltNam in pltVec){
    
    pltNamSplt = strsplit(pltNam, "_")[[1]]
    category = pltNamSplt[1]
    nam = pltNamSplt[2]
    
    mdlBool = grepl("Mdl$",category)
    
    if(!mdlBool){
      print(pltList[[category]][[nam]])
    }else{
      subnamPos = regexpr("[A-Z]", nam)
      subnam = substring(nam,subnamPos)
      subnam = paste0(tolower(sub("^(.).*", "\\1", subnam)), sub("^.(.*)", "\\1", subnam))
      
      nam = substring(nam,1,subnamPos-1)
      
      print(pltList[[category]][[nam]][[subnam]])
    }
  }
  dev.off()
}

findPlts = function(filepath){
  latext = readLines(filepath)
  picLines = grep("includesvg",latext)
  
  pics = sub("^.*\\{(.*)\\}.*$", "\\1", latext[picLines])
  
  return(pics)
}

excludedScaling = function(df, excltime=NULL, exclgroups=c()){
  exclBool = logical(nrow(df))
  
  if(!is.null(excltime)) exclBool = exclBool | df$times < excltime
  
  if(!length(exclgroups) == 0) exclBool = exclBool | (df$groups %in% exclgroups)
  
  newRange = range(df[!exclBool, "value"])
  
  newRangeexpand = (newRange[2]-newRange[1]) * 0.05
  
  newRange = newRange + c(-1, 1) * newRangeexpand
  
  return(newRange)
}

smootDat = function(varDat, time){
  res = list()
  
  varSM = apply(varDat, 2, function(column){
    sm.spline(x=time, y=column)}
  )
  
  res$smthObj = varSM #save smoothing object
  
  varSM = sapply(varSM, function(x){
    x$ysmth
  })
  
  row.names(varSM) = row.names(varDat)
  
  res$smthDat = varSM
  
  return(res)
}

inflateTimes = function(times, fac){
  timediff = times[-1] - times[-length(times)]
  
  newTimeDiff = rep(timediff/fac, each = fac)
  
  newTime = c(times[1], newTimeDiff)
  
  newTime = sapply(2:length(newTime), function(x){
    newTime[x] <<-newTime[x]+newTime[x-1]
  })
  
  newTime[c(rep(F, fac-1), T)] = times[-1]
  
  newTime = c(times[1], newTime)
  
  return(newTime)
}

findDrops = function(mat, times, groups, slopeThresh, trans = "lin", minl = 3, timeThresh = 1){
  matPeakDpseg = apply(mat, 2, function(matCol){
    
    if(trans == "lin"){
      y=matCol
    }else if(trans == "log"){
      y=log(scatCol)
    }else{simpleError(message = "unknown value for 'trans'")}
    
    dpseg::dpseg(x=times, y=y, verb=F, type = "r2", minl = minl)
  })
  
  
  matPeaks = lapply(matPeakDpseg, function(matDpseg){
    
    #ignore anything before t= timeThresh
    matDpseg = matDpseg$segments[matDpseg$segments$x2>timeThresh,]
    
    #arbtrary threshold to exclude small slopes
    matDpseg = matDpseg[matDpseg$slope < slopeThresh,c("x1", "x2", "start", "end")]
    
    
    if(nrow(matDpseg) == 0) return(NULL)
    
    
    start = matDpseg$start
    end = matDpseg$end
    
    #find the rows, which are not following to another one (those are the beginning of a new range)
    bool = start[-1] != end[-length(end)]
    bool = c(T, bool)
    
    #rows which end a range
    boolEnd = c(bool[-1], T)
    
    
    if(length(bool) == 1) if(!bool) return(NULL)
    
    t1 = matDpseg$x1[bool]
    t2 = matDpseg$x2[boolEnd]
    
    n1 = matDpseg$start[bool]
    n2 = matDpseg$end[boolEnd]
    
    res = data.frame(t1, t2, n1, n2)
    
    return(res)
  })
  
  matPeaksDf = data.frame()
  for(z in 1:length(matPeaks)){
    mat = matPeaks[[z]]
    if(is.null(mat)) next
    
    mat = data.frame(mat)
    
    mat$groups = groups[z]
    
    mat$column = z
    
    matPeaksDf = rbind(matPeaksDf, mat)
  }
  
  return(matPeaksDf)
}

updatePackages = function(packs="platexpress", platexpressSource='local_mod', PC=F){
  
  if(PC){INPATHHeader = "D://"}else{INPATHHeader = "C://Users/"}
  
  ## load devtools
  require(devtools)
  
  ## Thomas Petzoldt's `growthrates` package, development version
  if("growthrates" %in% packs | packs == "all")
    install_github('tpetzoldt/growthrates')
  
  ## in-house package to parse and analyze platereader data
  if("platexpress" %in% packs | packs == "all"){
    if(platexpressSource=='git'){
      devtools::install_github('raim/platexpress')
    }else if(platexpressSource=='local_mod'){
      install.packages(paste0(INPATHHeader,'tpeng/OneDrive/OD__Universität/BA__Axmann/R/platexpressModified'), repos = NULL, type = 'source')
    }else if(platexpressSource=='local'){
      install.packages(paste0(INPATHHeader,'tpeng/OneDrive/OD__Universität/BA__Axmann/git_platexpress'), repos = NULL, type = 'source')
    }
  }
  
  ## in-house package adding growth models for use with `growthrates`
  if("growthmodels" %in% packs | packs == "all")
    devtools::install_git('https://gitlab.com/raim/growthmodels.git', subdir = 'pkg', quiet = FALSE)
  
  ## in-house package to segment growth curves into linear pieces
  if("growthphases" %in% packs | packs == "all")
    devtools::install_git('https://gitlab.com/raim/growthphases.git', subdir = 'pkg', quiet = FALSE)
  
}


#=========================================================================================

if(PC){INPATHHeader = "D://"}else{INPATHHeader = "C://Users/"}

INPATH = paste0(INPATHHeader,INPATH)

DATPATH <- file.path(INPATH,"data")
RESPATH <- file.path(INPATH,"results") # results

setwd(RESPATH)


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




#======================================= analysis =======================================

dats = list()
plts = list("glc"=list(),
            "ace"=list(),
            "glcBM"=list(),
            "glcMdl" = list())
dpsegs = list("glc"=list(),
              "ace"=list())
coefs = list("glc"=list(),
             "ace"=list())
smooths = list("glc"=list(),
               "ace"=list())

#=================================== get biomass data ===================================
source(paste0(INPATHHeader,'tpeng/OneDrive/OD__Universität/BA__Axmann/R/biolector_eval_biomassG.R'), echo=TRUE)



i= c("glucose", "glc")
#i= c("acetate", "ace")

for(i in list(c("glucose", "glc"), c("acetate", "ace"))){
  CSource =i[1]
  CsourceShort = i[2]
  
  if(CsourceShort =="glc"){
    skipWellList = c(paste0(toupper(letters[4:6]), rep(1:7,3)))
    Amounts = "amount"
    PltColors = "color"
    bgCols = colorRampPalette(c("white", "royalblue3"))(7)
    xrng=c(0,20)
    
    scatterPeakThresh = -0.1
    scatterPeakMinL = 5
    
  }else if(CsourceShort =="ace"){
    skipWellList = c(paste0(toupper(letters[1:3]), rep(1:7,3)), "E7") #additionally E7, sinceit shows too large scatter values
    Amounts = "ace.amount"
    PltColors = "ace.color"
    bgCols = colorRampPalette(c("white", "firebrick2"))(7)
    xrng=c(0,50)
    
    #scatterPeakThresh = -0.02
    #scatterPeakMinL = 5 # slopes are long and gentle
    scatterPeakThresh = -0.015
    scatterPeakMinL = 15 # slopes are long and gentle
  }

  #read data
  dat <- readExperiment(data.file, type="BioLectorPro", time.conversion=1/3600,
                        layout=layout.file, 
                        fields=c("strain","glc","ace","aTc"), 
                        afields=c("glc","ace","aTc"),
                        blank.id="blank",blank.data=c("Biomass","Riboflavine"),
                        skip.wells = skipWellList,
                        group1=CsourceShort, group2 = c(CsourceShort,Amounts),
                        group2.color=PltColors)
  
  dat <- prettyData(dat, yids=c(ribof="Riboflavine",O2="DO(Pst3)",
                                scatter="Biomass", pH="pH(HP8)",
                                NADH="NADH - NADPH"))
  
  #format data
  dat$layout[[Amounts]] <- (dat$layout[[Amounts]]/max(dat$layout[[Amounts]])) * MaxConc[CsourceShort] #convert to C-mmol
  dat$layout[[Amounts]][dat$layout$strain=="blank"] <- NA
  dat$layout$color[dat$layout$strain=="blank"] <- "#999999"
  
  
  #cut off long stationary phase
  # xStat = findSatGrouped(dat, "scatter", Amounts, npoints = 10)
  
  dat <- cutData(dat, xrng=xrng)
  nTime = length(dat$Time)
  concs = dat$layout[[Amounts]][order(dat$layout$well)]
  concsUniq = sort(unique(concs))
  times = dat$Time
  concCols = c("#BEBEBE", bgCols[-1]) #for coloring points: 0 C-mmol is grey, not white
  bgColsDf = data.frame(col=bgCols, groups = sort(unique(dat$layout[,Amounts])))
  
  
  #raw plots
  viewGroups(dat, yids="scatter", show.ci95=TRUE, lwd.orig=0)
  plts[[CsourceShort]][["scatter"]]= recordPlot()
  viewGroups(dat, yids="ribof", show.ci95=TRUE, lwd.orig=0)
  plts[[CsourceShort]][["ribof"]]= recordPlot()
  viewGroups(dat, yids="O2", show.ci95=TRUE, lwd.orig=0)
  plts[[CsourceShort]][["O2"]]= recordPlot()
  viewGroups(dat, yids="pH", show.ci95=TRUE, lwd.orig=0)
  plts[[CsourceShort]][["pH"]]= recordPlot()
  viewGroups(dat, yids="NADH", show.ci95=TRUE, lwd.orig=0)
  plts[[CsourceShort]][["NADH"]]= recordPlot()
  
  
  
  
  #ggplots
  
  bgColsGG = geom_rect(data = bgColsDf, aes(fill = col), fill= bgCols, xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf, alpha = 0.25)
  
  timeStartGG = geom_vline(xintercept = settleTime, linetype="dotted")
  
  labelFacetGG = scale_y_continuous(sec.axis = sec_axis(name = sprintf("%s concentration [C-mmol/l]\n", CSource), trans =~.*0, labels=NULL, breaks=NULL))
  
  
  
  # #find peaks via dpseg and mark them in plots
  # scatterMean = sapply(concsUniq, function(x){
  #   rowMeans(getData(dat, "scatter")[,concs == x], na.rm = T)
  # })
  # 
  # #smoothing necessary?
  # scatterMean = smootDat(scatterMean, dat$Time)$smthDat
  # # scatterMean = apply(scatterMean, 2, ma, 3)
  # 
  # scatterPeaksDf = findDrops(mat = scatterMean,
  #                            times = dat$Time,
  #                            groups = concsUniq,
  #                            slopeThresh = scatterPeakThresh,
  #                            minl = scatterPeakMinL,
  #                            timeThresh = 1)
  # 
  # #more
  # 
  # scatterPeaksGG = geom_rect(data = scatterPeaksDf, aes(xmin = t1, xmax = t2), fill= "grey", ymin = -Inf,ymax = Inf, alpha = 0.4)
  
  
  #find Peaks and their times via dpseg
  scatterPeaksAllDf = 
    findDrops(mat = smootDat(getData(dat, "scatter"), dat$Time)$smthDat, #use smoothed scatter for peak-determiantion
              times = dat$Time,
              groups = concs,
              slopeThresh = scatterPeakThresh,
              minl = scatterPeakMinL,
              timeThresh = 1)
  
  scatterPeaksAllDf$well = dat$layout$well[order(dat$layout$well)][scatterPeaksAllDf$column]
  scatterPeaksAllDf$row = sub("([A-Z]+)([1-9])", "\\1", scatterPeaksAllDf$well)
  scatterPeaksAllDf$col = sub("([A-Z]+)([1-9])", "\\2", scatterPeaksAllDf$well)
  
  scatterPeaksAllDf_appear = table(scatterPeaksAllDf$column)
  
  #show peaks in wells
  scatter = getData(dat,"scatter")
  
  scatterPlateDf = ggplotDf(scatter, concs, times)
  scatterPlateDf = scatterPlateDf[complete.cases(scatterPlateDf),]
  scatterPlateDf$row = sub("([A-Z]+)([1-9])", "\\1", scatterPlateDf$labels)
  scatterPlateDf$col = sub("([A-Z]+)([1-9])", "\\2", scatterPlateDf$labels)
  
  bgColsPlate = data.frame(row=if(CsourceShort == "glc"){c("A","B","C")}else{c("D","E","F")},
                           col=rep(1:7, each=3),
                           cols=rep(bgCols, each=3))
  
  scatterPlatePlt = ggplot(scatterPlateDf)+
    geom_line(aes(x=times, y=value, group=as.factor(labels)))+
    facet_grid(row~col)+
    scale_x_continuous(expand=c(0,0))+
    labs(x="time [h]", y="scatter [AU]")+
    geom_rect(data = bgColsPlate, aes(fill = cols), fill= bgColsPlate$cols, xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf, alpha = 0.25)+
    geom_rect(data = scatterPeaksAllDf, aes(xmin = t1, xmax = t2), fill= "grey", ymin = -Inf,ymax = Inf, alpha = 0.4)+
    theme_classic()
  plts[[CsourceShort]][["scatterPeaksPlateGG"]] = scatterPlatePlt
  
  
  #summarize Peaks for groups by taking largest ranges for all peaks
  scatterPeaksDf = lapply(unique(scatterPeaksAllDf$groups), function(conc){
    submat = scatterPeaksAllDf[scatterPeaksAllDf$groups==conc, c("t1", "t2", "n1", "n2", "groups")]
    
    peakMat = submat[1,]
    
    for(sRow in 2:(nrow(submat))){
      subRow = submat[sRow,]
      brokenBool = F
      
      for(pRow in 1:nrow(peakMat)){ #do the ranges overlap?
        underBool = all(peakMat[pRow, 1:2] < subRow[,1])
        overBool = all(peakMat[pRow, 1:2] > subRow[,2])
        
        if(!(underBool | overBool)){
          peakMat[pRow, 1] = min(peakMat[pRow, 1], subRow[,1])
          peakMat[pRow, 3] = min(peakMat[pRow, 3], subRow[,3])
          peakMat[pRow, 2] = max(peakMat[pRow, 2], subRow[,2])
          peakMat[pRow, 4] = max(peakMat[pRow, 4], subRow[,4])
          
          brokenBool = T
          break
        }
      }
      if(!brokenBool){
        peakMat = rbind(peakMat, subRow)
      }
      return(peakMat)
    }
  })
  
  scatterPeaksDf = do.call(rbind, scatterPeaksDf)
  
  #mark peaks in following plots
  scatterPeaksGG = geom_rect(data = scatterPeaksDf, aes(xmin = t1, xmax = t2), fill= "grey", ymin = -Inf,ymax = Inf, alpha = 0.4)
  
  
  #only take the peaks beginning times for boxplots
  scatterPeaksPos = matrix(NA,length(scatterPeaksAllDf_appear),max(scatterPeaksAllDf_appear))
  row.names(scatterPeaksPos) = names(scatterPeaksAllDf_appear)
  
  for(k in names(scatterPeaksAllDf_appear)){
    columnPos = which(scatterPeaksAllDf$column == k)
    scatterPeaksPos[k,1:length(columnPos)] = scatterPeaksAllDf[columnPos, "t1"]
  }
  
  scatterPeaksPos = data.frame(scatterPeaksPos)
  colnames(scatterPeaksPos) = sub("X", "Peak",colnames(scatterPeaksPos))
  
  if(ncol(scatterPeaksPos)>=2){
    #also mark differences in peak times per well
    scatterPeaksPos_difs = sapply(2:ncol(scatterPeaksPos), function(colNum){
      scatterPeaksPos[,colNum] - scatterPeaksPos[,colNum-1]
    })
    
    #colnames(scatterPeaksPos_difs) = sapply(2:ncol(scatterPeaksPos), function(colNum){paste(colnames(scatterPeaksPos)[c(colNum-1,colNum)], collapse="_")})
    colnames(scatterPeaksPos_difs) = paste0("Delta*", colnames(scatterPeaksPos)[-1])
    
    scatterPeaksPos = cbind(scatterPeaksPos, scatterPeaksPos_difs)
  }
  
  
  #add the respectibe glucose level
  scatterPeaksPos$groups = sapply(names(scatterPeaksAllDf_appear), function(x){
    scatterPeaksAllDf[scatterPeaksAllDf$column == x, "groups"][1]
  })
  
  scatterPeaksPosDf = melt(scatterPeaksPos, id.vars = "groups")
  scatterPeaksPosDf = scatterPeaksPosDf[complete.cases(scatterPeaksPosDf),]
  
  
  bgColsBox = (concsUniq[-1] - concsUniq[-length(concsUniq)])/2
  bgColsBox = bgColsBox[c(1:length(bgColsBox), length(bgColsBox))]
  bgColsBoxDf = data.frame(xmin=concsUniq - bgColsBox, xmax=concsUniq+bgColsBox, cols= bgColsGG$data$col)
  
  scatterPeaksPosPlt = ggplot()+
    geom_point(data=scatterPeaksPosDf, aes(x=groups, y=value, group=groups))+
    geom_rect(data = bgColsBoxDf, aes(xmin =xmin, xmax=xmax, fill = cols), fill= rep(bgColsBoxDf$cols, length(unique(scatterPeaksPosDf$variable))),ymin = -Inf,ymax = Inf, alpha = 0.25)+
    geom_hline(yintercept = 0, linetype ="dashed")+
    labs(y="time of respective peak [h]", x="glucose concentration [C-mmol/l]")+
    scale_x_continuous(breaks=concsUniq, expand=c(0,0))+
    theme_classic()
  plts[[CsourceShort]][["scatterPeaksPosGG"]] = scatterPeaksPosPlt + if(length(unique(scatterPeaksPosDf$variable))>=2) {facet_grid(variable~., labeller = label_parsed)}
  
  
  
  
  #make ggplots for all measured variables and create smoothed counterparts
  for(j in list(c("scatter", " [AU]"),
                c("ribof", " [AU]"),
                c("O2", " [%]"),
                c("pH",""),
                c("NADH", " [AU]"))){
  
    varDat = getData(dat,j[1])
    
    pltDf = ggplotDf(varDat, concs, times)
    pltDf = pltDf[complete.cases(pltDf),]
    pltPlt = ggplot(pltDf)+
      geom_line(aes(x=times, y=value, group=as.factor(labels)))+
      facet_grid(groups~.)+
      scale_x_continuous(expand=c(0,0))+
      labs(x="time [h]", y=paste0(j, collapse=""))+
      bgColsGG+
      timeStartGG+
      labelFacetGG+
      scatterPeaksGG+
      theme_classic()
    
    if(j[1] %in% c("scatter", "ribof")) pltPlt = pltPlt + geom_hline(yintercept = 0, linetype ="dashed")
    
    plts[[CsourceShort]][[paste0(j[1],"GG")]] = pltPlt + if(CsourceShort == "ace" & j[1] == "pH"){coord_cartesian(ylim = c(5.7,7.5))}
    
    #smooth variable
    varSM = smootDat(varDat, dat$Time)
    
    smooths[[CsourceShort]][[j[1]]] = varSM$smthObj #save smoothing data
    
    dat <- addData(dat, paste0(j[1],"SM"), dat=varSM$smthDat)
  }
  
  # #find scatter peaks
  # scatterPeak = getData(dat, "scatter")
  # scatterPeak = apply(scatterPeak, 2, ma, 5)
  # #is the following value bigger (1) or smaller (-1)
  # scatterPeak = (rbind(scatterPeak[-1,], NA) - scatterPeak) 
  # scatterPeak = sign(scatterPeak)
  # 
  # #changes in this sign show peaks or valleys
  # scatterPeak = (rbind(NA, scatterPeak[-nrow(scatterPeak),]) + scatterPeak)
  # row.names(scatterPeak) = names(dat$Time)
  
  # scatterSM = smooths[[CsourceShort]][["scatter"]]
  # scatterDerif = sapply(scatterSM, function(smthdat){
  #   predict(smthdat, xarg = dat$Time ,nderiv = 1)
  #   })
  # scatterDerif = sapply(concsUniq, function(x){
  #   rowMeans(scatterDerif[,concs == x], na.rm = T)
  # })
  # 
  # 
  # scatterDerif2 = sapply(scatterSM, function(smthdat){
  #   predict(smthdat, xarg = dat$Time ,nderiv = 2)
  # })
  # scatterDerif2 = sapply(concsUniq, function(x){
  #   rowMeans(scatterDerif2[,concs == x], na.rm = T)
  # })
  
  
  
  scatter <- getData(dat, "scatter") # a.u., arbitrary unit
  gps <- 0.603 # mg/ml/scatter, biomass per scatter
  cpg <- 0.474 # g/g, gram carbon per gram biomass, Folsom&Carlson 2015
  mpc <- 1/12 # mol/g, mol carbon per gram carbon
  
  if(biomassEqivalent == "scatter"){
    #convert scatter to biomass
    biomass <- scatter * 1e3 * gps * cpg * mpc # mmol/L
    dat <- addData(dat, "biomass", dat=biomass)
  
    
  }else if(biomassEqivalent == "ribof"){
    #convert riboflavin to biomass
    ribof = getData(dat, "ribof")
    
    #find stable and linear ("saturated") areas near the end of measurement for scatter and riboflavin for each group of C-molarity
    SatBeginnings = data.frame(
      ribof = findSatGrouped(dat, "ribof", Amounts),
      scatter = findSatGrouped(dat, "scatter", Amounts)
    )
    
    SatBeginnings = apply(SatBeginnings, 1, max) #take the later saturation for each group
    
    #Each group of Concentrations evaluated at smaller saturation range (ribof or scatter)
    #SatBeginnings = rep(SatBeginnings,each = 3)[order(dat$layout$well)]
    
    #Alternatively: All groups evaluated in smallest saturation range
    SatBeginnings = rep(max(SatBeginnings), 3*length(SatBeginnings))[order(dat$layout$well)]
    
    
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
    
    plts[[CsourceShort]][["ribof_scatterGG"]]= ribofScatterPlt
    
    biomass <- ribof/ ribofScatterScale * 1e3 * gps * cpg * mpc   # C-mmol/L, scale riboflavin to scater and convert to biomass with scatter scaling
    dat <- addData(dat, "biomass", dat=biomass)
  }
  
  
  #raw plot
  viewGroups(dat, yids="biomass", show.ci95=TRUE, lwd.orig=0)
  plts[[CsourceShort]][["biomass"]]= recordPlot()
  
  #smooth biomass
  biomassSM = smootDat(biomass, dat$Time)
  
  smooths[[CsourceShort]][["biomass"]] = biomassSM$smthObj #save smoothing data
  
  dat <- addData(dat, "biomassSM", dat=biomassSM$smthDat)
  
  
  #ggplot
  biomassDf = ggplotDf(biomass, concs, times)
  biomassDf = biomassDf[complete.cases(biomassDf),]
  biomassPlt = ggplot(biomassDf)+
    geom_line(aes(x=times, y=value, group=as.factor(labels)))+
    facet_grid(groups~.)+
    geom_hline(yintercept = 0, linetype ="dashed")+
    scale_x_continuous(expand=c(0,0))+
    labs(x="time [h]", y="biomass [C-mmol/l]")+
    bgColsGG+
    timeStartGG+
    labelFacetGG+
    scatterPeaksGG+
    theme_classic()
  plts[[CsourceShort]][["biomassGG"]] = biomassPlt
  
  
  #get yield
  cmol <- dat$layout[[Amounts]]
  names(cmol) <- dat$layout$well
  
  cmol <- cmol[dat$layout$strain!="blank"] #exclude blanks
  
  SatBeginnings = findSatGrouped(dat, "biomass", Amounts) #find a good beginning for the saturation
  
  SatTime = dat$Time[max(SatBeginnings)]
  
  bm <- apply(getData(dat, "biomass", xrng=c(SatTime,max(dat$Time))),2,mean)[names(cmol)] #biomass between highest saturation time and and last time point
  
  yieldLmFull <- bestFrontLm(x = cmol, y = bm)
  yieldLm  = yieldLmFull$lm
  yieldCoef = coef(yieldLm)
  yield = yieldCoef[2]
  
  coefs[[CsourceShort]][["yield"]] = yield #save yield
  
  
  yieldDf = data.frame(cmol, bm)
  
  
  yieldPlt_label = c(sprintf("yield: %.3f", yield), sprintf(": %.3f", summary(yieldLm)$r.squared))
  yieldPlt_label = list(bquote(atop(.(yieldPlt_label[1]),~~~~R^2*.(yieldPlt_label[2]))))
  
  
  yieldPlt = ggplot()+
    geom_point(data = yieldDf, aes(x=cmol, y=bm))+
    geom_abline(intercept = yieldCoef[1] ,slope=yield, color="red", linetype="dashed")+
    labs(x=sprintf("%s concentration [C-mmol/l]",CSource), y="biomass [C-mmol/l]")+
    geom_label(label=yieldPlt_label, aes(x=max(cmol),y=0.5*max(bm)), hjust="right", color="red", parse=T)+
    scale_x_continuous(breaks = concsUniq)+
    theme_classic()
  
  plts[[CsourceShort]][["yieldGG"]]= yieldPlt
  
  
  #Oxygen
  O2 = getData(dat, "O2")
  O2 = O2/max(O2)*100
  
  O2Conc = apply(O2, 2, O2PercToConc) #mmol/l
  
  dat <- addData(dat, "O2conc", dat=O2Conc)
  
  
  # O2AddRate = apply(O2Conc, 2, O2concToRate, c_max = O2PercToConc(100), kLa = 260) #mmol/(l h)
  O2AddRate = (O2PercToConc(100) - O2Conc) * 260 #(c_max - c) * kLa (kLa = Oxygen transfer-coefficient, given by aeration setup) (mmol/(l h))
  
  O2Change = (rbind(O2Conc[-1,], NA) - O2Conc) #O2 difference of two consecutive datapoints (mmol/l)
  row.names(O2Change) = row.names(O2Conc)
  O2Change = O2Change/ c(dat$Time[-1] - dat$Time[-nTime], NA) #devide by time difference (mmol/(l h))
  
  O2MinusRate = O2Change - O2AddRate # difference in real and theoretical rate is consumption rate (mmol (l h))
  O2MinusRate = O2MinusRate *10^(-3) # mmol/(ml h) (*1 ml) -> mmol/h
  
  biomassG = biomass/cpg/ mpc * 10^(-6) #g/ml (*1ml) -> g
  
  O2Rate = O2MinusRate/biomassG # O2 consumpion rate of the bacteria (mmol/(h g))
  
  # max(O2Rate[!is.infinite(O2Rate)])
  # boxplot(O2Rate[!is.infinite(O2Rate)])
  # 
  # plot(y = O2Rate, x=matrix(dat$Time, nrow = nrow(O2RateNorm), ncol = ncol(O2RateNorm)))
  
  dat <- addData(dat, "O2Rate", dat=O2Rate)
  
  viewGroups(dat, yids="O2Rate", show.ci95=TRUE, lwd.orig=0)
  plts[[CsourceShort]][["O2Rate"]]= recordPlot()
  
  
  O2RateDf = ggplotDf(O2Rate, concs, times)
  O2RateDf = O2RateDf[complete.cases(O2RateDf),]
  
  O2RatePlt = ggplot(O2RateDf)+
    geom_line(aes(x=times, y=value, group=as.factor(labels)))+
    facet_grid(groups~.)+
    geom_hline(yintercept = 0, linetype ="dashed")+
    scale_x_continuous(expand=c(0,0))+
    labs(x="time [h]", y=expression("normalized"~O[2]~"change rate [mmol/(g h)]"))+
    bgColsGG+
    timeStartGG+
    labelFacetGG+
    scatterPeaksGG+
    theme_classic()
  
  O2RatePlt  = O2RatePlt + coord_cartesian(ylim = excludedScaling(O2RateDf, excltime = settleTime)) #Only after 0.5h usaful values, so restrict to these values scale
  plts[[CsourceShort]][["O2RateGG"]] = O2RatePlt
  
  
  
  #pH
  pH = getData(dat,"pH")
  
  HPlus = 10^(-pH) #mmol/ml
  HPlusRate = (rbind(HPlus[-1,], NA) - HPlus) #O2 difference of two consecutive datapoints (mmol/ml) (*1ml) -> mmol
  row.names(HPlusRate) = row.names(HPlus)
  HPlusRate = HPlusRate/ c(dat$Time[-1] - dat$Time[-nTime], NA) #devide by time difference (mmol/h)
  HPlusRate = HPlusRate/biomassG #H+ change per g biomass  mmol/(h g)
  
  dat <- addData(dat, "HPlusrate", dat=HPlusRate)
  
  HPlusRateDf = ggplotDf(HPlusRate, concs, times)
  HPlusRateDf = HPlusRateDf[complete.cases(HPlusRateDf),]
  
  HPlusRatePlt = ggplot(HPlusRateDf)+
    geom_line(aes(x=times, y=value, group=as.factor(labels)))+
    facet_grid(groups~.)+
    geom_hline(yintercept = 0, linetype ="dashed")+
    scale_x_continuous(expand=c(0,0))+
    labs(x="time [h]", y=bquote("normalized"~H^"+"*~"concentration change rate [mmol/(g h)])"))+
    bgColsGG+
    timeStartGG+
    labelFacetGG+
    scatterPeaksGG+
    theme_classic()
  
  HPlusRatePlt = HPlusRatePlt + coord_cartesian(ylim = c(0.001,-0.0005)) #only after 3h useful values
  plts[[CsourceShort]][["HPlusRateGG"]] = HPlusRatePlt
  
  #Hplus vs O2
  HPlsO2Df = cbind(HPlusRateDf, O2 = O2RateDf$value)
  
  HPlsO2Plt = ggplot(HPlsO2Df[HPlsO2Df$times >= settleTime,])+
    geom_point(aes(x=O2, y=value))+
    facet_grid(groups~.)+
    #geom_hline(yintercept = 0, linetype ="dashed")+
    geom_hline(yintercept = 0, linetype ="dashed")+
    labs(x=expression("normalized"~O[2]~"change rate [mmol/(g h)]"), y=bquote("normalized"~H^"+"*~"concentration change rate [mmol/(g h)])"))+
    bgColsGG+
    timeStartGG+
    labelFacetGG+
    theme_classic()
  
  plts[[CsourceShort]][["HPlusRate_O2GG"]] = HPlsO2Plt
  
  
  
  #==================================== dpseg analysis ====================================
  
  #growth rates in biomass data
  
  #biomassSM <- apply(biomass, 2, ma, 5)
  
  ## inspect individual well
  # id <- "B5"
  # if ( nutrient=="ace" ) id <- "E5"
  # plot(biomass[,id])
  # lines(biomassSM[,id],col=2, type="b", cex=.5)
  
  #dat <- addData(dat, ID="biomassSM", biomassSM)
  
  dpseg <- dpseg_plate(dat, "biomassSM",P=.0001)
  dpsegs[[CsourceShort]][["biomass"]] = dpseg
  
  ## add to plate express object to inspect
  mdat <- addModel(dpseg, dat, ID="biomassSM.dpseg")
  viewPlate(mdat, yids=c("biomassSM","biomassSM.dpseg"), log="y")
  viewGroups(mdat, yids=c("biomassSM.dpseg"),log="y")
  
  ## add growth rates
  mdat <- addModel(dpseg, mdat, ID="biomassSM.mu.dpseg", add.slopes=TRUE, col="#FF0000")
  
  ## inspect a single well:
  
  ### NOTE: requires knowledge of the data structures
  ### use RStudio data browser to inspect
  # id <- "A1"
  # if ( nutrient=="ace" ) id <- "E5"
  # 
  # par(mfcol=c(1,1), mai=c(.5,.5,.1,.5),mgp=c(1.3,.4,0),tcl=-.25)
  # plot(mdat$biomassSM$data[,id],log="y",cex=.5)
  # lines(mdat$biomassSM.dpseg$data[,id], lwd=2)
  # par(new=TRUE)
  # plot(mdat$biomass.mu.dpseg$data[,id], axes=FALSE, type="l",ylab=NA)
  # abline(v=dpseg[[id]]$segments[,"start"])
  # abline(v=dpseg[[id]]$segments[,"end"], col=2)
  # axis(4)
  # mtext("growth rate, h-1", 4, par("mgp")[1])
  
  ## inspect to get growth rates
  
  viewPlate(mdat, yids=c("biomassSM","biomassSM.mu.dpseg"))
  plts[[CsourceShort]][["growthratesBiomassPlate"]]= recordPlot()
  
  viewGroups(mdat, yids=c("biomassSM.mu.dpseg"), show.ci95=FALSE, lwd.orig=0,
             #ylim=c(-.01,.5),
             xlab="time, h", ylab=expression("biomass growth rate"~mu*","~h^-1),
             embed=TRUE,no.par=TRUE, g1.legend=FALSE, g2.legend=FALSE)
  abline(h=0, lwd=2)
  plts[[CsourceShort]][["growthratesBiomassAll"]]= recordPlot()
  
  
  #ggplot
  growthratesBiomassDf = ggplotDf(getData(mdat, "biomassSM.mu.dpseg"), concs, times)
  growthratesBiomassDf = growthratesBiomassDf[complete.cases(growthratesBiomassDf),]
  
  growthratesBiomassPlt = ggplot(growthratesBiomassDf)+
    geom_line(aes(x=times, y=value, group=as.factor(labels)))+
    facet_grid(groups~.)+
    geom_hline(yintercept = 0, linetype ="dashed")+
    scale_x_continuous(expand=c(0,0))+
    labs(x="time [h]", y=bquote("biomass growth rate "~mu~" ["~h^{-1}~"]"))+
    bgColsGG+
    timeStartGG+
    labelFacetGG+
    scatterPeaksGG+
    theme_classic()
  
  plts[[CsourceShort]][["growthratesBiomassGG"]] = growthratesBiomassPlt + coord_cartesian(ylim = excludedScaling(growthratesBiomassDf, excltime = settleTime))
  
  
  
  #growth rates in raw riboflavin data
  
  #ribofSM <- apply(ribof, 2, ma, 5)
  
  ## inspect individual well
  # id <- "B5"
  # if ( nutrient=="ace" ) id <- "E5"
  # plot(ribof[,id])
  # lines(ribofSM[,id],col=2, type="b", cex=.5)
  
  #dat <- addData(dat, ID="ribofSM", ribofSM)
  
  
  dpseg <- dpseg_plate(dat, "ribofSM",P=.0001)
  dpsegs[[CsourceShort]][["ribof"]] = dpseg
  
  ## add to plate express object to inspect
  mdat <- addModel(dpseg, mdat, ID="ribofSM.dpseg")
  viewPlate(mdat, yids=c("ribofSM","ribofSM.dpseg"), log="y")
  viewGroups(mdat, yids=c("ribofSM.dpseg"),log="y")
  
  ## add growth rates
  mdat <- addModel(dpseg, mdat, ID="ribofSM.mu.dpseg", add.slopes=TRUE, col="#FF0000")
  
  ## inspect a single well:
  
  ### NOTE: requires knowledge of the data structures
  ### use RStudio data browser to inspect
  # id <- "B5"
  # if ( nutrient=="ace" ) id <- "E5"
  # 
  # par(mfcol=c(1,1), mai=c(.5,.5,.1,.5),mgp=c(1.3,.4,0),tcl=-.25)
  # plot(mdat$ribofSM$data[,id],log="y",cex=.5)
  # lines(mdat$ribofSM.dpseg$data[,id], lwd=2)
  # par(new=TRUE)
  # plot(mdat$ribofSM.mu.dpseg$data[,id], axes=FALSE, type="l",ylab=NA)
  # abline(v=dpseg[[id]]$segments[,"start"])
  # abline(v=dpseg[[id]]$segments[,"end"], col=2)
  # axis(4)
  # mtext("growth rate, h-1", 4, par("mgp")[1])
  
  ## inspect to get growth rates
  viewPlate(mdat, yids=c("ribofSM","ribofSM.mu.dpseg"))
  plts[[CsourceShort]][["growthratesRibofPlate"]]= recordPlot()
  
  
  viewGroups(mdat, yids=c("ribofSM.mu.dpseg"), show.ci95=FALSE, lwd.orig=0,
             #ylim=c(-.01,.1),
             xlab="time, h", ylab=expression("riboflavine growth rate"~mu*","~h^-1),
             embed=TRUE,no.par=TRUE, g1.legend=FALSE, g2.legend=FALSE)
  abline(h=0, lwd=2)
  plts[[CsourceShort]][["growthratesRibofAll"]]= recordPlot()
  
  
  #ggplot
  growthratesRibofDf = ggplotDf(getData(mdat, "ribofSM.mu.dpseg"), concs, times)
  growthratesRibofDf = growthratesRibofDf[complete.cases(growthratesRibofDf),]
  
  growthratesRibofPlt = ggplot(growthratesRibofDf)+
    geom_line(aes(x=times, y=value, group=as.factor(labels)))+
    facet_grid(groups~.)+
    geom_hline(yintercept = 0, linetype ="dashed")+
    scale_x_continuous(expand=c(0,0))+
    labs(x="time [h]", y=bquote("riboflavine growth rate "~mu~" ["~h^{-1}~"]"))+
    bgColsGG+
    timeStartGG+
    labelFacetGG+
    scatterPeaksGG+
    theme_classic()
  
  plts[[CsourceShort]][["growthratesRibofGG"]] = growthratesRibofPlt + coord_cartesian(ylim = excludedScaling(growthratesRibofDf, excltime = settleTime))
  
  
  #growth rates in raw scatter data
  
  #scatterSM <- apply(scatter, 2, ma, 5)
  
  ## inspect individual well
  # id <- "B5"
  # if ( nutrient=="ace" ) id <- "E5"
  # plot(scatter[,id])
  # lines(scatterSM[,id],col=2, type="b", cex=.5)
  
  #dat <- addData(dat, ID="scatterSM", scatterSM)
  
  dpseg <- dpseg_plate(dat, "scatterSM",P=.0001)
  dpsegs[[CsourceShort]][["scatterSM"]] = dpseg
  
  ## add to plate express object to inspect
  mdat <- addModel(dpseg, mdat, ID="scatterSM.dpseg")
  viewPlate(mdat, yids=c("scatterSM","scatterSM.dpseg"), log="y")
  viewGroups(mdat, yids=c("scatterSM.dpseg"),log="y")
  
  ## add growth rates
  mdat <- addModel(dpseg, mdat, ID="scatterSM.mu.dpseg", add.slopes=TRUE, col="#FF0000")
  
  ## inspect a single well:
  
  ### NOTE: requires knowledge of the data structures
  ### use RStudio data browser to inspect
  # id <- "B5"
  # if ( nutrient=="ace" ) id <- "E5"
  # 
  # par(mfcol=c(1,1), mai=c(.5,.5,.1,.5),mgp=c(1.3,.4,0),tcl=-.25)
  # plot(mdat$scatterSM$data[,id],log="y",cex=.5)
  # lines(mdat$scatterSM.dpseg$data[,id], lwd=2)
  # par(new=TRUE)
  # plot(mdat$scatterSM.mu.dpseg$data[,id], axes=FALSE, type="l",ylab=NA)
  # abline(v=dpseg[[id]]$segments[,"start"])
  # abline(v=dpseg[[id]]$segments[,"end"], col=2)
  # axis(4)
  # mtext("growth rate, h-1", 4, par("mgp")[1])
  
  ## inspect to get growth rates
  viewPlate(mdat, yids=c("scatterSM","scatterSM.mu.dpseg"))
  plts[[CsourceShort]][["growthratesScatterPlate"]]= recordPlot()
  
  
  viewGroups(mdat, yids=c("scatterSM.mu.dpseg"), show.ci95=FALSE, lwd.orig=0,
             #ylim=c(-.01,.1),
             xlab="time, h", ylab=expression("scatter growth rate"~mu*","~h^-1),
             embed=TRUE,no.par=TRUE, g1.legend=FALSE, g2.legend=FALSE)
  abline(h=0, lwd=2)
  plts[[CsourceShort]][["growthratesScatterAll"]]= recordPlot()
  
  
  #ggplot
  growthratesScatterDf = ggplotDf(getData(mdat, "scatterSM.mu.dpseg"), concs, times)
  growthratesScatterDf = growthratesScatterDf[complete.cases(growthratesScatterDf),]
  
  growthratesScatterPlt = ggplot(growthratesScatterDf)+
    geom_line(aes(x=times, y=value, group=as.factor(labels)))+
    facet_grid(groups~.)+
    geom_hline(yintercept = 0, linetype ="dashed")+
    scale_x_continuous(expand=c(0,0))+
    labs(x="time [h]", y=bquote("scatter growth rate "~mu~" ["~h^{-1}~"]"))+
    bgColsGG+
    timeStartGG+
    labelFacetGG+
    scatterPeaksGG+
    theme_classic()
  
  plts[[CsourceShort]][["growthratesScatterGG"]] = growthratesScatterPlt + coord_cartesian(ylim = excludedScaling(growthratesScatterDf, excltime = settleTime))
  
  
  ## TODO - platexpress: fix orig colors!
  # colors.fixed <- FALSE
  # if ( colors.fixed ) {
  #   par(mfcol=c(2,1))
  #   viewGroups(mdat, yids=c("scatterSM"), log="y",embed=TRUE,
  #              show.ci95=FALSE, lwd.orig=1)
  #   viewGroups(mdat, yids=c("mu.dpseg"), embed=TRUE, ylim=c(-.01,.22),
  #              g2.legend=FALSE, show.ci95=FALSE, lwd.orig=1)
  #   abline(h=0, lwd=2)
  
  growtratesO2Df = cbind(growthratesBiomassDf, O2 = O2RateDf[O2RateDf$times %in% unique(growthratesBiomassDf$times), "value"])
  
  growtratesO2Plt = ggplot(growtratesO2Df[growtratesO2Df$times >= settleTime,])+
    geom_point(aes(x=O2, y=value))+
    facet_grid(groups~.)+
    #geom_hline(yintercept = 0, linetype ="dashed")+
    geom_hline(yintercept = 0, linetype ="dashed")+
    labs(x=expression("normalized"~O[2]~"change rate [mmol/(g h)]"), y=bquote("biomass growth rate "~mu~" ["~h^{-1}~"]"))+
    bgColsGG+
    labelFacetGG+
    theme_classic()
  
  plts[[CsourceShort]][["growthrates_O2GG"]] = growtratesO2Plt
  
  
  #model maximal growhrates against
  #dpsegs[[Amounts]]$biomass
  
  #save data
  dats[[CsourceShort]] <- mdat
  plts[[CsourceShort]]$mods = list("bgColsGG" = bgColsGG,
                                   "timeStartGG" = timeStartGG,
                                   "labelFacetGG" = labelFacetGG,
                                   "scatterPeaksGG" = scatterPeaksGG)
}


savePlts(plts)

save.image(file = paste0(INPATH, "/biolector_evalDAT.RData"))


#==================================== modelling ====================================

# dat = dats$glc
# 
# dfg <- data2growthrates(dat, "biomass", plate=dat$layout)
# 
# #initial values
# yinit <- c(s=270, y=1)
# poinit <- c(phi=0, sin=0, mumax=.25, K=5, Y=0.35, d=0)