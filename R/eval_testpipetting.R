library(reshape2)
library(pracma)
library(ggplot2)

#functions
integrate_points = function(y, increment = 1){
  left = y[1:length(y)-1]
  right = y[2:length(y)]
  
  rect = sum(left)
  
  tri = 0.5 * sum(right - left)
  
  res = (rect + tri) * increment
  return(res)
}

#read matrix
pipMat = read.csv2("C://Users/tpeng/OneDrive/OD__Universität/BA__Axmann/testpipetting1_data.csv", header=T, skip=10, stringsAsFactors = F)

#restructure matrix
pipMat_wl = as.numeric(pipMat[1,-1:-2])

pipMat = pipMat[-1,]

pipMat_wells = pipMat$Well
pipMat_type = pipMat$Content

pipMat = pipMat[,-1:-2]

colnames(pipMat) = pipMat_wl
row.names(pipMat) = pipMat_wells


#define layout
pipMat_layout = rep(c("R_base", "R_base","R_grad","R_grad","R_grad", NA,
                      "H_base", "H_base","H_grad","H_grad","H_grad", NA),8)

pipMat_meth = c()
for(i in 1:length(pipMat_layout)){
  x = pipMat_layout[i]
  splt = strsplit(x,"_")[[1]]
  pipMat_meth = rbind(pipMat_meth,splt)
}

pipMat_conc_base = rep(50,8)
pipMat_conc_grad = linspace(0,100,8)

pipMat_conc = cbind(pipMat_conc_base,
                    pipMat_conc_base,
                    pipMat_conc_grad,
                    pipMat_conc_grad,
                    pipMat_conc_grad,
                    NA)

pipMat_conc = cbind(pipMat_conc,pipMat_conc)

pipMat_conc = as.vector(t(pipMat_conc))


#integrate matrix
pipMat_int = apply(pipMat,1, integrate_points)


#result matrix
pipMat_res = data.frame(integral = pipMat_int, type = pipMat_layout, conc = pipMat_conc)
pipMat_res = cbind(pipMat_res, pipMat_meth)
colnames(pipMat_res)[c(4,5)] = c("experimentor", "method")

pipMat_res = pipMat_res[complete.cases(pipMat_res),]

#plot
pipMat_baseplt =
ggplot(pipMat_res[pipMat_res$method=="base",], aes(x=experimentor,y=integral))+
  geom_boxplot()+
  geom_point()+
  ggtitle("Baseline")

pipMat_gradplt =
ggplot(pipMat_res[pipMat_res$method=="grad",], aes(x=conc, y=integral, group=conc))+
  geom_boxplot()+
  geom_point()+
  facet_grid(experimentor~.)+
  ggtitle("Gradient")

pipMat_baseplt
pipMat_gradplt
