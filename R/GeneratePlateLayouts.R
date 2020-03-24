strain = "Ec_DH5alphaZ1"

ConcRange = c(0,0.36)

# ===========================================================================================================================

plateSpecsGrad = round(linspace(ConcRange[1],ConcRange[2],7), 3)

plateSpecsGradGlc = c(rep(plateSpecsGrad, 3), rep(0, 21))
plateSpecsGradAce = c(rep(0, 21), rep(plateSpecsGrad, 3))

plateLayout = sprintf("%s;Glc:%.3f;Ace:%.3f;aTc:0", strain, plateSpecsGradGlc, plateSpecsGradAce)
plateLayout = matrix(plateLayout, nrow = 6, ncol = 7, byrow = T)
plateLayout = cbind(plateLayout, sprintf("blank;Glc:%.3f;Ace:%.3f;aTc:0", c(0,0.18,0.36,0,0,0), c(0,0,0,0,0.18,0.36)))

dimnames(plateLayout) = list(LETTERS[1:6], 1:8)

write.csv2(plateLayout, file = "C://Users/tpeng/OneDrive/OD__Universität/BA__Axmann/data/25_Ecoli_2020_REDUCTION-1_layout.csv")
