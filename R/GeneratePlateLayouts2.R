strain = "Ec_DH5alphaZ1"

# ===========================================================================================================================

plateSpecsGradGlc = rep(0.36, 48)
plateSpecsGradAce = 0

plateLayout = sprintf("%s;Glc:%.3f;Ace:%.3f;aTc:0", strain, plateSpecsGradGlc, plateSpecsGradAce)
plateLayout = matrix(plateLayout, nrow = 6, ncol = 8, byrow = T)
plateLayout[4:6,8] = sub(strain, "blank", plateLayout[4:6,8])

dimnames(plateLayout) = list(LETTERS[1:6], 1:8)

write.csv2(plateLayout, file = "C://Users/tpeng/OneDrive/OD__Universität/BA__Axmann/data/27_Ecoli_2020_REDUCTION-1_layout.csv")
