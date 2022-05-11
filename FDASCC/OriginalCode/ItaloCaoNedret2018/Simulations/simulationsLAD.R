#Peak ----
## ## ## ## ## ## ## #### ## ## ##
# setwd("Final Version/")

set.seed(90210)

cont = seq(0.05, 0.5, by = 0.05)
sampleSize = c(30, 50, 100, 200, 300)
# sampleSize = c(30, 50, 100)
# sampleSize = c(100)
resultsPeak <- matrix(0, ncol = length(cont),nrow = length(sampleSize))
colnames(resultsPeak) <- cont
rownames(resultsPeak) <- sampleSize

source("RobustMethod/LAD/peakSimulationLAD_smoothFirst.R")

for (ii in 1:length(sampleSize)) {
  for(jj in 1:length(cont)) {
    resultsPeak[ii, jj] <- peakSimLADSmoothFirst(contamination = cont[jj], n = sampleSize[ii], strength = 5) 
    resultsPeak[ii,jj]
    print(jj)
  }
  print(ii)
}

write.csv(resultsPeak, "resultsPeak.csv")
# resultsPeak = read.csv("results/LAD/resultsPeak", header = TRUE, row.names = 1)


#Bump ----
## ## ## ## ## ## ## #### ## ## ##
set.seed(90210)
resultsBump <- matrix(0, ncol = length(cont),nrow = length(sampleSize))
colnames(resultsBump) <- cont
rownames(resultsBump) <- sampleSize

source("RobustMethod/LAD/bumpSimulationLAD_smoothFirst.R")

for (ii in 1:length(sampleSize)) {
  for(jj in 1:length(cont)) {
    resultsBump[ii, jj] <- bumpSimLADSmoothFirst(contamination = cont[jj], n = sampleSize[ii], strength = 5) 
    resultsBump[ii,jj]
    print(jj)
  }
  print(ii)
}

write.csv(resultsBump, "resultsBump.csv")
# resultsBump = read.csv("results/LAD/resultsBump", header = TRUE, row.names = 1)


#Step ----
## ## ## ## ## ## ## #### ## ## ##
set.seed(90210)
resultsStep <- matrix(0, ncol = length(cont),nrow = length(sampleSize))
colnames(resultsStep) <- cont
rownames(resultsStep) <- sampleSize

source("RobustMethod/LAD/stepSimulationLAD_smoothFirst.R")

for (ii in 1:length(sampleSize)) {
  for(jj in 1:length(cont)) {
    resultsStep[ii, jj] <- stepSimLADSmoothFirst(contamination = cont[jj], n = sampleSize[ii], strength = 5) 
    resultsStep[ii,jj]
    print(jj)
  }
  print(ii)
}

write.csv(resultsStep, "resultsStep.csv")
# resultsStep = read.csv("results/LAD/resultsStep", header = TRUE, row.names = 1)


#Symmetric ----
## ## ## ## ## ## ## #### ## ## ##
set.seed(90210)
resultsSym <- matrix(0, ncol = length(cont),nrow = length(sampleSize))
colnames(resultsSym) <- cont
rownames(resultsSym) <- sampleSize

source("RobustMethod/LAD/symmetricSimulationLAD_smoothFirst.R")

for (ii in 1:length(sampleSize)) {
  for(jj in 1:length(cont)) {
    resultsSym[ii, jj] <- symmetricSimLADSmoothFirst(contamination = cont[jj], n = sampleSize[ii], strength = 5) 
    resultsSym[ii,jj]
    print(jj)
  }
  print(ii)
}

write.csv(resultsSym, "resultsSym.csv")
# resultsSym = read.csv("results/LAD/resultsSym", header = TRUE, row.names = 1)


#Translation ----
## ## ## ## ## ## ## #### ## ## ##
set.seed(90210)
resultsTranslation <- matrix(0, ncol = length(cont),nrow = length(sampleSize))
colnames(resultsTranslation) <- cont
rownames(resultsTranslation) <- sampleSize

source("RobustMethod/LAD/translationSimulationLAD_smoothFirst.R")

for (ii in 1:length(sampleSize)) {
  for(jj in 1:length(cont)) {
    resultsTranslation[ii, jj] <- translationSimLADSmoothFirst(contamination = cont[jj], n = sampleSize[ii], strength = 5) 
    resultsTranslation[ii,jj]
    print(jj)
  }
  print(ii)
}

write.csv(resultsTranslation, "resultsTranslation.csv")
# resultsTranslation = read.csv("results/LAD/resultsTranslation", header = TRUE, row.names = 1)
