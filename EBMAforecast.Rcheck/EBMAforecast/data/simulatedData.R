# # create test data for: logit
# data(calibrationSample)
# data(testSample)
# data("presidentialForecast")
# 
# # should be .5 for the first two models
# testLogit <- rbinom(length(calibrationSample[,1]), 1, prob=.5*(calibrationSample[,1]+calibrationSample[,2]))
# save(testLogit, file="EBMAforecast/data/simulatedLogitData.rda")
# 
# this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],
#                                       .outcomeCalibration=testObserved,
#                                       .predTest=testSample[,c("LMER", "SAE", "GLM")],
#                                       .outcomeTest=testSample[,"Insurgency"],
#                                       .modelNames=c("LMER", "SAE", "GLM"))
# 
# this.ensemble <- calibrateEnsemble(this.ForecastData, model="logit", tol=0.0001, maxIter=25000, exp=3)
# this.ensemble@modelWeights
# 
# # create test data for: normal
# n <- 20
# test.forecasts <- data.frame(matrix(rep(t(presidentialForecast[,c(1:6)]),n),
#                          ncol=ncol(presidentialForecast[,c(1:6)]), byrow=TRUE))
# 
# # draw observations from different models
# testObserved1 <- rnorm((n*15)/2, test.forecasts[,1], 1)
# testObserved2 <- rnorm((n*15)/2, test.forecasts[,2], 1)
# testObserved <- c(testObserved1, testObserved2)
# testNorm <- cbind(test.forecasts, testObserved)
# 
# save(testNorm, file="EBMAforecast/data/simulatedNormData.rda")
# 
# test.ForecastData<-makeForecastData(.predCalibration=test.forecasts[c(1:(n*15)-1),],
#                                     .outcomeCalibration=testObserved[c(1:(n*15)-1)],
#                                     .predTest=test.forecasts[(n*15),],
#                                     .outcomeTest=testObserved[(n*15)],
#                                     .modelNames=c("Campbell", "Lewis-Beck","EWT2C2","Fair","Hibbs","Abramowitz"))
# thisEnsemble<-calibrateEnsemble(test.ForecastData, model="normal", useModelParams=FALSE, tol=0.000000001)
# 
# thisEnsemble@modelWeights
# testthat