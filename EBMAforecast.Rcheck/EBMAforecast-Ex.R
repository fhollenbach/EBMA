pkgname <- "EBMAforecast"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "EBMAforecast-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('EBMAforecast')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("EBMAforecast")
### * EBMAforecast

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: EBMAforecast
### Title: EBMAforecast
### Aliases: EBMAforecast EBMAforecast-package

### ** Examples

## Not run: 
##D demo(EBMAforecast)
##D demo(presForecast)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("EBMAforecast", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ForecastData")
### * ForecastData

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ForecastData-class
### Title: An ensemble forecasting data object
### Aliases: ForecastData-class

### ** Examples

## Not run: 
##D  
##D data(calibrationSample)
##D data(testSample)
##D 
##D this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],
##D .outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],
##D .outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
##D 
##D ### to acces individual slots in the ForecastData object
##D getPredCalibration(this.ForecastData)
##D getOutcomeCalibration(this.ForecastData)
##D getPredTest(this.ForecastData)
##D getOutcomeTest(this.ForecastData)
##D getModelNames(this.ForecastData)
##D 
##D ### to assign individual slots, use set functions
##D 
##D setPredCalibration(this.ForecastData)<-calibrationSample[,c("LMER", "SAE", "GLM")]
##D setOutcomeCalibration(this.ForecastData)<-calibrationSample[,"Insurgency"]
##D setPredTest(this.ForecastData)<-testSample[,c("LMER", "SAE", "GLM")]
##D setOutcomeTest(this.ForecastData)<-testSample[,"Insurgency"]
##D setModelNames(this.ForecastData)<-c("LMER", "SAE", "GLM")
## End(Not run)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ForecastData", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PrintShow")
### * PrintShow

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print,SummaryForecastData-method
### Title: Print and Show methods for forecast data
### Aliases: print,SummaryForecastData-method print,
###   SummaryForecastData-method show, show,SummaryForecastData-method
###   print,ForecastData-method show,ForecastData-method

### ** Examples

## Not run: 
##D  data(calibrationSample)
##D 
##D data(testSample) 
##D 
##D this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],
##D .outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],
##D .outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
##D 
##D this.ensemble <- calibrateEnsemble(this.ForecastData, model="logit", tol=0.001,exp=3)
##D 
##D summary.object <- summary(this.ensemble, period="calibration") 
##D print(summary.object)
##D show(summary.object)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PrintShow", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("SummarizePlot")
### * SummarizePlot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary,FDatFitLogit-method
### Title: Summarize and Plot Ensemble models
### Aliases: summary,FDatFitLogit-method summary,F DatFitLogit-method
###   summary, FDatFitNormal-method plot, FDatFitLogit-method
###   summary,FDatFitNormal-method plot,FDatFitLogit-method
###   plot,FDatFitNormal-method

### ** Examples

## Not run: 
##D  data(calibrationSample)
##D data(testSample) 
##D 
##D this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],
##D .outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],
##D .outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
##D 
##D this.ensemble <- calibrateEnsemble(this.ForecastData, model="logit", tol=0.001,exp=3)
##D 
##D summary(this.ensemble, period="calibration") 
##D 
##D summary(this.ensemble, period="test",showCoefs=FALSE)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("SummarizePlot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calibrateEnsemble")
### * calibrateEnsemble

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calibrateEnsemble
### Title: Calibrate an ensemble Bayesian Model Averaging model
### Aliases: calibrateEnsemble fitEnsemble, ForecastData-method
###   ForecastDataLogit-method ForecastDataNormal-method,
###   FDatFitLogit-class, ForecastDataLogit-class,
###   ForecastDataNormal-class, FDatFitNormal-class calibrateEnsemble,

### ** Examples

## Not run: 
##D data(calibrationSample)
##D 
##D data(testSample)
##D 
##D this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],
##D .outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],
##D .outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
##D 
##D this.ensemble <- calibrateEnsemble(this.ForecastData, model="logit", tol=0.001, exp=3)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calibrateEnsemble", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compareModels")
### * compareModels

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compareModels
### Title: Function for comparing multiple models based on predictive
###   performance
### Aliases: compareModels compareModels,

### ** Examples

## Not run: 
##D data(calibrationSample)
##D 
##D data(testSample) 
##D 
##D this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],
##D .outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],
##D .outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
##D 
##D this.ensemble <- calibrateEnsemble(this.ForecastData, model="logit", tol=0.001, exp=3)
##D 
##D compareModels(this.ensemble,"calibration")
##D 
##D compareModels(this.ensemble,"test") 
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compareModels", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("makeForecastData")
### * makeForecastData

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: makeForecastData
### Title: Build a ensemble forecasting data object
### Aliases: makeForecastData

### ** Examples


## Not run: 
##D data(calibrationSample)
##D data(testSample)
##D this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],
##D .outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],
##D .outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
##D 
##D ### to acces individual slots in the ForecastData object
##D getPredCalibration(this.ForecastData)
##D getOutcomeCalibration(this.ForecastData)
##D getPredTest(this.ForecastData)
##D getOutcomeTest(this.ForecastData)
##D getModelNames(this.ForecastData)
##D 
##D ### to assign individual slots, use set functions
##D 
##D setPredCalibration(this.ForecastData)<-calibrationSample[,c("LMER", "SAE", "GLM")]
##D setOutcomeCalibration(this.ForecastData)<-calibrationSample[,"Insurgency"]
##D setPredTest(this.ForecastData)<-testSample[,c("LMER", "SAE", "GLM")]
##D setOutcomeTest(this.ForecastData)<-testSample[,"Insurgency"]
##D setModelNames(this.ForecastData)<-c("LMER", "SAE", "GLM")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("makeForecastData", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("setFunctions")
### * setFunctions

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: setPredCalibration<-
### Title: "Set" functions
### Aliases: setPredCalibration<- setFunctions setPredCalibration
###   setOutcomeCalibration setPredTest setOutcomeTest setModelNames
###   setPredTest<- setOutcomeCalibration<- setOutcomeTest<-
###   setModelNames<-

### ** Examples


## Not run: 
##D data(calibrationSample)
##D data(testSample)
##D setPredCalibration(this.ForecastData)<-calibrationSample[,c("LMER", "SAE", "GLM")]
##D setOutcomeCalibration(this.ForecastData)<-calibrationSample[,"Insurgency"]
##D setPredTest(this.ForecastData)<-testSample[,c("LMER", "SAE", "GLM")]
##D setOutcomeTest(this.ForecastData)<-testSample[,"Insurgency"]
##D setModelNames(this.ForecastData)<-c("LMER", "SAE", "GLM")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("setFunctions", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
