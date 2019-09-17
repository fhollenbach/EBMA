#'@include summary_plot.R
NULL

##
#' Print and Show methods for forecast data
#'
#' Functions to print and show the contents of a data object of the class 'ForecastData' or 'SummaryForecastData'.
#'
#' @param object An object of the class 'ForecastData' or 'SummaryForecastData'.
#' @param x An object of the class 'ForecastData' or 'SummaryForecastData'.
#' @param digits An integer specifying the number of significant digits to print. The default is 3.
#' @param ... Not implemented
#'
#'
#' @return NULL
#' 
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2015). Calibrating ensemble forecasting models with sparse data in the social sciences.   \emph{International Journal of Forecasting}. \bold{31}: 930-942.
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2012). Improving Predictions Using Ensemble Bayesian Model Averaging. \emph{Political Analysis}. \bold{20}: 271-291.
#'
#' @author  Michael D. Ward <\email{michael.d.ward@@duke.edu}> and Jacob M. Montgomery <\email{jacob.montgomery@@wustl.edu}> and Florian M. Hollenbach <\email{florian.hollenbach@@tamu.edu}>
#'
#' @examples \dontrun{ data(calibrationSample)
#'
#' data(testSample) 
#' 
#' this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],
#' .outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],
#' .outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
#' 
#' this.ensemble <- calibrateEnsemble(this.ForecastData, model="logit", tol=0.001,exp=3)
#' 
#' summary.object <- summary(this.ensemble, period="calibration") 
#' print(summary.object)
#' show(summary.object)
#'}
#' @rdname PrintShow
#' @export
setMethod(
  f="print",
  signature="SummaryForecastData",
  definition=function(x, digits=3, ...){
    print(x@summaryData, na.print="", digits=digits)
  }
)

#' @rdname PrintShow
#' @export
setMethod(
  f="show",
  signature="SummaryForecastData",
  definition=function(object){
    print(object@summaryData, na.print="", digits=3)
  }
)

#' @rdname PrintShow
#' @export
setMethod(
  f="print",
  signature="ForecastData",
  definition=function(x, digits=3, ...){
    cat("* Prediction Calibration = \n"); 
    if(length(x@predCalibration)>0)
    {print(x@predCalibration, na.print="", digits=digits);}
    else{print("Nothing Here")}
    cat("* Prediction Test = \n"); 
    if(length(x@predTest)>0)
    {print(x@predTest, na.print="", digits=digits);}
    else{print("Nothing Here")}
    cat("* Outcome Calibration = \n");
    if(length(x@outcomeCalibration)>0)
    {print(x@outcomeCalibration, na.print="", digits=digits);}
    else{print("Nothing Here")}
    cat("* Outcome Test = \n");
    if(length(x@outcomeTest)>0)
    {print(x@outcomeTest, na.print="", digits=digits);}
    else{print("Nothing Here")}
    cat("* Model Names = \n ");print(x@modelNames, na.print="");
  }
)

#' @rdname PrintShow
#' @export
setMethod(
  f="show",
  signature="ForecastData",
  definition=function(object){
    if (length(object@predCalibration)==0) {
      cat("* Prediction Calibration = \n"); 
      if(length(object@predCalibration)>0)
      {print(object@predCalibration, na.print="", digits=1);}
      else{print("Nothing Here")}
      cat("* Prediction Test = \n"); 
      if(length(object@predTest)>0)
      {print(object@predTest, na.print="", digits=1);}
      else{print("Nothing Here")}
      cat("* Outcome Calibration = \n");
      if(length(object@outcomeCalibration)>0)
      {print(object@outcomeCalibration, na.print="", digits=1);}
      else{print("Nothing Here")}
      cat("* Outcome Test = \n");
      if(length(object@outcomeTest)>0)
      {print(object@outcomeTest, na.print="", digits=1);}
      else{print("Nothing Here")}
      cat("* Model Names = \n ");print(object@modelNames, na.print="");
    }
    else{
      nrowCal=min(10,nrow(object@predCalibration))
      nrowTest=min(10,nrow(object@predTest))
      cat("* Prediction Calibration = \n"); 
      if(length(object@predCalibration)>0)
      {print(object@predCalibration[1:nrowCal,1:ncol(object@predCalibration),1], na.print="", digits=2);}
      else{print("Nothing Here")}
      cat("* Prediction Test = \n"); 
      if(length(object@predTest)>0)
      {print(object@predTest[1:nrowTest,1:ncol(object@predTest),1], na.print="", digits=2);}
      else{print("Nothing Here")}
      cat("* Outcome Calibration = \n");
      if(length(object@outcomeCalibration)>0)
      {print(print(object@outcomeCalibration[1:nrowCal]),na.print="", digits=2);}
      else{print("Nothing Here")}
      cat("* Outcome Test = \n");
      if(length(object@outcomeTest)>0)
      {print(object@outcomeTest[1:nrowTest], na.print="", digits=2);}
      else{print("Nothing Here")}
      cat("* Model Names = \n ");print(object@modelNames,na.print="");
    }
  }
)
