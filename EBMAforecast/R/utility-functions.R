#' "Set" functions
#'
#' To assign individual slots, use set functions
#' @method set setPredCalibration
#' @method set setOutcomeCalibration
#' @method set setPredTest
#' @method set setOutcomeTest
#' @method set setModelNames
#' 
#' @param object The object to which values are assigned.
#' @param value Values to be assigned.
#' 
#' @examples
#' 
#' \dontrun{
#' data(calibrationSample)
#' data(testSample)
#' setPredCalibration(this.ForecastData)<-calibrationSample[,c("LMER", "SAE", "GLM")]
#' setOutcomeCalibration(this.ForecastData)<-calibrationSample[,"Insurgency"]
#' setPredTest(this.ForecastData)<-testSample[,c("LMER", "SAE", "GLM")]
#' setOutcomeTest(this.ForecastData)<-testSample[,"Insurgency"]
#' setModelNames(this.ForecastData)<-c("LMER", "SAE", "GLM")
#' }
#' 
#' @return A data object of the class 'ForecastData' with the following slots: 
#' \item{predCalibration}{An array containing the predictions of all component models for all observations in the calibration period.} 
#' \item{predTest}{An array containing the predictions of all component models for all observations in the test period.}
#' \item{outcomeCalibration}{A vector containing the true values of the dependent variable for all observations in the calibration period.} 
#' \item{outcomeTest}{A vector containing the true values of the dependent variable for all observations in the test period.}
#' \item{modelNames}{A character vector containing the names of all component models.  If no model names are specified, names will be assigned automatically.}
#' 
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2012). Improving Predictions Using Ensemble Bayesian Model Averaging. \emph{Political Analysis}. \bold{20}: 271-291.
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2015). Calibrating ensemble forecasting models with sparse data in the social sciences. \emph{International Journal of Forecasting}. \bold{31}:930–942.#'
#' 
#' @author  Michael D. Ward <\email{michael.d.ward@@duke.edu}> and Jacob M. Montgomery <\email{jacob.montgomery@@wustl.edu}> and Florian M. Hollenbach <\email{florian.hollenbach@@tamu.edu}>
#'
#' @rdname setFunctions
setGeneric("setPredCalibration<-",function(object,value){standardGeneric("setPredCalibration<-")})

#' @rdname setFunctions
#' @export
setReplaceMethod(
	f="setPredCalibration",
	signature="ForecastData",
	definition=function(object,value){
          if(is(value, "data.frame")){value <- as.matrix(value)}
          if(is(value, "matrix")){value <- array(value, dim=c(nrow(value), ncol(value), 1))}
          object@predCalibration = value
          validObject(object)
          return(object)
	}
)

#' @rdname setFunctions
setGeneric("setPredTest<-",function(object,value){standardGeneric("setPredTest<-")})

#' @rdname setFunctions
#' @export
setReplaceMethod(
	f="setPredTest",
	signature="ForecastData",
	definition=function(object,value){
          if(is(value, "data.frame")){value <- as.matrix(value)}
          if(is(value, "matrix")){value <- array(value, dim=c(nrow(value), ncol(value), 1))}
          object@predTest<- value
		validObject(object)
		return(object)
	}
)

#' @rdname setFunctions
setGeneric("setOutcomeCalibration<-",function(object,value){standardGeneric("setOutcomeCalibration<-")})

#' @rdname setFunctions
#' @export
setReplaceMethod(
	f="setOutcomeCalibration",
	signature="ForecastData",
	definition=function(object,value){
		object@outcomeCalibration <- value
		validObject(object)
		return(object)
	}
)

#' @rdname setFunctions
setGeneric("setOutcomeTest<-",function(object,value){standardGeneric("setOutcomeTest<-")})

#' @rdname setFunctions
#' @export
setReplaceMethod(
	f="setOutcomeTest",
	signature="ForecastData",
	definition=function(object,value){
		object@outcomeTest<-value
		validObject(object)
		return(object)
	}
)

#' @rdname setFunctions
setGeneric("setModelNames<-",function(object,value){standardGeneric("setModelNames<-")})

#' @rdname setFunctions
#' @export
setReplaceMethod(
	f="setModelNames",
	signature="ForecastData",
	definition=function(object,value){
		object@modelNames <-value
		colnames(object@predCalibration)<-value
		colnames(object@predTest)<-value
		validObject(object)
		return(object)
	}
)




