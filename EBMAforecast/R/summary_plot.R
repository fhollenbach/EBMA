#' @include forecastData.R
NULL

#' @rdname SummarizePlot
#' @importFrom separationplot separationplot

#' @export
setClass(Class="SummaryForecastData",
         representation=representation(
           summaryData="matrix"        
         ),
         prototype=prototype(
           summaryData=matrix(NA, nrow=0, ncol=0)
         ),
)

#' Summarize and Plot Ensemble models
#'
#' These functions summarize and plot Ensemble models that have been fit previously by the user.
#'
#' @param object An object of the subclass "FDatFitLogit" or "FDatFitNormal".
#' @param x An object of class "FDatFitLogit" or "FDatFitNormal"
#' @param period The period for which the summary should be provided or plot produced, either "calibration" or "test". The default is "calibration".
#' @param fitStatistics A vector naming statistics that should be calculated.  Possible values for objects in the "FDatFitLogit" subclass include "auc", "brier", "percCorrect", "pre". Possible values for objects in the "FDatFitNormal" subclass include "rmse" and "mae." The default is for all statistics to be calculated. Additional metrics will be made available in a future release of this package.
#' @param threshold The threshold used to calculate when a "positive" prediction is made for a model.  The default is 0.5. Not used for objects of the "FDatFitNormal" subclass.
#' @param baseModel A vector containing predictions used to calculate proportional reduction of error ("pre"). The default is 0. Not used for objects of the "FDatFitNormal" subclass.
#' @param showCoefs A logical indicating whether model coefficients from the ensemble should be shown. The default is TRUE. 
#' @param subset The row names or numbers for the observations the user wishes to plot. The default is the first row.  Only implemented for the subclass "FDatFitNormal".
#' @param mainLabel A vector strings to appear at the top of each predictive posterior plot.  Only implemented for the subclass "FDatFitNormal"
#' @param xLab The label for the x-axis. Only implemented for the subclass "FDatFitNormal"
#' @param yLab The label for the y-axis.  Only implemented for the subclass "FDatFitNormal"
#' @param cols A vector containing the color for plotting the predictive PDF of each component model forecast. The default is a unique color for each PDF. Only implemented for the subclass "FDatFitNormal" 
#' @param ... Not implemented

#'
#' @return Either a plot or a data object of the class 'SummaryForecastData'. The data object has the following slots:
#' \item{summaryData}{Under the default, the function produces a matrix containing one row for each model plus one row for the EBMA forecast.  The first column is always the model weights assigned to the component models.  The second and third columns display the model parameters for the transformation of the component models.  The remaining columns are the requested fit statistics for all models, as calculated by the \code{copareModels} function.  If \code{showCoefs=FALSE}, then the model parameters will not be shown.}
#' @author  Michael D. Ward <\email{michael.d.ward@@duke.edu}> and Jacob M. Montgomery <\email{jacob.montgomery@@wustl.edu}> and Florian M. Hollenbach <\email{florian.hollenbach@@tamu.edu}>
#'
#' @references Raftery, A. E., T. Gneiting, F. Balabdaoui and M. Polakowski. (2005). Using Bayesian Model Averaging to calibrate forecast ensembles. \emph{Monthly Weather Review}. \bold{133}:1155--1174.
#' @references Greenhill, B., M.D. Ward, A. Sacks. (2011). The Separation Plot: A New Visual Method For Evaluating the Fit of Binary Data. \emph{American Journal of Political Science}.\bold{55}: 991--1002.
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2012). Improving Predictions Using Ensemble Bayesian Model Averaging. \emph{Political Analysis}. \bold{20}: 271-291.
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2015). Calibrating ensemble forecasting models with sparse data in the social sciences. \emph{International Journal of Forecasting}. \bold{31}:930–942.
#' @examples \dontrun{ data(calibrationSample)
#' data(testSample) 
#' 
#' this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],
#' .outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],
#' .outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
#' 
#' this.ensemble <- calibrateEnsemble(this.ForecastData, model="logit", tol=0.001,exp=3)
#' 
#' summary(this.ensemble, period="calibration") 
#'
#' summary(this.ensemble, period="test",showCoefs=FALSE)
#'}
#'
#' @rdname SummarizePlot
#' @export
setMethod(
  f="summary",
  signature="FDatFitLogit",
  definition=function(object,
                      period="calibration",
                      fitStatistics=c("brier", "auc", "perCorrect", "pre"),
                      threshold=.5,
                      baseModel=0,
                      showCoefs=TRUE,
                      ...){
    out <- compareModels(object, .period=period, .fitStatistics=fitStatistics, .threshold=threshold, .baseModel=baseModel)@fitStatistics
    if(showCoefs){
      coefs <- data.matrix(as.data.frame(t(plyr::aaply(object@modelParams, 1:2, function(x) {mean(x, na.rm=TRUE)}))))
      coefs <- rbind(c(NA,NA), coefs)
      out <- cbind(coefs, out)
    }
    # Adding weights column
    W <- object@modelWeights
    W <- c(NA, W)
    out <- cbind(W, out)
    
    rownames(out) <- c("EBMA", object@modelNames)
    new("SummaryForecastData", summaryData=out)
  }
)

#' @rdname SummarizePlot
#' @export
setMethod(
  f="summary",
  signature="FDatFitNormal",
  definition=function(object,
                      period="calibration",
                      fitStatistics=c("rmse", "mae"),
                      threshold=.5,
                      baseModel=0,
                      showCoefs=TRUE,
                      ...){
    
    out <- compareModels(object, .period=period, .fitStatistics=fitStatistics, .threshold=threshold, .baseModel=baseModel)@fitStatistics
    
    if(showCoefs){
      coefs <- data.matrix(as.data.frame(t(plyr::aaply(object@modelParams, 1:2, function(x) {mean(x, na.rm=TRUE)}))))
      coefs <- rbind(c(NA,NA), coefs)
      out <- cbind(coefs, out)
    }
    # Adding weights column
    W <- object@modelWeights
    W <- c(NA, W)
    out <- cbind(W, out)
    
    
    rownames(out) <- c("EBMA", object@modelNames)
    new("SummaryForecastData", summaryData=out)
  }
)

#' @rdname SummarizePlot
#' @export
setMethod(
  f="plot",
  signature="FDatFitLogit",
  definition=function(x, period="calibration",  subset=1, 
                      mainLabel="", xLab="", yLab="", cols=1, ...){ #everything behind period, is not used for logit and is a hack so we don't get a warning about unused documented functions
    nDraw=1
    numModels <- length(x@modelWeights)+1
    modelNames <- c("EBMA", x@modelNames)
    if(period=="calibration"){
      .pred <- x@predCalibration; .actual <- x@outcomeCalibration
    }
    else{
      .pred <- x@predTest; .actual <- x@outcomeTest
    }
    par(mgp=c(1, 0, 0), lend = 2, mar=c(1,0,1,0), mfrow=c(numModels, 1))
    for (i in 1:numModels){
      .miss <- is.na(.pred[,i, nDraw])
      separationplot::separationplot(pred=as.vector(.pred[!.miss,i, nDraw]), actual=as.vector(.actual[!.miss]), heading=modelNames[i], newplot=F)
    }
  }
)

#' @rdname SummarizePlot
#' @export
setMethod(
  f="plot",
  signature="FDatFitNormal",
  definition=function(x, period="calibration",  subset=1,
                      mainLabel=paste("Observation", subset), xLab="Outcome", yLab="Posterior Probability", cols=2:(length(x@modelNames)+1), ... )
  {
    if(x@method == "EM"){
      thisDraw=1
    
      if(period=="calibration"){
        .nMod <- length(x@modelNames)
        .pred <- matrix(x@predCalibration[subset,,thisDraw], ncol=.nMod+1);  colnames(.pred) <- c("EBMA", x@modelNames)
        .actual <- x@outcomeCalibration[subset]
      } else{
        .nMod <- length(x@modelNames)
        .pred <- matrix(x@predTest[subset,,thisDraw], ncol=.nMod+1); colnames(.pred) <- c("EBMA", x@modelNames)
        .actual <- x@outcomeTest
      }
      
      .sd <- sqrt(x@variance)
      
      if (length(subset)>1){
        for (j in 1:nrow(.pred)){
          .means <- .pred[j,x@modelNames]
          .miss <- is.na(.means)
          .nModThis <- sum(!.miss)
          .means <- .means[!.miss]
          
          .xMin <- min(.means)-2.5*.sd;  .xMax <- max(.means)+2.5*.sd
          .xAxis <- seq(.xMin, .xMax, by=.01);  .yAxis <- matrix(NA, .nModThis, length(.xAxis)) 
          W <- x@modelWeights[!.miss]
          for(i in 1:.nModThis){ .yAxis[i,] <- dnorm(.xAxis, mean=.means[i], sd=.sd)*W[i] }
          .totals <- colSums(.yAxis)
          plot(NULL, xlim=c(.xMin, .xMax), ylim=c(0,max(.totals)), main=mainLabel[j], xlab=xLab, ylab=yLab)
          for(i in 1:.nModThis){lines(.xAxis, .yAxis[i,], type="l", lty=2,  col=cols[i])}
          lines(.xAxis, colSums(.yAxis), lwd=2)
          rug(.means);  rug(.pred[j,"EBMA"], lwd=3)
          abline(v=.actual[j], lty=3)
        }
      } else {
        .means <- .pred[,x@modelNames]
        .miss <- is.na(.means)
        .nModThis <- sum(!.miss)
        .means <- .means[!.miss]
        .xMin <- min(.means)-2.5*.sd
        .xMax <- max(.means)+2.5*.sd
        .xAxis <- seq(.xMin, .xMax, by=.01)
        .yAxis <- matrix(NA, .nModThis, length(.xAxis)) 
        W <- x@modelWeights[!.miss]
        for(i in 1:.nModThis){.yAxis[i,] <- dnorm(.xAxis, mean=.means[i], sd=.sd)*W[i]}
        .totals <- colSums(.yAxis)
        plot(NULL, xlim=c(.xMin, .xMax), ylim=c(0,max(.totals)), main=mainLabel, xlab=xLab, ylab=yLab)
        for(i in 1:.nModThis){lines(.xAxis, .yAxis[i,], type="l", lty=2, col=cols[i])}
        lines(.xAxis, colSums(.yAxis))
        rug(.means); rug(.pred[,"EBMA"], lwd=3)
        abline(v=.actual, lty=3)
      }
    }
    if(x@method == "gibbs"){
      thisDraw <- 1
      
      if(period=="calibration"){
        .nMod <- length(x@modelNames)-1 ### not EBMA pred
        .pred <- matrix(x@predCalibration[subset,-which(colnames(x@predCalibration)=="EBMA"),thisDraw], ncol=.nMod+1);  colnames(.pred) <- c(x@modelNames)
        .actual <- x@outcomeCalibration[subset]
        .posteriorW <- x@posteriorWeights
      } else{
        .nMod <- length(x@modelNames)-1
        .pred <- matrix(x@predTest[subset,-which(colnames(x@predCalibration)=="EBMA"),thisDraw], ncol=.nMod+1); colnames(.pred) <- c(x@modelNames)
        .actual <- x@outcomeTest
        .posteriorW <- x@posteriorWeights
        }
      
      .sd <- sqrt(x@variance)
      
      if (length(subset)>1){
        for (j in 1:nrow(.pred)){
          .means <- .pred[j,x@modelNames]
          .miss <- is.na(.means)
          .nModThis <- sum(!.miss)
          .means <- .means[!.miss]
          
          .xMin <- min(.means)-2.5*.sd;  .xMax <- max(.means)+2.5*.sd
          #.xAxis <- seq(.xMin, .xMax, by=.01);  
          #.yAxis <- matrix(NA, .nModThis, length(.xAxis)) 
          .yAxis <- matrix(NA, .nModThis, 1000) ### matrix to then draw predictive density for each model
          .maxima <- rep(NA, .nModThis) ### vector to record maxima of each density
          W <- x@modelWeights[!.miss]
          for(i in 1:.nModThis){ 
            .yAxis[i,] <- rnorm(1000, mean=.means[i], sd=.sd) ### draw predictive distribution for each component model
            .maxima[i] <- max(density(.yAxis[i,])$y) ### record maximum of each density
          }
          .posteriorPred <- apply(.yAxis,2, FUN =  function(x){.posteriorW%*%as.matrix(x)})
          if(x@predType == "posteriorMedian"){.posteriorSummary <- median(.posteriorPred)}
          if(x@predType == "posteriorMean"){.posteriorSummary <- mean(.posteriorPred)}
          #.totals <- colSums(.yAxis)
          plot(NULL, xlim=c(.xMin, .xMax), ylim=c(0,max(c(density(.posteriorPred)$y,.maxima))), main=mainLabel[j], xlab=xLab, ylab=yLab)
          for(i in 1:.nModThis){lines(density(.yAxis[i,]), type="l", lty=2,  col=cols[i])}
          lines(density(.posteriorPred), lwd=2)
          rug(.means);  
          rug(.posteriorSummary, lwd=3)
          abline(v=.actual[j], lty=3)
        }
      } else {
        .means <- .pred[,x@modelNames]
        .miss <- is.na(.means)
        .nModThis <- sum(!.miss)
        .means <- .means[!.miss]
        .xMin <- min(.means)-2.5*.sd
        .xMax <- max(.means)+2.5*.sd
        .yAxis <- matrix(NA, .nModThis, 1000) ### matrix to then draw predictive density for each model
        .maxima <- rep(NA, .nModThis) ### vector to record maxima of each density
        W <- x@modelWeights[!.miss]
        for(i in 1:.nModThis){ 
          .yAxis[i,] <- rnorm(1000, mean=.means[i], sd=.sd) ### draw predictive distribution for each component model
          .maxima[i] <- max(density(.yAxis[i,])$y) ### record maximum of each density
        }
        .posteriorPred <- apply(.yAxis,2, FUN =  function(x){.posteriorW%*%as.matrix(x)})
        if(x@predType == "posteriorMedian"){.posteriorSummary <- median(.posteriorPred)}
        if(x@predType == "posteriorMean"){.posteriorSummary <- mean(.posteriorPred)}
        #.totals <- colSums(.yAxis)
        plot(NULL, xlim=c(.xMin, .xMax), ylim=c(0,max(c(density(.posteriorPred)$y,.maxima))), main=mainLabel, xlab=xLab, ylab=yLab)
        for(i in 1:.nModThis){lines(density(.yAxis[i,]), type="l", lty=2,  col=cols[i])}
        lines(density(.posteriorPred), lwd=2)
        rug(.means);  
        rug(.posteriorSummary, lwd=3)
        abline(v=.actual, lty=3)
      }
    }
  }
)





