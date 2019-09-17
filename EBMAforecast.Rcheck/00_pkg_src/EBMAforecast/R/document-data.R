#' Sample data Presidential Election
#'
#' This includes the data for the presidential election forecasting example in Montgomery, Hollenbach and Ward (2012). The data ranges from 1952 to 2008 and includes predictions for the six different component models included in the Ensemble model. Users may split the sample into calibration and test sample.
#'
#' The variables included in the dataset are:
#' \itemize{
#' \item\code{Campbell} Predictions of Campbell's ``Trial-Heat and Economy Model'' (Campbell 2008).
#' \item\code{Abramowitz} Predictions of Abramowitz's ``Time for Change Model'' (Abramowitz 2008).
#' \item\code{Hibbs} Predictions for the ``Bread and Peace Model'' created by Douglas Hibbs (2008).
#' \item\code{Fair} Forecasts from Fair's presidential vote share model (2010).
#' \item\code{Lewis-Beck/Tien} Predictions from the ``Jobs Model Forecast''	by Michael Lewis-Beck and Charles Tien (2008).
#' \item\code{EWT2C2} Predictions from the model in Column 2 in Table 2 by Erickson and Wlezien (2008).
#' \item\code{Actual} The true values of the dependent variable, i.e. the incumbent-party voteshare in each presidential election in the sample.
#'}
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2012). Improving Predictions Using Ensemble Bayesian Model Averaging.  \emph{Political Analysis}. \bold{20}: 271-291.
#'
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2015). Calibrating ensemble forecasting models with sparse data in the social sciences. \emph{International Journal of Forecasting}. \bold{31}:930–942.
#'
#' @references Campbell, James E. 2008. The trial-heat forecast of the 2008 presidential vote: Performance and value considerations in an open-seat election.  \emph{PS: Political Science & Politics} \bold{41}:697-701.
#'
#' @references Hibbs, Douglas A. 2000. Bread and peace voting in U.S presidential elections.  \emph{Public Choice} \bold{104}:149-180.
#'
#' @references Fair, Ray C. 2010. Presidential and Congressional vote-share equations: November 2010 update. Working Paper. Yale University.
#'
#' @references Lewis-Beck, Michael S. and Charles Tien. 2008. The job of president and the jobs model forecast: Obama for '08?  \emph{PS: Political Science & Politics} \bold{41}:687-690.
#'
#' @references Erikson, Robert S. and Christopher Wlezien. 2008. Leading economic indicators, the polls, and the presidential vote.  \emph{PS: Political Science & Politics} \bold{41}:703-707.
#'
#' @rdname presidentialForecast
#' @docType data
"presidentialForecast"


#' Sample data Insurgency Predictions
#'
#' This includes the data for the predictions of insurgencies in 29 countries for 2010. 
#' 
#' The predictions included in the dataset are:
#' \itemize{
#' \item\code{LMER} Predictions from a generalized linear mixed effects model using a logistic link function and including a randomeffects term for lagged GDP per capita and the lagged number of conflictual events involving the United States in the country of interest. 
#' \item\code{SAE} Predictions from a one model developed as part of the ICEWS project and was designed by Strategic Analysis Enterprises.
#' \item\code{GLM} Predictions from a crude logistic model that includes only population size, GDP growth (both lagged 3 months), the number of minority groups at risk in the country, and a measure of anocracy supplied in the Polity IV data set.
#'}
#' More detail about each model can be found in Mongomery et al. (2012)
#' 
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2012). Improving Predictions Using Ensemble Bayesian Model Averaging.  \emph{Political Analysis}. \bold{20}: 271-291.
#'
#' @examples \dontrun{
#' data(calibrationSample)
#' data(testSample)
#'
#' this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],
#'.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],
#' .outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
#' initW <- rep(1/3,3)
#' 
#' this.ensemble.em <- calibrateEnsemble(this.ForecastData, model="logit", tol=0.001)
#'
#' this.ensemble.gibbs <- calibrateEnsemble(this.ForecastData, model="logit", method = "gibbs")
#'}
#'
#' @rdname InsurgencyPredictions
#' @docType data
"calibrationSample"

#' @rdname InsurgencyPredictions
"testSample"

