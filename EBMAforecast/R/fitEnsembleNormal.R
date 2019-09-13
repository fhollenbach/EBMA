#' @include forecastData.R
NULL

#' @useDynLib EBMAforecast
#' @importFrom Rcpp sourceCpp


#' @importFrom plyr alply aaply
setMethod(f="fitEnsemble",
          signature(.forecastData="ForecastDataNormal"),
          definition=function(.forecastData, tol = sqrt(.Machine$double.eps),
            maxIter=1e6,
            method="EM",
            exp=numeric(),
            useModelParams = TRUE,
            predType="posteriorMedian",
            const=0,
            W = c(),
            whichW=1)
          {
            # Creating blank store matrix
            store.W <- matrix()
            # Check if W is vector or matrix
            # Matrix
            if(is.matrix(W)){
              if(dim(W)[2] != dim(.forecastData@predCalibration)[2]){
                stop("The number of initial model weights must be of length of the number of predictive models included.")}
              for(i in 1:nrow(W)){
                if(sum(W[i,]) != 1){
                  stop("Each set of initial model weights must sum to 1.")}
              }
            }
            # Vector
            if(is.null(dim(W))){
              #check wether W is of right length and sums to 1
              if(length(W) != dim(.forecastData@predCalibration)[2] & is.null(W)==FALSE){
                stop("Vector of initial model weights must be of length of the number of predictive models included.")}
              if(sum(W) != 1 & is.null(W)==FALSE){
                stop("Vector of initial model weights must sum to 1.")}
            }

            .predictCal <- function(x){
              .rawPred <- predict(x)
              .outPred <- rep(NA, nObsCal)
              .outPred[as.numeric(names(.rawPred))] <- .rawPred
              return(.outPred)
            }

            .modelFitter <- function(preds){
              thisModel <- lm(outcomeCalibration~preds)
              return(list(thisModel, cook=cooks.distance(thisModel)))
            }

            .predictTest <- function(x, i){
              .rawPred <- predict(.models[[i]], newdata=data.frame(preds=x))
              .outPred <- rep(NA, nObsTest)
              .outPred[as.numeric(names(.rawPred))] <- .rawPred
              return(.outPred)
            }

            .ebmaMedian<-function(W, x, sdVec){
              .x <- x[!is.na(x)]
              .W <- W[!is.na(x)]
              .sdVec <- sdVec[!is.na(x)]

              ebmaCdf<-function(z, .x, .sdVec, .W){
                sum(.W*pnorm(z, mean=.x, sd=.sdVec))
              }
              low <- min(.x-6*.sdVec)
              up <- max(.x+6*.sdVec)
              out <- uniroot(function(z){ebmaCdf(z, .x=.x, .sdVec=.sdVec, .W=.W)-.5}
                             , lower = low
                             , upper = up
              )

              out$root
            }


            ##Extract data
            predCalibration <- .forecastData@predCalibration; outcomeCalibration <- .forecastData@outcomeCalibration
            predTest <- .forecastData@predTest; outcomeTest <- .forecastData@outcomeTest
            .testPeriod <- length(predTest)>0
            modelNames <- .forecastData@modelNames

            ## Set constants
            nMod <-  ncol(predCalibration); nDraws <- dim(predCalibration)[3]
            nObsCal <- nrow(predCalibration); nObsTest <- nrow(predTest)
            ZERO<-1e-4

            ## Fit Models
            if(useModelParams==TRUE){
              .models <- plyr::alply(predCalibration, 2:3, .fun=.modelFitter)
              for(i in 1:nMod){
                if(any(unname(.models[[i]][[2]]) > 0.5)){
                  cat("Problematic Cook's Distances (> 0.5) \n", "Model", names(.models[i]), ":",
                      which(unname(.models[[i]][[2]]) > 0.5), "\n")
                  warning("Problematic Cook's Distances (> 0.5), see above output (under 'this.ensemble').")
                }
                .models[[i]] <- .models[[i]][[1]]
              }
            }

            ## Extract needed info
            if(nDraws==1 & useModelParams==TRUE){
              predCalibrationAdj <- aperm(array(plyr::laply(.models, .predictCal), dim=c(nMod, nObsCal, nDraws)), c(2,1,3))
              modelParams <- aperm(array(plyr::laply(.models, coefficients), dim=c(nMod, 2, nDraws)), c(2,1,3))
            }
            if(nDraws>1 & useModelParams==TRUE){ # This code is in development for exchangeability
              predCalibrationAdj <- aperm(plyr::aaply(.models, 1:2, .predictCal), c(3,1,2))
              modelParams <- aperm(plyr::aaply(.models, 1:2, coefficients), c(3,1,2))
            }
            if(useModelParams==FALSE){
              predCalibrationAdj <- predCalibration
              modelParams <- array(c(0,1), dim=c(2,nMod,nDraws))
            }
            calResiduals <- outcomeCalibration-predCalibrationAdj
            calResiduals2 <- calResiduals^2

            dimnames(modelParams) <- list(c("Constant", "Predictor"), modelNames, 1:nDraws)
            dimnames(calResiduals) <- dimnames(calResiduals2) <-dimnames(predCalibrationAdj) <- list(1:nObsCal, modelNames, 1:nDraws)

            sigma2<-1

            if(is.matrix(W)){
              W2 <- W[whichW,]
              # Matrix to store all posterior weights in
              store.W <- matrix(data=NA, nrow=dim(W)[1], ncol=dim(W)[2])
              colnames(store.W) <- modelNames
              for(i in 1:dim(W)[1]){
                out  = emNorm(outcomeCalibration, matrix(predCalibrationAdj[,,1],ncol=nMod),matrix(calResiduals2[,,1],ncol=nMod), W, tol, maxIter, const, sigma2)
                if (out$Iterations==maxIter){print("WARNING: Maximum iterations reached")}
                W <- out$W*rowSums(!colSums(predCalibration, na.rm=T)==0); names(W) <- modelNames
                sigma2 = out$Sigma2
                LL = out$LL
                iter = out$Iterations
                store.W[i,] <- W
              }
              # Calculating mean absolute difference of posterior weights
              store.MAD <- matrix(data=NA, nrow=1, ncol=dim(store.W)[2])
              colnames(store.MAD) <- modelNames
              for(i in 1:dim(store.W)[2]){
                dif <- abs(store.W[1, i] - store.W[2, i])
                out <- (mean(dif, na.rm = TRUE))
                store.MAD[, i] <- out
              }
              W <- W2
              # Error if any mean absolute difference of posterior weights > 0.0001
              if(any(store.MAD > 0.0001)){
                warning("The mean absolute difference between the sets of posterior weights is above 0.0001.
                  The posterior EBMA prediction is only based on the first set of weights.")
                }
              }



            ## Set initial values for parameters
            if(is.null(W)){
            W <- rep(1/(nMod), nMod) ; names(W) <- modelNames
            }


            ### call to rcpp for em
            if(is.null(dim(W))){
			        out  = emNorm(outcomeCalibration, matrix(predCalibrationAdj[,,1],ncol=nMod),matrix(calResiduals2[,,1],ncol=nMod), W, tol, maxIter, const, sigma2)
              if (out$Iterations==maxIter){print("WARNING: Maximum iterations reached")}
              W <- out$W*rowSums(!colSums(predCalibration, na.rm=T)==0); names(W) <- modelNames
              sigma2 = out$Sigma2
              LL = out$LL
              iter = out$Iterations
            }

            ## Merge the EBMA forecasts for the calibration sample onto the predCalibration matrix
            .flatPreds <- plyr::aaply(predCalibrationAdj, c(1,2), function(x) {mean(x, na.rm=TRUE)})
            .sdVec <- rep(sqrt(sigma2), nMod)

            if (predType=="posteriorMean"){
              bmaPred <- array(plyr::aaply(.flatPreds, 1, function(x) {sum(x* W, na.rm=TRUE)}), dim=c(nObsCal, 1,nDraws))
              bmaPred <-  bmaPred/array(t(W%*%t(1*!is.na(.flatPreds))), dim=c(nObsCal, 1, nDraws))
              bmaPred[,,-1] <- NA
            }

            if (predType=="posteriorMedian"){
              .altQBMAnormal <- function(x){
                .x <- x[!is.na(x)]
                .W <- W[!is.na(x)]
                ..sdVec <- .sdVec[!is.na(x)]
                .ebmaMedian(.W, .x, ..sdVec)
              }
             bmaPred <- array(plyr::aaply(.flatPreds, 1, .altQBMAnormal),  dim=c(nObsCal, 1,nDraws))
             bmaPred[,,-1] <- NA
            }
            cal <- abind::abind(bmaPred, .forecastData@predCalibration, along=2); colnames(cal) <- c("EBMA", modelNames)



            if(.testPeriod){
              if(useModelParams==TRUE){
                predTestAdj <- array(NA, dim=c(nObsTest, nMod, nDraws))
                for (k in 1:nMod){
                 for (j in 1:nDraws){
                   predTestAdj[,k,j] <- .predictTest(predTest[,k,j], i=k)
                   }
                }
              }
              if(useModelParams==FALSE){predTestAdj <- predTest}
              .flatPredsTest <- matrix(plyr::aaply(predTestAdj, c(1,2), function(x) {mean(x, na.rm=TRUE)}), ncol=nMod)

              if (predType=="posteriorMean"){
                bmaPredTest <-array(plyr::aaply(.flatPredsTest, 1, function(x) {sum(x* W, na.rm=TRUE)}), dim=c(nObsTest, 1,nDraws))
                bmaPredTest <-  bmaPredTest/array(t(W%*%t(1*!is.na(.flatPredsTest))), dim=c(nObsTest, 1, nDraws))
                bmaPredTest[,,-1] <- NA
              }

              if (predType=="posteriorMedian"){
                .altQBMAnormal <- function(x){
                  .x <- x[!is.na(x)]
                  .W <- W[!is.na(x)]
                  ..sdVec <- .sdVec[!is.na(x)]
                  .ebmaMedian( .W, .x, ..sdVec)
                }
                bmaPredTest <- array(plyr::aaply(.flatPredsTest, 1, .altQBMAnormal),  dim=c(nObsTest, 1,nDraws))
                bmaPredTest[,,-1] <- NA
              }

              test <- abind::abind(bmaPredTest, .forecastData@predTest, along=2);  colnames(test) <- c("EBMA", modelNames)
            }
            if(!.testPeriod){{test <- .forecastData@predTest}}
            if(useModelParams==FALSE){.models = list()}

            new("FDatFitNormal",
                predCalibration=cal,
                outcomeCalibration=outcomeCalibration,
                predTest=test,
                outcomeTest=.forecastData@outcomeTest,
                modelNames=modelNames,
                modelWeights=W,
                useModelParams = useModelParams,
                modelParams=modelParams,
                variance=sigma2,
                logLik=LL,
                exp=exp,
                tol=tol,
                maxIter=maxIter,
                predType=predType,
                method=method,
                iter=iter,
                model="normal",
                modelResults = .models,
                call=match.call(),
                posteriorWeights=store.W
                )
          }
          )

