#' @include forecastData.R
NULL

#' @useDynLib EBMAforecast
#' @importFrom Rcpp sourceCpp


#' @importFrom plyr alply aaply
setMethod(f="fitEnsemble",
          signature(.forecastData="ForecastDataNormal"),
          definition=function(.forecastData,
            tol = sqrt(.Machine$double.eps),
            maxIter=1e6,
            method="EM",
            exp=numeric(),
            useModelParams = TRUE,
            predType="posteriorMedian",
            const=0,
            W = rep(1/dim(.forecastData@predCalibration)[2],dim(.forecastData@predCalibration)[2]),
            iterations= 40000,
            burns = 20000,
            thinning = 20,
            modelPriors = rep(1/dim(.forecastData@predCalibration)[2],dim(.forecastData@predCalibration)[2])
            )
          {
            if(iterations < (burns)){
              stop("Number of iterations is smaller than the burnin. Increase the number of iterations or decrease the burnin.")
            }
            if(method == "gibbs" & any(is.na(.forecastData@predCalibration))){
              stop("Missing values in the calibration set are currently only allowed with EM estimation.")
            }
            if(method == "gibbs"){
              cat("Model weights estimated using gibbs sampling")
            }
            if(method == "EM"){
              cat("Model weights estimated using EM algorithm")
            }
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
            # # Vector
            if(is.null(dim(W))){
              #check wether W is of right length and sums to 1
              if(length(W) != dim(.forecastData@predCalibration)[2] & is.null(W)==FALSE){
                stop("Vector of initial model weights must be of length of the number of predictive models included.")}
              if(sum(W) != 1 & is.null(W)==FALSE){
                stop("Vector of initial model weights must sum to 1.")}
            }

            ##Extract data
            predCalibration <- .forecastData@predCalibration; 
            outcomeCalibration <- .forecastData@outcomeCalibration
            predTest <- .forecastData@predTest; 
            outcomeTest <- .forecastData@outcomeTest
            .testPeriod <- length(predTest)>0
            modelNames <- .forecastData@modelNames
            
            ## Set constants
            nMod <-  ncol(predCalibration); 
            nDraws <- dim(predCalibration)[3]
            nObsCal <- nrow(predCalibration); 
            nObsTest <- nrow(predTest)
            ZERO<-1e-4
            
            
            ## Set initial values for parameters
            if(is.null(W)){
              W <- rep(1/(nMod), nMod) 
              names(W) <- modelNames
            }
            
            
            ## Set model priors if unspecified
            if(is.null(modelPriors)){
              modelPriors <- rep(1, 1/nMod)
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


           

            ## Fit Models
            if(useModelParams==TRUE){
              .models <- plyr::alply(predCalibration, 2:3, .fun=.modelFitter)
              for(i in 1:nMod){
                if(any(unname(.models[[i]][[2]]) > 0.5)){
                  cat("WARNING: Problematic Cook's Distances (> 0.5) \n", "Model", names(.models[i]), ":",
                      which(unname(.models[[i]][[2]]) > 0.5), "\n")
                  #warning("Problematic Cook's Distances (> 0.5), see above output (under 'this.ensemble').")
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
           
            #### code block for method EM
            if(method == "EM"){
              postPredCal <- postPredTest <- W.mat <-  matrix() ### empty slots only used in gibbs
              Sigma.mat <- numeric()
              ### if W is not a matrix, call to rcpp for em
              if(is.null(dim(W))){
                out  = emNorm(outcomeCalibration, matrix(predCalibrationAdj[,,1],ncol=nMod),matrix(calResiduals2[,,1],ncol=nMod), W, tol, maxIter, const, sigma2)
                if (out$Iterations==maxIter){print("WARNING: Maximum iterations reached")}
                W <- out$W*rowSums(!colSums(predCalibration, na.rm=T)==0); names(W) <- modelNames
                sigma2 = out$Sigma2
                LL = out$LL
                iter = out$Iterations
              }
              # This is calibration based on all sets of starting weights
              # Calibrating for as many times as are different weights if is a matrix of weights input,
              # checking to see if the weights significantly differ
              if(is.matrix(W)){
                # Matrix to store all posterior weights in
                store.W <- matrix(data=NA, nrow=dim(W)[1], ncol=dim(W)[2])
                store.LL <- rep(NA, dim(W)[1])
                store.iter <- rep(NA, dim(W)[1])
                colnames(store.W) <- modelNames
                
                for(i in 1:dim(W)[1]){
                  vectorW <- W[i,]
                  
                  out  = emNorm(outcomeCalibration, matrix(predCalibrationAdj[,,1],ncol=nMod),matrix(calResiduals2[,,1],ncol=nMod), vectorW, tol, maxIter, const, sigma2)
                  if (out$Iterations==maxIter){print("WARNING: Maximum iterations reached")}
                  vectorW <- out$W*rowSums(!colSums(predCalibration, na.rm=T)==0); 
                  names(vectorW) <- modelNames
                  sigma2 = out$Sigma2
                  store.LL[i] = out$LL
                  store.iter[i] = out$Iterations
                  store.W[i,] <- vectorW
                }
                # Calculating mean absolute difference of posterior weights
                combs <- gtools::combinations(dim(W)[1], 2)
                store.MAD <- rep(NA,dim(combs)[1])
                for(i in 1:dim(combs)[1]){
                  store.MAD[i] <- mean(abs(store.W[combs[i,1], ] - store.W[combs[i,2], ]), na.rm = T)
                }
                # Warning if any mean absolute difference of posterior weights > 0.0001
                if(any(store.MAD > 0.0001)){
                  cat("WARNING: The mean absolute difference between the sets of posterior weights is above 0.0001. \n
                    The posterior EBMA prediction is only based on the last set of weights.")
                  }
                #### Use weights from last row of starting weights matrix 
                W <- store.W[dim(W)[1],]
                LL = store.LL[dim(W)[1]]
                iter = store.iter[dim(W)[1]]
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
              
              if(sum(W)<=.99 || sum(W)>1.01){
                cat("WARNING: Model weights do not sum to approximately one. Something might be wrong.")
              }
            }

            if(method == "gibbs"){
              LL <-  iter <- numeric() ### empty slots, only used for EM
              
              ### if W is not a matrix, call to rcpp for em
              if(is.null(dim(W))){
                out  = GibbsNormal(outcomeCalibration, matrix(predCalibrationAdj[,,1],ncol=nMod), W, modelPriors, sigma2, iterations, burns, thinning)
                W.mat <- out$W
                Sigma.mat = out$Sigma
              }
              # This is calibration based on all sets of starting weights
              # Calibrating for as many times as are different weights if is a matrix of weights input,
              # checking to see if the weights significantly differ
              if(is.matrix(W)){
                # Matrix to store all posterior weights in
                store.W.median <- matrix(data=NA, nrow=dim(W)[1], ncol=dim(W)[2])
                colnames(store.W.median) <- modelNames
                
                for(i in 1:dim(W)[1]){
                  vectorW <- W[i,]
                  out  = GibbsNormal(outcomeCalibration, matrix(predCalibrationAdj[,,1],ncol=nMod), vectorW, alpha = modelPriors, sigma = sigma2, iterations, burns, thinning)
                  W.mat <- out$W
                  store.W.median[i,] <- apply(W.mat, 2, FUN=median)
                  }
                ### save posterior weights based on last set of initial weights
                W.mat <- W.mat
                Sigma.mat <- out$Sigma
                
                # Calculating mean absolute difference of posterior weights
                combs <- gtools::combinations(dim(W)[1], 2)
                store.MAD <- rep(NA,dim(combs)[1])
                for(i in 1:dim(combs)[1]){
                  store.MAD[i] <- mean(abs(store.W.median[combs[i,1], ] - store.W.median[combs[i,2], ]), na.rm = T)
                  }
                # Error if any mean absolute difference of posterior weights > 0.0001
                if(any(store.MAD > 0.0001)){
                  cat("WARNING: The mean absolute difference between the sets of median posterior weights is above 0.0001. \n
                          The posterior EBMA prediction is only based on the last set of weights.")
                  }
              }
              
              ## Merge the EBMA forecasts for the calibration sample onto the predCalibration matrix
              .flatPreds <- plyr::aaply(predCalibrationAdj, c(1,2), function(x) {mean(x, na.rm=TRUE)})
              postPredCal <- matrix(data=NA, nrow=dim(predCalibrationAdj)[1], ncol=dim(W.mat)[1])
              for(i in 1:dim(W.mat)[1]){
                bmaPred <- array(plyr::aaply(.flatPreds, 1, function(x) {sum(x* W.mat[i,], na.rm=TRUE)}), dim=c(nObsCal, 1,nDraws))
                bmaPred <-  bmaPred/array(t(W.mat[i,]%*%t(1*!is.na(.flatPreds))), dim=c(nObsCal, 1, nDraws))
                bmaPred[,,-1] <- NA
                postPredCal[,i] <- bmaPred[,1,]
              }
              ### median or mean weights for results (depending on predType) and prediction
              if(predType == "posteriorMedian"){
                cat("Predictive performance statistics and vector of model weights based on posterior median.")
                W <- apply(W.mat, 2, FUN=median)
                sigma2 <- median(Sigma.mat)
                bmaPred[,1,] <- apply(postPredCal, 1, FUN=median)
              }
              if(predType == "posteriorMean"){
                cat("Predictive performance statistics and vector of model weights based on posterior mean.")
                W <- apply(W.mat, 2, FUN=mean)
                sigma2 <- mean(Sigma.mat)
                bmaPred[,1,] <- apply(postPredCal, 1, FUN=mean)
              }

              cal <- abind::abind(bmaPred, .forecastData@predCalibration, along=2); colnames(cal) <- c("EBMA", modelNames)
              if(sum(apply(W.mat, 2, mean))<=.99 || sum(apply(W.mat, 2, mean))>1.01){
                cat("WARNING: The mean posterior model weights do not sum to approximately one. Something might be wrong.")
              }
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
              
              if(method == "EM"){
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
              if(method == "gibbs"){
                postPredTest <- matrix(data=NA, nrow=dim(predTestAdj)[1], ncol=dim(W.mat)[1])
                for(i in 1:dim(W.mat)[1]){
                  bmaPredTest <-array(plyr::aaply(.flatPredsTest, 1, function(x) {sum(x* W.mat[i,], na.rm=TRUE)}), dim=c(nObsTest, 1,nDraws))
                  bmaPredTest <-  bmaPredTest/array(t(W.mat[i,]%*%t(1*!is.na(.flatPredsTest))), dim=c(nObsTest, 1, nDraws))
                  bmaPredTest[,,-1] <- NA
                  postPredTest[,i] <- bmaPredTest[,1,]
                }
                if(predType == "posteriorMean"){
                  bmaPredTest[,1,] <- apply(postPredTest, 1, FUN=mean)
                }
                if(predType == "posteriorMedian"){
                  bmaPredTest[,1,] <- apply(postPredTest, 1, FUN=median)
                }
                test <- abind::abind(bmaPredTest, .forecastData@predTest, along=2);  colnames(test) <- c("EBMA", modelNames)
              }

              test <- abind::abind(bmaPredTest, .forecastData@predTest, along=2);  colnames(test) <- c("EBMA", modelNames)

            }
            if(!.testPeriod){{test <- .forecastData@predTest}
              {postPredTest <- matrix()}
              }
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
                posteriorWeights = W.mat,
                posteriorSigma = Sigma.mat,
                posteriorPredCalibration = postPredCal,
                posteriorPredTest = postPredTest,
                call=match.call()
                )
          }
          )

