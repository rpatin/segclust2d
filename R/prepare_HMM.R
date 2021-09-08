#' Prepare HMM output for proper comparison plots
#'
#' \code{prepare_HMM}
#' @param data data
#' @param hmm.model hmm.model
#' @param diag.var diag.var
#' @param order.var order.var
#' @export
#' @examples 
#' \dontrun{
#' # Example taken from moveHMM package.
#' # 1. simulate data
#' # define all the arguments of simData
#' nbAnimals <- 1
#' nbStates <- 2
#' nbCovs <- 2
#' mu<-c(15,50)
#' sigma<-c(10,20)
#' angleMean <- c(pi,0)
#' kappa <- c(0.7,1.5)
#' stepPar <- c(mu,sigma)
#' anglePar <- c(angleMean,kappa)
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' zeroInflation <- FALSE
#' obsPerAnimal <- c(50,100)
#' 
#' data <- moveHMM::simData(nbAnimals=nbAnimals,nbStates=nbStates,
#'                 stepDist=stepDist,angleDist=angleDist,
#'                 stepPar=stepPar,anglePar=anglePar,nbCovs=nbCovs,
#'                 zeroInflation=zeroInflation,
#'                 obsPerAnimal=obsPerAnimal)
#' 
#' ### 2. fit the model to the simulated data
#' # define initial values for the parameters
#' mu0 <- c(20,70)
#' sigma0 <- c(10,30)
#' kappa0 <- c(1,1)
#' stepPar0 <- c(mu0,sigma0) # no zero-inflation, so no zero-mass included
#' anglePar0 <- kappa0 # the angle mean is not estimated,
#' # so only the concentration parameter is needed
#' formula <- ~cov1+cos(cov2)

#' m <- moveHMM::fitHMM(data=data,nbStates=nbStates,stepPar0=stepPar0,
#'          anglePar0=anglePar0,formula=formula,
#'          stepDist=stepDist,angleDist=angleDist,angleMean=angleMean)
#'          
#' ### 3. Transform into a segmentation-class object
#' res.hmm <- prepare_HMM(data=data, 
#' hmm.model = m, diag.var = c("step","angle"))
#' ### 4. you can now apply the same function than for segclust2d outputs
#' plot(res.hmm)
#' segmap(res.hmm)
#' }

prepare_HMM <- function(
  data, hmm.model = NULL, 
  diag.var, order.var = diag.var[1]
  ){
  outputsHMM <- 
    stat_segm_HMM(
      data = data,hmm.model = hmm.model,
      diag.var = diag.var, order.var = order.var)
  names(outputsHMM) <- c("segments","states")

  segmented <- list("data" = data,
                    "type" = "behavior",
                    "seg.type" = "HMM",
                    "outputs" = outputsHMM,
                    "model.hmm" = hmm.model,
                    "Diagnostic variables" = diag.var,
                    "Order variable" = order.var)
  class(segmented) <- "segmentation"
  return(segmented)
}

#' Get segment statistic for HMM model
#'
#' \code{stat_segm_HMM}
#' @inheritParams prepare_HMM

stat_segm_HMM <- function(data, hmm.model = NULL, diag.var, order.var = NULL){
  df.segm <- prep_segm_HMM(data,hmm.model)
  data$indice <- seq_len(nrow(data))
  df.states <- calc_stat_states(data,df.segm,diag.var,order.var)
  return(list(df.segm,df.states))
}

#' Internal function for HMM
#'
#' \code{prep_segm_HMM}
#' @export
#' @inheritParams prepare_HMM

prep_segm_HMM <- function(data,hmm.model){
  cluster <- moveHMM::viterbi(hmm.model)
  proba_states <- as.data.frame(moveHMM::stateProbs(hmm.model))
  nstates <- dim(proba_states)[2]
  colnames(proba_states) <- paste("state",1:nstates,sep="")
  df.segm <- data.frame("begin"=1:(nrow(data)-1),
                        "end"=2:(nrow(data)),
                        "state"=cluster[1:(nrow(data)-1)])
  df.segm <- cbind(df.segm,proba_states[1:(nrow(data)-1),])
  return(df.segm)
}
