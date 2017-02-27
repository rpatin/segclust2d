#' Segmentation of movement data
#'
#' Segmentation of movement data. Actually implementation exists for
#' \code{moveHMM} and \code{segTraj} methods. The function calls either
#' \code{\link{segmentation_HMM}} or \code{\link{segmentation_picard}}
#' @param data the data.frame with the different variable
#' @param seg.var names of the variables used for segmentation (if picard)
#' @param diag.var names of the variables on which statistics are calculated
#' @param order.var names of the variable with which states are ordered
#' @param type type of model either 'hmm' or 'picard'
#' @param ... other arguments passed to segmentation_HMM and segmentation_picard
#' @return  a \code{\link{segmentation}} object
#'
#' @examples
#' segmentation(data,diag.var=c("dist","angle"),order.var='dist',type='hmm',hmm.model=mod1.hmm)
#' @export

segmentation <- function(data, seg.var, diag.var = seg.var, order.var = dplyr::first(diag.var), type = 'picard', ...){
  if(type == 'HMM'){
    seg <- segmentation_hmm(data, ...)
    return(seg)
  } else if (type == 'picard'){
    seg <- segmentation_picard(data,seg.var = seg.var, diag.var = diag.var, order.var = order.var, ...)
    return(seg)
  } else {
    stop("type of segmentation not recognized. Either \'HMM\' or \'picard\'")
  }
}

#' Segmentation Function for HMM
#' @param hmm.param parameters used for fitting the hmm model

segmentation_hmm <- function(data, nbStates, mu0, sigma0, zeromass0, angleMean0, kappa0, stepDist = "gamma", angleDist = "vm", coordNames = c("x","y"), coord.type =c("UTM"),timecol='expectTime'){

  evalstr <- paste("data2=dplyr::select(data,",timecol,",",coordNames[1],",",coordNames[2],")",sep="")
  eval(parse(text=evalstr))

  traj.hmm <- moveHMM::prepData(data2,type=coord.type,coordNames = coordNames)

  stepPar0 <- c(mu0,sigma0,zeromass0)
  anglePar0 <- c(angleMean0,kappa0)

  mod.hmm <- moveHMM::fitHMM(data=traj.hmm, nbStates = nbStates, stepDist = stepDist, stepPar0=stepPar0, angleDist = angleDist, anglePar0 = anglePar0)
  outputsHMM <- stat_segm(data,diag.var = diag.var, order.var = order.var, model.type='hmm',hmm.model = mod.hmm)
  names(outputsHMM) <- c("segments","states")

  segmented <- list("data" = data,
                    "type" = "HMM",
                    "outputs" = outputsHMM,
                    "model.hmm" = mod.hmm,
                    "Diagnostic variables" = diag.var,
                    "Order variable" = order.var)
  class(segmented) <- "segmentation"
  return(segmented)
}

#' Segmentation Function for Picard/segTraj
#' @param picard.type either 'hybrid' or 'dynprog'
#' @param scale.variable if variable needs to be scaled for segmentation
#' @param nclass number of class for hybrid_simultanee
#' @param Kmax maximum number of segments
#' @param lmin minimum size of segments
#' @param sameSigma whether segment should have equal variance or not.

segmentation_picard <- function(data, seg.var = NULL, diag.var = NULL, order.var = NULL, picard.type = 'DynProg', scale.variable = F, nclass = NULL, Kmax = NULL, lmin = NULL, sameSigma=F){

  if(is.null(seg.var)) stop("seg.var missing for picard segmentation")

  if( length(seg.var) == 1 ){
    dat <- t(data[,rep(seg.var,2)])
  } else if ( length(seg.var) == 2 ) {
    dat <- t(data[,seg.var])
  } else {
    stop("seg.var must contains either one or two column names")
  }

  if ( picard.type == 'DynProg' ) {
    CostLoc <- segTraj::Gmean_simultanee(dat, lmin = lmin)
    res.DynProg <- segTraj::DynProg(CostLoc, Kmax)

    outputs <- lapply(1:Kmax,function(k){
      out <- stat_segm(data, diag.var, order.var, model.type='picard', picard.param = res.DynProg, picard.type = 'dynprog', picard.nseg=k)
      names(out) <- c("segments","states")
      return(out)
      })
    names(outputs) <- paste(1:Kmax, "segments")
    segmented <- list("data" = data,
                      "type" = "picard",
                      "picard.type" = "DynProg",
                      "outputs" = outputs,
                      "likelihood" = data.frame(nseg=1:Kmax,likelihood=-res.DynProg$J.est),
                      "Segmented variables" = seg.var,
                      "Diagnostic variables" = diag.var,
                      "Order variable" = order.var)
    class(segmented) <- "segmentation"
    return(segmented)
  } else if ( picard.type == 'hybrid_simultanee' ) {
    res <- segTraj::hybrid_simultanee(dat, P = nclass, Kmax = Kmax, lmin = lmin, sameSigma = sameSigma)
    outputs <- lapply(nclass:Kmax,function(k){
      out <- stat_segm(data, diag.var, order.var, model.type='picard', picard.param = res$param[[k]], picard.type = 'hybrid')
      names(out) <- c("segments","states")
      return(out)
    })
    names(outputs) <- paste(nclass:Kmax, "segments")
    segmented <- list("data" = data,
                      "type" = "picard",
                      "picard.type" = "hybrid_simultanee",
                      "outputs" = outputs,
                      "likelihood" = data.frame(nseg=1:Kmax,likelihood = c(res$Linc)),
                      "picard.param" = res$param,
                      "Segmented variables" = seg.var,
                      "Diagnostic variables" = diag.var,
                      "Order variable" = order.var)
    class(segmented) <- "segmentation"
    return(segmented)

  } else {
    stop("picard.type must be either \'DynProg\' or \'hybrid_simultanee\'")
  }
}
