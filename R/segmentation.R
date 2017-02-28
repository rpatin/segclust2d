#' Segmentation of movement data
#'
#' Segmentation of movement data. Actually implementation exists for
#' \code{moveHMM} and \code{segTraj} methods. The function calls either
#' \code{\link{segmentation_HMM}} or \code{\link{segmentation_picard}}
#' @param data the data.frame with the different variable
#' @param type type of model either 'hmm' or 'picard'
#' @inheritParams segmentation_hmm
#' @inheritParams segmentation_picard
#' @return  a \code{\link{segmentation-class}} object
#'
#' @examples
#' segmentation(data,diag.var=c("dist","angle"),order.var='dist',type='hmm',hmm.model=mod1.hmm)
#' @export

segmentation <- function(data, seg.var, diag.var = seg.var, order.var = dplyr::first(diag.var), type = 'picard', picard.type = NULL, scale.variable = F, ...){
  if(type == 'HMM'){
    seg <- segmentation_hmm(data, ...)
    return(seg)
  } else if (type == 'picard'){

    if(is.null(seg.var)) stop("seg.var missing for picard segmentation")

    if( length(seg.var) == 1 ){
      dat <- t(data[,rep(seg.var,2)])
    } else if ( length(seg.var) == 2 ) {
      dat <- t(data[,seg.var])
    } else {
      stop("seg.var must contains either one or two column names")
    }
    if(scale.variable) {
      dat[1,]<- scale(dat[1,])
      dat[2,]<- scale(dat[2,])
    }

    seg <- segmentation_picard(data,seg.var = seg.var, diag.var = diag.var, order.var = order.var, dat=dat, ...)

    return(seg)
  } else {
    stop("type of segmentation not recognized. Either \'HMM\' or \'picard\'")
  }
}

#' Segmentation Function for HMM
#' Wrapper for moveHMM::fitHMM
#' @param data data.frame
#' @param nbStates number of states fitted
#' @param mu0 a priori mean of dist
#' @param sigma0 a priori sd of dist
#' @param zeromass0 a priori zeromass of dist
#' @param angleMean0 a priori mean of angle
#' @param kappa0 a priori kappa of angle
#' @param angleDist type of HMM
#' @param coordNames names of coordinates in data
#' @param coord.type whether coordinates are in UTM or lat/lon
#' @param timecol column names for time in data
#' @param hmm.param parameters used for fitting the hmm model
#' @return a segmentation object

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

#' Segmentation Function for Picard/segTraj segmentation
#' @param data data.frame
#' @param seg.var names of the variables used for segmentation (either one or two names)
#' @param diag.var names of the variables on which statistics are calculated
#' @param order.var names of the variable with which states are ordered
#' @param picard.type whether mode should be segmentation only (DynProg),
#'   clustering-segmentation with a given number of class (hybrid_simultanee) or
#'   with an unknown number of class (variable_class).
#' @param scale.variable whether variables need to be scaled for segmentation
#' @param Kmax maximum number of segments
#' @param lmin minimum size of segments
#' @param nclass number of class (hybrid_simultanee)
#' @param nclass.max maximum number of class (variable_class)
#' @param sameSigma whether segment should have equal variance or not.
#'   (hybrid_simultanee & variable_class)

segmentation_picard <- function(data, seg.var = NULL, diag.var = NULL, order.var = NULL, scale.variable = F, nclass.max = NULL, Kmax = NULL, lmin = NULL, sameSigma=F, dat=NULL, picard.type='DynProg'){
  if(picard.type == 'DynProg'){
    seg <- segmentation_picard_dynprog()
  } else if ( picard.type == 'hybrid_simultanee' ) {
    seg <- segmentation_picard_hybrid(data,seg.var = seg.var, diag.var = diag.var, order.var = order.var,  dat=dat, ...)
  } else if ( picard.type == 'variable_class') {
    seg <- segmentation_picard_variable_class(data,seg.var = seg.var, diag.var = diag.var, order.var = order.var, dat=dat, ...)
  } else {
    stop("picard.type must be either \'DynProg\', \'hybrid_simultanee\' or \'variable_class\' ")
  }
  return(seg)
}

#' Segmentation Function for Picard/segTraj segmentation only mode
#' @rdname segmentation_picard

segmentation_picard_dynprog <- function(data, seg.var = NULL, diag.var = NULL, order.var = NULL, scale.variable = F, Kmax = NULL, lmin = NULL, dat=NULL){

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
                    "Order variable" = order.var,
                    "param"= list("lmin"=lmin,
                                  "Kmax"=Kmax))
  class(segmented) <- "segmentation"
  return(segmented)
}

#' Segmentation Function for Picard/segTraj clustering-segmentation mode
#' @rdname segmentation_picard

segmentation_picard_hybrid <- function(data, seg.var = NULL, diag.var = NULL, order.var = NULL, scale.variable = F, nclass = NULL, Kmax = NULL, lmin = NULL, sameSigma=F, dat=NULL){

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
                    "Order variable" = order.var,
                    "param"= list("lmin"=lmin,
                                  "Kmax"=Kmax,
                                  "nclass"=nclass))
  class(segmented) <- "segmentation"
  return(segmented)
}

#' Segmentation Function for Picard/segTraj testing both clustering-segmentation and segmentation-only
#' @rdname segmentation_picard

segmentation_picard_variable_class <- function(data, seg.var = NULL, diag.var = NULL, order.var = NULL, scale.variable = F, nclass.max = NULL, Kmax = NULL, lmin = NULL, sameSigma=F, dat=NULL){

  segmented <- list("data" = data,
                    "type" = "picard",
                    "picard.type" = "variable_class",
                    "outputs" = list(),
                    "likelihood" = NULL,
                    "picard.param" = list(),
                    "Segmented variables" = seg.var,
                    "Diagnostic variables" = diag.var,
                    "Order variable" = order.var,
                    "param"= list("lmin"=lmin,
                                  "Kmax"=Kmax,
                                  "nclass.max"=nclass.max))

  class(segmented) <- "segmentation"

  # DynProg segmentation nclass=0
  CostLoc <- segTraj::Gmean_simultanee(dat, lmin = lmin)
  res.DynProg <- segTraj::DynProg(CostLoc, Kmax)

  outputs <- lapply(1:Kmax,function(k){
    out <- stat_segm(data, diag.var, order.var, model.type='picard', picard.param = res.DynProg, picard.type = 'dynprog', picard.nseg=k)
    names(out) <- c("segments","states")
    return(out)
  })
  names(outputs) <- paste("0 class -",1:Kmax, "segments")

  likelihood <- data.frame(nseg=1:Kmax,likelihood=-res.DynProg$J.est,nclass=0)

  segmented$outputs <- c(segmented$outputs,outputs)
  segmented$likelihood <- likelihood

  if(nclass.max < 2){stop("nclass.max must be >= 2")} else {
    for(P in 2:nclass.max){
      res <- segTraj::hybrid_simultanee(dat, P = P, Kmax = Kmax, lmin = lmin, sameSigma = sameSigma)
      outputs <- lapply(P:Kmax,function(k){
        out <- stat_segm(data, diag.var, order.var, model.type='picard', picard.param = res$param[[k]], picard.type = 'hybrid')
        names(out) <- c("segments","states")
        return(out)
      })
      names(outputs) <- paste(P,"class -",P:Kmax, "segments")
      likelihood = data.frame(nseg=1:Kmax,likelihood = c(res$Linc),nclass=P)
      picard.param <- list(res$param)
      names(picard.param) <- paste(P,"class")
      segmented$likelihood <- rbind(segmented$likelihood,likelihood)
      segmented$picard.param <- c(segmented$picard.param,picard.param)
      segmented$outputs <- c(segmented$outputs,outputs)
    }
  }
  return(segmented)
}
