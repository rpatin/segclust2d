#' Segmentation/Clustering of movement data - Generic function
#'
#' Joint Segmentation/Clustering of movement data. Method available for data.frame, move and ltraj object
#' @param x data used for segmentation. Supported : data.frame, 2-columns matrix, move object, ltraj object
#' @param type type of segmentation. Either "home-range" or "behavior"
#' @param seg.var for behavioral segmentation : names of the variables used for segmentation (either one or two names)
#' @param diag.var for behavioral segmentation : names of the variables on which statistics are calculated
#' @param order.var for behavioral segmentation : names of the variable with which states are ordered
#' @param scale.variable for behavioral segmentation : whether variables need to be scaled for segmentation
#' @param Kmax maximum number of segments. Default to 10.
#' @param lmin minimum size of segments. Default to length of time series/Kmax/2.
#' @inheritParams segmentation.data.frame
#' @inheritParams segmentation.move
#' @inheritParams segmentation.ltraj
#' @return  a \code{\link{segmentation-class}} object
#'
#' @examples
#' segmentation(data,diag.var=c("dist","angle"),order.var='dist',type='hmm',hmm.model=mod1.hmm)
#' @export

segclust <- function (x, ...) {
  UseMethod("segclust", x)
}



#' Segmentation/Clustering function for data.frame
#' @param coord.names for home-range segmentation and data.frame, names of coordinates. Default x and y. Not needed for move and ltraj objects
#' @param ncluster number or list of cluster to be chosen.
#' @rdname segclust
#' @export

segclust.data.frame <- function(x, Kmax = 10, lmin = Kmax/2, ncluster = 2, type = "home-range", scale.variable = F, seg.var = NULL, diag.var = seg.var, order.var = seg.var[1], coord.names = c("x","y"),S=0.75,sameSigma = F){

  if(type == "home-range"){
    dat <- t(x[,coord.names])
    seg.var = coord.names
    diag.var = coord.names
    order.var = coord.names[1]
  } else if ( type == "behavior" ){
    if(is.null(seg.var)) stop("seg.var missing for behavioral segmentation")
    if( length(seg.var) == 1 ){
      dat <- t(x[,rep(seg.var,2)])
    } else if ( length(seg.var) == 2 ) {
      dat <- t(x[,seg.var])
    } else {
      stop("seg.var must contains either one or two column names")
    }
  } else {
    stop("type must be either home-range or behavior")
  }

  segmented <- segclust_internal(x, seg.var = seg.var, diag.var = diag.var, order.var = order.var, scale.variable = scale.variable, Kmax = Kmax, ncluster = ncluster, lmin = lmin, dat=dat,data.type = "data.frame",S=S,type=type,sameSigma = sameSigma)
  return(segmented)
}

#' Segmentation/Clustering function for move objects
#' @rdname segclust
#' @export


segclust.Move <- function(x, Kmax = 10, lmin = Kmax/2, ncluster = 2, type = "home-range", scale.variable = F, seg.var = NULL, diag.var = seg.var, order.var = seg.var[1], coord.names = c("x","y"),S=0.75,sameSigma = F){

  if(type == "home-range"){
    dat <- t(coordinates(x))
    seg.var = coord.names
    diag.var = coord.names
    order.var = coord.names[1]
    x.df = x@data
    x.df[,coord.names[1]] <- dat[1,]
    x.df[,coord.names[2]] <- dat[2,]
  } else if ( type == "behavior" ){
    x.df = x@data
    if(is.null(seg.var)) stop("seg.var missing for behavioral segmentation")
    if( length(seg.var) == 1 ){
      dat <- t(x.df[,rep(seg.var,2)])
    } else if ( length(seg.var) == 2 ) {
      dat <- t(x.df[,seg.var])
    } else {
      stop("seg.var must contains either one or two column names")
    }
  } else {
    stop("type must be either home-range or behavior")
  }

  segmented <- segclust_internal(x.df, seg.var = seg.var, diag.var = diag.var, order.var = order.var, scale.variable = scale.variable, Kmax = Kmax, ncluster = ncluster, lmin = lmin, dat=dat,data.type = "data.frame",S=S,type=type,sameSigma = sameSigma)
  return(segmented)
}

#' Segmentation/Clustering function for ltraj objects
#' @rdname segclust
#' @export


segclust.ltraj <- function(x, Kmax = 10, lmin = Kmax/2, ncluster = 2, type = "home-range", scale.variable = F, seg.var = NULL, diag.var = seg.var, order.var = seg.var[1], coord.names = c("x","y"),S=0.75,sameSigma = F){

  if(type == "home-range"){
    tmp <- x[[1]]
    dat <- t(cbind(tmp$x,tmp$y))
    if(any(is.na(tmp$x))) stop("Please filter NA from ltraj object")
    seg.var = coord.names
    diag.var = coord.names
    order.var = coord.names[1]
    x.df =  adehabitatLT::infolocs(x)[[1]]
    x.df[,coord.names[1]] <- dat[1,]
    x.df[,coord.names[2]] <- dat[2,]
  } else if ( type == "behavior" ){
    x.df = adehabitatLT::infolocs(x)[[1]]
    if(is.null(seg.var)) stop("seg.var missing for behavioral segmentation")
    if( length(seg.var) == 1 ){
      dat <- t(x.df[,rep(seg.var,2)])
    } else if ( length(seg.var) == 2 ) {
      dat <- t(x.df[,seg.var])
    } else {
      stop("seg.var must contains either one or two column names")
    }
  } else {
    stop("type must be either home-range or behavior")
  }

  segmented <- segclust_internal(x.df, seg.var = seg.var, diag.var = diag.var, order.var = order.var, scale.variable = scale.variable, Kmax = Kmax, ncluster = ncluster, lmin = lmin, dat=dat,data.type = "data.frame",S=S,type=type,sameSigma = sameSigma)
  return(segmented)
}

#' Internal segmentation/clustering function

segclust_internal <- function(x, seg.var = NULL, diag.var = NULL, order.var = NULL, scale.variable = scale.variable, Kmax = NULL, ncluster = NULL, lmin = NULL, dat=NULL, data.type = NULL,S=NULL,type=NULL,sameSigma = NULL){

  if(scale.variable) {
    dat[1,]<- scale(dat[1,])
    dat[2,]<- scale(dat[2,])
  }


  segmented <- list("data" = x,
                    "type" = type,
                    "seg.type" = "segclust",
                    "outputs" = list(),
                    "likelihood" = NULL,
                    "picard.param" = list(),
                    "Segmented variables" = seg.var,
                    "Diagnostic variables" = diag.var,
                    "Order variable" = order.var,
                    "Kopt.BIC" = rep(NA,max(ncluster)),
                    "ncluster.BIC" = NULL,
                    "BIC" = NULL,
                    "param"= list("lmin"=lmin,
                                  "Kmax"=Kmax,
                                  "ncluster"=ncluster))
  class(segmented) <- "segmentation"

  # DynProg segmentation nclass=0 - Initialisation
  CostLoc <- Gmean_simultanee(dat, lmin = lmin,sameVar = sameSigma)
  res.DynProg <- DynProg(CostLoc, Kmax)
  # CostLoc <- segTraj::Gmean_simultanee(dat, lmin = lmin)
  # res.DynProg <- segTraj::DynProg(CostLoc, Kmax)

  outputs <- lapply(1:Kmax,function(k){
    out <- stat_segm(x, diag.var, order.var, param = res.DynProg, nseg=k, seg.type = "segmentation")
    names(out) <- c("segments","states")
    return(out)
  })
  names(outputs) <- paste("0 class -",1:Kmax, "segments")

  likelihood <- data.frame(nseg=1:Kmax,likelihood=-res.DynProg$J.est,nclass=0)
  # dfBIC <- calc_BIC(likelihood,ncluster = 1:Kmax, nseg = 1:Kmax)
  # dfBIC$ncluster = 0
  segmented$outputs <- c(segmented$outputs,outputs)
  segmented$likelihood <- likelihood

  for(P in ncluster){
      # res <- segTraj::hybrid_simultanee(dat, P = P, Kmax = Kmax, lmin = lmin, sameSigma = sameSigma)
      res <- hybrid_simultanee(dat, P = P, Kmax = Kmax, lmin = lmin, sameSigma = sameSigma)
      outputs <- lapply(P:Kmax,function(k){
        out <- stat_segm(x, diag.var, order.var, param = res$param[[k]], seg.type = 'segclust')
        names(out) <- c("segments","states")
        return(out)
      })
      names(outputs) <- paste(P,"class -",P:Kmax, "segments")
      likelihood = data.frame(nseg=1:Kmax,likelihood = c(res$Linc),nclass=P)
      tmpBIC <- calc_BIC(likelihood$likelihood,ncluster = P, nseg = 1:Kmax, n = dim(x)[2])
      param <- list(res$param)
      names(param) <- paste(P,"class")
      segmented$likelihood <- rbind(segmented$likelihood,likelihood)
      segmented$BIC <- rbind(segmented$BIC,tmpBIC)
      segmented$param <- c(segmented$param,param)
      segmented$outputs <- c(segmented$outputs,outputs)
      segmented$Kopt.BIC[P] <- which.max(tmpBIC$BIC)
  }
  tmp <- dplyr::filter(segmented$BIC,ncluster > 0)
  segmented$ncluster.BIC <- tmp[which.max(tmp$BIC),"ncluster"]

  return(segmented)
}