#' Segmentation/Clustering of movement data - Generic function
#'
#' Joint Segmentation/Clustering of movement data. Method available for
#' data.frame, move and ltraj objects. The algorithm finds the optimal
#' segmentation for a given number of cluster and segments using an iterated
#' alternation of a Dynamic Programming algorithm and an
#' Expectation-Maximization algorithm. Among the different segmentation found,
#' the best one can be chosen using the maximum of a BIC penalized likelihood.
#'
#' @param ncluster number of cluster into which segments should be grouped. Can
#'   be a vector if one want to test several number of clusters.
#' @inheritParams segmentation
#' @inheritParams segclust.data.frame
#' @inheritParams segclust.Move
#' @inheritParams segclust.ltraj
#' @param ... additional parameters given to \code{\link{segclust_internal}}.
#' @return  a \code{\link{segmentation-class}} object
#'
#' @examples
#' \dontrun{segclust(data, Kmax = 20, lmin = 10, ncluster = 2:4, seg.var =
#' c("dist","angle"))}
#' @export

segclust <- function (x, ...) {
  UseMethod("segclust", x)
}



#' Segmentation/Clustering function for data.frame
#' @rdname segclust
#' @export

segclust.data.frame <- function(x, Kmax, lmin, ncluster, type = "behavior", seg.var = NULL, diag.var = seg.var, order.var = seg.var[1], coord.names = c("x","y"),...){
  if(type == "home-range" & missing(seg.var)){
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

  segmented <- segclust_internal(x, seg.var = seg.var, diag.var = diag.var, order.var = order.var, Kmax = Kmax, ncluster = ncluster, lmin = lmin, dat=dat, type=type, ... )
  return(segmented)
}

#' Segmentation/Clustering function for move objects
#' @rdname segclust
#' @export


segclust.Move <- function(x, Kmax, lmin, ncluster, type = "behavior", seg.var = NULL, diag.var = seg.var, order.var = seg.var[1], coord.names = c("x","y"), ...){
  if(!requireNamespace("move", quietly = TRUE))
    stop("move package required for calling segclust on a Move object.")
  if(type == "home-range"){
    if(!requireNamespace("sp", quietly = TRUE))
      stop("sp package required for calling segclust (home-range) on a Move object.")
    dat <- t(sp::coordinates(x))
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

  segmented <- segclust_internal(x.df, seg.var = seg.var, diag.var = diag.var, order.var = order.var, Kmax = Kmax, ncluster = ncluster, lmin = lmin, dat=dat, type=type, ...)
  return(segmented)
}

#' Segmentation/Clustering function for ltraj objects
#' @rdname segclust
#' @export


segclust.ltraj <- function(x, Kmax, lmin, ncluster, type = "behavior", seg.var = NULL, diag.var = seg.var, order.var = seg.var[1], coord.names = c("x","y"), ...){
  if(!requireNamespace("adehabitatLT", quietly = TRUE))
    stop("adehabitatLT package required for calling segclust on a ltraj object.")

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

  segmented <- segclust_internal(x.df, seg.var = seg.var, diag.var = diag.var, order.var = order.var, Kmax = Kmax, ncluster = ncluster, lmin = lmin, dat=dat, type=type, ...)
  return(segmented)
}

#' Internal segmentation/clustering function
#' @param ... additional arguments given to \code{\link{chooseseg_lavielle}}
#' @inheritParams segclust
#' @inheritParams segmentation_internal

segclust_internal <- function(x, seg.var = NULL, diag.var = NULL, order.var = NULL, scale.variable = NULL, Kmax, ncluster = NULL, lmin = NULL, dat=NULL, type=NULL, sameSigma = F, subsample_over = 1000, subsample_by = NA, subsample = TRUE, ...){

  if(missing(Kmax)){
    Kmax = floor(dim(dat)[2]/lmin)
    message(paste("Unspecified Kmax, taking maximum possible value : Kmax = ",Kmax,". Think about reducing Kmax if running is too slow"))
  }
  if(subsample){
    x_nrow <- nrow(x)
    tmp <- subsample(x,subsample_over, subsample_by)
    x <- tmp$x
    subsample_by <- tmp$by
  }
  else{
    subsample_by <- 1
  }


  dat <- dat[,!is.na(x$subsample_ind)]
  if(missing(scale.variable)){
    message("Rescaling variables")
    scale.variable <- T
  }

  if(scale.variable) {
    dat[1,]<- scale(dat[1,])
    dat[2,]<- scale(dat[2,])
  }


  lmin <- floor(lmin/subsample_by)
  if(subsample_by > 1){
    message(paste("Adjusting lmin to subsampling. New lmin divided by",subsample_by,"and set to",lmin,"."))
  }
  if(lmin < 1){
    stop("lmin should be > 1")
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

  # DynProg segmentation ncluster=0 - Initialisation
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

  likelihood <- data.frame(nseg=1:Kmax,likelihood=-res.DynProg$J.est,ncluster=0)
  # dfBIC <- calc_BIC(likelihood,ncluster = 1:Kmax, nseg = 1:Kmax)
  # dfBIC$ncluster = 0
  segmented$outputs <- c(segmented$outputs,outputs)
  segmented$likelihood <- likelihood

  for(P in ncluster){
    # res <- segTraj::hybrid_simultanee(dat, P = P, Kmax = Kmax, lmin = lmin, sameSigma = sameSigma)
    res <- hybrid_simultanee(dat, P = P, Kmax = Kmax, lmin = lmin, sameSigma = sameSigma, ...)
    outputs <- lapply(P:Kmax,function(k){
      out <- stat_segm(x, diag.var, order.var, param = res$param[[k]], seg.type = 'segclust')
      names(out) <- c("segments","states")
      return(out)
    })
    names(outputs) <- paste(P,"class -",P:Kmax, "segments")
    likelihood = data.frame(nseg=1:Kmax,likelihood = c(res$Linc),ncluster=P)
    tmpBIC <- calc_BIC(likelihood$likelihood,ncluster = P, nseg = 1:Kmax, n = dim(x)[1])
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
