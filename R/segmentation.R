#' Segmentation of movement data - Generic function
#'
#' Segmentation of movement data. No clustering. Method available for data.frame,  move and ltraj object
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

segmentation <- function (x, ...) {
  UseMethod("segmentation", x)
}


#' Segmentation function for data.frame
#' @param coord.names for home-range segmentation and data.frame, names of coordinates. Default x and y. Not needed for move and ltraj objects
#' @rdname segmentation
#' @export

segmentation.data.frame <- function(x, Kmax = 10, lmin = Kmax/2, type = "home-range", scale.variable = F, seg.var = NULL, diag.var = seg.var, order.var = seg.var[1], coord.names = c("x","y"),S=0.75,sameSigma = F){

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

  segmented <- segmentation_internal(x, seg.var = seg.var, diag.var = diag.var, order.var = order.var, scale.variable = scale.variable, Kmax = Kmax, lmin = lmin, dat=dat,data.type = "data.frame",S=S,type=type,sameSigma = sameSigma)
  return(segmented)
}


#' Segmentation function for move objects
#' @rdname segmentation
#' @export

segmentation.Move <- function(x, Kmax = 10, lmin = Kmax/2, type = "home-range", scale.variable = F, seg.var = NULL, diag.var = seg.var, order.var = seg.var[1], coord.names = c("coords.x1","coords.x2"),S=0.75,sameSigma = F){

  if(type == "home-range"){
    dat <- t(coordinates(datamove))
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

  segmented <- segmentation_internal(x.df, seg.var = seg.var, diag.var = diag.var, order.var = order.var, scale.variable = scale.variable, Kmax = Kmax, lmin = lmin, dat=dat,data.type = "move",S=S,type=type,sameSigma = sameSigma)
  return(segmented)
}

#' Segmentation function for ltraj objects
#' @rdname segmentation
#' @export

segmentation.ltraj <- function(x, Kmax = 10, lmin = Kmax/2, type = "home-range", scale.variable = F, seg.var = NULL, diag.var = seg.var, order.var = seg.var[1], coord.names = c("x","y"),S=0.75,sameSigma = F){

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

  segmented <- segmentation_internal(x.df, seg.var = seg.var, diag.var = diag.var, order.var = order.var, scale.variable = scale.variable, Kmax = Kmax, lmin = lmin, dat=dat,data.type = "move",S=S,type=type,sameSigma = sameSigma)
  return(segmented)
}


#' Internal segmentation function

segmentation_internal <- function(x, seg.var = NULL, diag.var = NULL, order.var = NULL, scale.variable = scale.variable, Kmax = NULL, lmin = NULL, dat=NULL, data.type = NULL,S=NULL,type=NULL,sameSigma = NULL){

  if(scale.variable) {
    dat[1,]<- scale(dat[1,])
    dat[2,]<- scale(dat[2,])
  }

  CostLoc <- Gmean_simultanee(dat, lmin = lmin,sameVar = sameSigma)
  res.DynProg <- DynProg(CostLoc, Kmax)

  # CostLoc <- segTraj::Gmean_simultanee(dat, lmin = lmin)
  # res.DynProg <- segTraj::DynProg(CostLoc, Kmax)

  outputs <- lapply(1:Kmax,function(k){
    # print(k)
    out <- stat_segm(x, diag.var = diag.var, order.var = order.var, param = res.DynProg, seg.type = 'segmentation', nseg=k)
    names(out) <- c("segments","states")
    return(out)
  })
  names(outputs) <- paste(1:Kmax, "segments")
  output_lavielle <- chooseseg_lavielle(res.DynProg$J.est,S=S)
  segmented <- list("data" = x,
                    "type" = type,
                    "seg.type" = "segmentation",
                    "outputs" = outputs,
                    "likelihood" = data.frame(nseg=1:Kmax,likelihood=-res.DynProg$J.est),
                    "Segmented variables" = seg.var,
                    "Diagnostic variables" = diag.var,
                    "Order variable" = order.var,
                    "param"= list("lmin"=lmin,
                                  "Kmax"=Kmax),
                    "Kopt.lavielle" = output_lavielle["Kopt"][[1]],
                    "df.lavielle" = output_lavielle["lavielle"][[1]])
  class(segmented) <- "segmentation"
  return(segmented)
}
