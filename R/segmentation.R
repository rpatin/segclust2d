#' Segmentation of movement data - Generic function
#'
#' Segmentation of movement data. No clustering. Method available for
#' data.frame,  move and ltraj object. The algorithm finds for each number of
#' segment the optimal segmentation using a Dynamic Programming approach. The
#' number of segment is then chosen using Lavielle's (2005) procedure based on
#' locating rupture in the penalized likelihood.
#'
#' @param x data used for segmentation. Supported: data.frame, Move object,
#' ltraj object
#' @param type type of segmentation. Either "home-range" or "behavior". Changes
#'   default values of arguments order, scale.variable in the different
#'   functions used on the output. Default for segmentation: "home-range";
#'   default for segmentation/clustering : "behavior".
#' @param seg.var for behavioral segmentation: names of the variables used for
#'   segmentation (either one or two names).
#' @param diag.var for behavioral segmentation: names of the variables on which
#'   statistics are calculated.
#' @param order.var for behavioral segmentation: names of the variable with
#'   which states are ordered.
#' @param Kmax maximum number of segments.
#' @param lmin minimum length of segments.
#' @param ... additional parameters given to \code{\link{segmentation_internal}}
#' @inheritParams segmentation.data.frame
#' @inheritParams segmentation.Move
#' @inheritParams segmentation.ltraj
#' @return  a \code{\link{segmentation-class}} object
#'
#' @examples
#' df <-  test_data()$data
#' #' # data is a data.frame with column 'x' and 'y'
#' # Simple segmentation with automatic subsampling if data has more than 1000 rows:
#' res <- segmentation(df, Kmax = 30, lmin = 10, coord.names = c("x","y"), 
#' type = 'home-range')
#'  # Plot results
#'  plot(res)
#'  segmap(res)
#'  # check likelihood of alternative number of segment possible. There should
#'  # be a clear break if the segmentation is good
#'  plot_likelihood(res)
#' \dontrun{
#' # Advanced options:
#' # Run with automatic subsampling if df has more than 500 rows:
#' res <- segmentation(df, Kmax = 30, lmin = 10, coord.names = c("x","y"), 
#' type = 'home-range', subsample_over = 500)
#' 
#' # Run with subsampling by 2:
#' res <- segmentation(df, Kmax = 30, lmin = 10, coord.names = c("x","y"),
#'  type = 'home-range', subsample_by = 2)
#'  
#' # Disable subsampling:
#' res <- segmentation(df, Kmax = 30, lmin = 10, coord.names = c("x","y"), 
#' type = 'home-range', subsample = FALSE)
#' 
#' # Run on other kind of variables : 
#'  res <- segmentation(df, Kmax = 30, lmin = 10, seg.var = c("dist","angle"), 
#'  type = 'behavior')
#'  
#' # Automatic scaling of variables for segmentation 
#' (set a mean of 0 and a standard deviation of 1 for both variables)
#' 
#'  res <- segmentation(df, Kmax = 30, lmin = 10, seg.var = c("dist","angle"),
#'   type = 'behavior', scale.variable = TRUE)
#'  
#' }
#' @export

segmentation <- function(x, ...) {
  UseMethod("segmentation", x)
}


#' Segmentation function for data.frame
#' @param coord.names for home-range segmentation and data.frame, names of coordinates. Default x and y. Not needed for move and ltraj objects
#' @rdname segmentation
#' @export

segmentation.data.frame <- function(x, Kmax, lmin, type = "home-range", seg.var, diag.var = seg.var, order.var = seg.var[1], coord.names = c("x","y"), ...){

  if(type == "home-range"){
    if(missing(seg.var)){
      seg.var <- coord.names
    }
    dat <- t(x[,seg.var])
    seg.var = seg.var
    diag.var = seg.var
    order.var = seg.var[1]
  } else if ( type == "behavior" ){
    if(is.null(seg.var)) stop("Please provide seg.var for a behavioral segmentation")
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
  message("Segmentation on ",paste(seg.var,collapse = " and "))
  segmented <- segmentation_internal(x, seg.var = seg.var, diag.var = diag.var, order.var = order.var,  Kmax = Kmax, lmin = lmin, dat=dat, type=type, ...)
  return(segmented)
}


#' Segmentation function for move objects
#' @rdname segmentation
#' @export

segmentation.Move <- function(x, Kmax, lmin, type = "home-range", seg.var = NULL, diag.var = seg.var, order.var = seg.var[1], coord.names = c("coords.x1","coords.x2"), ...){

  if(!requireNamespace("move", quietly = TRUE))
    stop("move package required for calling segmentation on a Move object.")

  if(type == "home-range"){
    if(!requireNamespace("sp", quietly = TRUE))
      stop("sp package required for calling segmentation (home-range) on a Move object.")
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

  segmented <- segmentation_internal(x.df, seg.var = seg.var, diag.var = diag.var, order.var = order.var, Kmax = Kmax, lmin = lmin, dat=dat, type=type, ...)
  return(segmented)
}

#' Segmentation function for ltraj objects
#' @rdname segmentation
#' @export

segmentation.ltraj <- function(x, Kmax, lmin, type = "home-range", seg.var = NULL, diag.var = seg.var, order.var = seg.var[1], coord.names = c("x","y"), ...){
  if(!requireNamespace("adehabitatLT", quietly = TRUE))
    stop("adehabitatLT package required for calling segmentation on a ltraj object.")

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

  segmented <- segmentation_internal(x.df, seg.var = seg.var, diag.var = diag.var, order.var = order.var, Kmax = Kmax, lmin = lmin, dat=dat,type=type, ...)
  return(segmented)
}


#' Internal segmentation function
#' @param x data.frame with observations
#' @param dat bivariate data (one signal per row)
#' @param sameSigma does segments have same variance ?
#' @param scale.variable should variables be standardized ? (reduced and centered)
#' @param subsample_over over which size should subsampling begin (depending on
#'   speed and memory limitations)
#' @param subsample if FALSE disable subsample
#' @param subsample_by override subsample_over to adjust manually subsampling
#' @param ... additionnal parameters given to \link{chooseseg_lavielle}
#' @inheritParams segmentation
#'
#' @inheritParams chooseseg_lavielle
#' @export

segmentation_internal <- function(x, seg.var = NULL, diag.var = NULL, order.var = NULL, scale.variable = NULL, Kmax, lmin = NULL, dat=NULL, type=NULL, sameSigma = F, subsample_over = 10000, subsample = TRUE, subsample_by = NA, ...){



  x_nrow <- nrow(x)
  if(missing(Kmax)){
    Kmax = floor(x_nrow/lmin)
    message(paste("Unspecified Kmax, taking maximum possible value: Kmax = ",Kmax,". Think about reducing Kmax if running is too slow"))
  }
  if(subsample){
    x_nrow <- nrow(x)
    tmp <- subsample(x,subsample_over, subsample_by)
    x <- tmp$x
    subsample_by <- tmp$by
    dat <- dat[,!is.na(x$subsample_ind)]
  } else {
    subsample_by <- 1
    x$subsample_ind <- 1:nrow(x)
  }
  
  if(missing(scale.variable) & type == 'behavior'){
    message("Rescaling variables")
    scale.variable <- TRUE
  } else {
    scale.variable <- FALSE
  }

  if(scale.variable) {
    dat[1,]<- scale(dat[1,])
    dat[2,]<- scale(dat[2,])
  }
  lmin <- max(floor(lmin/subsample_by),2)
  Kmax <- min(Kmax, floor( x_nrow / ( lmin*subsample_by ) ) )
  
  # check that there are no repetitions in the series
  if(check_repetition(dat, lmin)){
    stop("There are repetitions of identical values in the time series larger than lmin, cannot estimate variance for such segment. This is potentially caused by interpolation of missing values or rounding of values.")
  }
  dat[1,] <- dat[1,]-mean(dat[1,])
  dat[2,] <- dat[2,]-mean(dat[2,])
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
  output_lavielle <- chooseseg_lavielle(res.DynProg$J.est, ...)
  # stationarity = test_stationarity(dat,outputs,Kmax)
  # seg_var = test_var(dat,outputs,Kmax)
  # seg_mean = test_mean(dat,outputs,Kmax)

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
