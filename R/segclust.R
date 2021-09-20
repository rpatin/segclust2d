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
#' #' @examples
#' df <-  test_data()$data
#' #' # data is a data.frame with column 'x' and 'y'
#' # Simple segmentation with automatic subsampling 
#' # if data has more than 1000 rows:
#' res <- segclust(df,
#'  Kmax = 15, lmin = 10, ncluster = 2:4, 
#'  seg.var = c("x","y"))
#'  # Plot results
#'  plot(res)
#'  segmap(res, coord.names = c("x","y"))
#'  # check penalized likelihood of 
#'  # alternative number of segment possible. 
#'  # There should be a clear break if the segmentation is good
#'  plot_BIC(res)
#' \dontrun{
#' # Advanced options:
#' # Run with automatic subsampling if df has more than 500 rows:
#' res <- segclust(df, Kmax = 30, lmin = 10, ncluster = 2:4,
#'                 seg.var = c("x","y"), subsample_over = 500)
#' # Run with subsampling by 2:
#' res <- segclust(df, Kmax = 30, lmin = 10, ncluster = 2:4,
#'                 seg.var = c("x","y"), subsample_by = 2)
#' # Disable subsampling:
#' res <- segclust(df, Kmax = 30, lmin = 10, 
#'                 ncluster = 2:4, seg.var = c("x","y"), subsample = FALSE)
#' # Disabling automatic scaling of variables for segmentation (standardazing
#' # the variables) :
#'  res <- segclust(df, Kmax = 30, lmin = 10,
#'                  seg.var = c("dist","angle"), scale.variable = FALSE)
#' }
#' @export

segclust <- function (x, ...) {
  UseMethod("segclust", x)
}



#' Segmentation/Clustering function for data.frame
#' @rdname segclust
#' @export

segclust.data.frame <-
  function(x, ...) {
  segmented <-
    segclust_internal(x, ... )
  return(segmented)
}

#' Segmentation/Clustering function for move objects
#' @rdname segclust
#' @export


segclust.Move <-
  function(x, ...){
    if(!requireNamespace("move", quietly = TRUE)){
      cli::cli_alert_danger(
        "move package not found.
        Please run install.packages('move')")
      stop("move package required for calling segclust on a Move object.")
    }
    
  segmented <- 
    segclust_internal(x,  ...)
  return(segmented)
}

#' Segmentation/Clustering function for ltraj objects
#' @rdname segclust
#' @export


segclust.ltraj <- function(x,  ...){
  
  if(!requireNamespace("adehabitatLT", quietly = TRUE))
    cli::cli_alert_danger(
      "adehabitatLT package not found.
        Please run install.packages('adehabitatLT')")
  stop("adehabitatLT package required for calling segclust() 
         on a ltraj object.")
  if (length(x) > 1){
    cli::cli_alert_warning(
      "segclust() cannot handle multi-individual ltraj objects")
    cli::cli_alert("running segclust only on the first individual")
  }
  
  segmented <-  segclust_internal(x, ...)
  return(segmented)
}

#' Internal segmentation/clustering function
#' @param ... additional arguments given to \code{\link{chooseseg_lavielle}}
#' @inheritParams segclust
#' @inheritParams segmentation_internal

segclust_internal <- 
  function(x,
           seg.var, diag.var, order.var, 
           Kmax, ncluster, lmin, 
           scale.variable, 
           sameSigma = FALSE,  ...){
    
    # Checking arguments ------------------------------------------------------
    
    
    cli::cli_h1("Checking arguments")
    # check deprecated argument 'type' and 'coord.names'
    argcheck_type_coord(...)
    # check seg.var
    tmp <- argcheck_seg.var(x, seg.var, is_segclust = TRUE)
    seg.var <- tmp$seg.var
    x <- tmp$x.df
    # format signal to be segmented for further functions
    # one signal per row
    dat <- t(x[,seg.var]) 
    
    # check lmin argumenet
    argcheck_lmin(lmin, is_segclust = TRUE)
    
    # check Kmax
    Kmax <- argcheck_Kmax(Kmax, lmin, dim(dat)[2])
    
    # check ncluster
    argcheck_ncluster(ncluster, Kmax)

    # check scale.variable
    scale.variable <- 
      argcheck_scale.variable(scale.variable, is_segclust = TRUE)
    
    # check diag.var
    diag.var <- 
      argcheck_diag.var(diag.var, seg.var)
    # check order.var
    order.var <- 
      argcheck_order.var(order.var, diag.var)
    
    # Subsampling and checks ----------------------------------
    cli::cli_h1("Preparing and checking data")
    
    cli::cli_h2("Subsampling")
    x_nrow <- nrow(x)
    x <- apply_subsampling(x, is_segclust = FALSE, ...)
    
    subsample_by <- attr(x,'subsample_by')
    if(subsample_by != 1){
      dat <- dat[,!is.na(x$subsample_ind)]
      lmin <- max(floor(lmin/subsample_by),5)
      cli::cli_alert_success(
        "Adjusting lmin to subsampling. 
        {cli::col_grey('Dividing lmin by ',
        subsample_by,', with a minimum of 5')}")
      cli::cli_alert("After subsampling, {cli::col_green('lmin = ', lmin)}. 
                    {cli::col_grey('Corresponding to lmin = ',lmin*subsample_by,
                     ' on the original time scale')}")
      
    }
    if(Kmax*lmin*subsample_by > x_nrow){
      Kmax <- min(Kmax, floor( x_nrow / ( lmin*subsample_by ) ) )
      cli::cli_alert_warning(
        "Adjusting Kmax so that lmin*Kmax < nrow(x). \\
      {cli::col_yellow('Kmax = ', Kmax)}")
    }
    
    cli::cli_h2("Scaling and final data check")
    
    if(scale.variable){
      cli::cli_alert_success(
        "Rescaling variables.
      {cli::col_grey('To desactivate, use scale.variable = FALSE')}")
      dat[1,] <- scale(dat[1,])
      dat[2,] <- scale(dat[2,])
    } else {
      cli::cli_alert_success(
        "No variable rescaling.
      {cli::col_grey('To activate, use scale.variable = TRUE')}")
      dat[1,] <- scale(dat[1,], center = TRUE, scale = FALSE)
      dat[2,] <- scale(dat[2,], center = TRUE, scale = FALSE)
    }
    
    
    
    
    # check that there are no repetitions in the series
    if(check_repetition(dat, lmin)){
      cli::cli_alert_danger(
        "Data have repetition of nearly-identical \\
        values longer than lmin. 
        {cli::col_grey('The algorithm cannot estimate variance \\
        for segment with repeated values. \\
        This is potentially caused by interpolation \\
         of missing values or rounding of values.')}
       {cli::symbol$arrow_right} Please check for repeated \\
        or very similar values of {seg.var}")
      
      stop("There are repetitions of identical values 
           in the time series larger than lmin.")
    } else {
      cli::cli_alert_success(
        "Data have no repetition of \\
        nearly-identical values larger than lmin")
    }
    cli::cli_h1("Running Segmentation/Clustering algorithm")
    cli::cli_alert_info("Running Segmentation/Clustering \\
                        with lmin = {lmin}, Kmax = {Kmax} \\
                        and ncluster = {deparse(ncluster)}")
    
    sb <- cli::cli_status(
      "{cli::symbol$arrow_right} Calculating initial \\
      segmentation without clustering")
    
    segmented <- list("data" = x,
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
    out <- stat_segm(x,
                     diag.var, order.var,
                     param = res.DynProg, 
                     nseg=k, seg.type = "segmentation")
    names(out) <- c("segments","states")
    return(out)
  })
  names(outputs) <- paste("0 class -",1:Kmax, "segments")
  
  likelihood <- data.frame(nseg=1:Kmax,likelihood=-res.DynProg$J.est,ncluster=0)
  # dfBIC <- calc_BIC(likelihood,ncluster = 1:Kmax, nseg = 1:Kmax)
  # dfBIC$ncluster = 0
  segmented$outputs <- c(segmented$outputs,outputs)
  segmented$likelihood <- likelihood
  cli::cli_alert_success("Initial segmentation with no cluster calculated.")
  
  for(P in ncluster){
    # res <- segTraj::hybrid_simultanee(dat, P = P, Kmax = Kmax, lmin = lmin,
    # sameSigma = sameSigma)
    cli::cli_h3("Segmentation/Clustering with ncluster = {P}")
    
    res <- hybrid_simultanee(dat, P = P,
                             Kmax = Kmax, lmin = lmin, 
                             sameSigma = sameSigma, ...)
    outputs <- lapply(P:Kmax,function(k){
      out <- 
        stat_segm(x,
                  diag.var, order.var, 
                  param = res$param[[k]], 
                  seg.type = 'segclust')
      names(out) <- c("segments","states")
      return(out)
    })
    names(outputs) <- paste(P,"class -",P:Kmax, "segments")
    likelihood <- data.frame(nseg=1:Kmax,likelihood = c(res$Linc),ncluster=P)
    tmpBIC <- calc_BIC(
      likelihood$likelihood,
      ncluster = P, nseg = 1:Kmax, 
      n = dim(x)[1])
    param <- list(res$param)
    names(param) <- paste(P,"class")
    segmented$likelihood <- rbind(segmented$likelihood,likelihood)
    segmented$BIC <- rbind(segmented$BIC,tmpBIC)
    segmented$param <- c(segmented$param,param)
    segmented$outputs <- c(segmented$outputs,outputs)
    segmented$Kopt.BIC[P] <- which.max(tmpBIC$BIC)
    cli::cli_alert_success(
      "Segmentation/Clustering with ncluster = {P} successfully calculated.
      BIC selected : nseg = {segmented$Kopt.BIC[P]}")
    
  }
  tmp <- dplyr::filter(segmented$BIC,ncluster > 0)
  segmented$ncluster.BIC <- tmp[which.max(tmp$BIC),"ncluster"]
  cli::cli_status_clear(id = sb)
  cli::cli_h1("Segmentation/Clustering results")
  cli::cli_alert_success(
    "Best segmentation/clustering estimated with \\
    {segmented$ncluster.BIC} cluster and \\
    {segmented$Kopt.BIC[segmented$ncluster.BIC]} segments according to BIC")
  cli::cli_text(cli::col_grey(
    '{cli::symbol$arrow_right} Number of cluster \\
    should preferentially be selected 
    according to biological knowledge. Exploring the BIC plot with plot_BIC()
    can also provide advice to select the number of cluster.'))
  cli::cli_text(cli::col_grey(
    '{cli::symbol$arrow_right} Once number of cluster is selected, \\
    the number of segment should be selected according to BIC.'))
  
  cli::cli_text(cli::col_grey(
    '{cli::symbol$arrow_right} Results of the segmentation/clustering
    may further be explored with plot() and segmap()'))
  
  return(segmented)
}
