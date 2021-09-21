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
#' @param ... additional parameters given to \code{\link{segmentation_internal}}
#' @return  a \code{\link{segmentation-class}} object
#'
#' @examples
#' df <-  test_data()$data
#' #' # data is a data.frame with column 'x' and 'y'
#' # Simple segmentation with automatic subsampling
#' # if data has more than 1000 rows:
#' res <- segmentation(df, Kmax = 30, lmin = 10, seg.var = c("x","y"))
#'  # Plot results
#'  plot(res)
#'  segmap(res)
#'  # check likelihood of alternative number of segment possible. There should
#'  # be a clear break if the segmentation is good
#'  plot_likelihood(res)
#' \dontrun{
#' # Advanced options:
#' # Run with automatic subsampling if df has more than 500 rows:
#' res <- segmentation(df, Kmax = 30, lmin = 10,
#'  seg.var = c("x","y"),  subsample_over = 500)
#' 
#' # Run with subsampling by 2:
#' res <- segmentation(df, Kmax = 30, lmin = 10, 
#' seg.var = c("x","y"), subsample_by = 2)
#'  
#' # Disable subsampling:
#' res <- segmentation(df, Kmax = 30, lmin = 10,
#'  seg.var = c("x","y"), subsample = FALSE)
#' 
#' # Run on other kind of variables : 
#'  res <- segmentation(df, Kmax = 30, lmin = 10, seg.var = c("dist","angle"))
#'  
#' # Automatic scaling of variables for segmentation 
#' (set a mean of 0 and a standard deviation of 1 for both variables)
#' 
#'  res <- segmentation(df, Kmax = 30, lmin = 10, 
#'  seg.var = c("dist","angle"), scale.variable = TRUE)
#'  
#' }
#' @export

segmentation <- function(x, ...) {
  UseMethod("segmentation", x)
}


#' Segmentation function for data.frame
#' @rdname segmentation
#' @export

segmentation.data.frame <-
  function(x, ...){
    
  segmented <-
      segmentation_internal(x, ...)
    return(segmented)
  }


#' Segmentation function for move objects
#' @rdname segmentation
#' @export

segmentation.Move <- 
  function(x,  ...){
    
    if(!requireNamespace("move", quietly = TRUE)){
      cli::cli_alert_danger(
        "move package not found.
        Please run install.packages('move')")
      stop("move package required for calling segclust on a Move object.")
    }
    
    segmented <-
      segmentation_internal(x, ...)
    return(segmented)
  }

#' Segmentation function for ltraj objects
#' @rdname segmentation
#' @export

segmentation.ltraj <- 
  function(x, ...){
    
    if(!requireNamespace("adehabitatLT", quietly = TRUE))
      cli::cli_alert_danger(
        "adehabitatLT package not found.
        Please run install.packages('adehabitatLT')")
    stop("adehabitatLT package required for calling segmentation 
         on a ltraj object.")
    if (length(x) > 1){
      cli::cli_alert_warning(
        "segmentation() cannot handle multi-individual ltraj objects")
      cli::cli_alert("running segmentation only on the first individual")
    }

      segmented <-
      segmentation_internal(
        x, ...)
    return(segmented)
  }


#' Internal segmentation function
#' @param x data.frame with observations
#' @param sameSigma does segments have same variance ?
#' @param seg.var names of the variables used for
#'   segmentation (either one or two names).
#' @param Kmax maximum number of segments.
#' @param lmin minimum length of segments.
#' @param diag.var names of the variables on which
#'   statistics are calculated.
#' @param order.var names of the variable with which states are ordered.
#' @param scale.variable TRUE or FALSE for automatic scaling of variables
#'  (reduction and  centering)
#' @param ... additional parameters given to \link{chooseseg_lavielle}
#'
#' @rdname segmentation
#' 
#' @export

segmentation_internal <-
  function(x,
           seg.var, diag.var, order.var,
           lmin, Kmax,
           scale.variable, 
           sameSigma = FALSE,
           ...){
   

# Checking arguments ------------------------------------------------------

    
    cli::cli_h1("Checking arguments")
    # check deprecated argument 'type' and 'coord.names'
    argcheck_type_coord(...)
    # check seg.var
    tmp <- argcheck_seg.var(x, seg.var, is_segclust = FALSE)
    seg.var <- tmp$seg.var
    x <- tmp$x.df
    # format signal to be segmented for further functions
    # one signal per row
    dat <- t(x[,seg.var]) 
    
    # check lmin argumenet
    argcheck_lmin(lmin, is_segclust = FALSE)
    
    # check Kmax
    Kmax <- argcheck_Kmax(Kmax, lmin, dim(dat)[2])
    
    # check scale.variable
    scale.variable <- 
      argcheck_scale.variable(scale.variable, is_segclust = FALSE)
    
    # check diag.var
    diag.var <- 
      argcheck_diag.var(diag.var, seg.var)
    # check order.var
    order.var <- 
      argcheck_order.var(order.var, diag.var)

# Subsampling and checks -----------------------------
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
      if(Kmax >= 2){
        cli::cli_alert_success(
          "Adjusting Kmax so that lmin*Kmax < nrow(x). \\
      {cli::col_green('Kmax = ', Kmax)}")
      } else {
        cli::cli_alert_danger(
          "lmin*Kmax > nrow(x) and Kmax cannot be adjusted. \\
          Please provide lower values for lmin")
        stop("lmin*Kmax > nrow(x)")
      }
      
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
      
      stop("There are repetitions of identical
           values in the time series larger than lmin.")
    } else {
      cli::cli_alert_success(
        "Data have no repetition of \\
        nearly-identical values larger than lmin")
    }
    
    cli::cli_h1("Running segmentation algorithm")
    cli::cli_alert_info("Running segmentation \\
                        with lmin = {lmin} and Kmax = {Kmax}")
    sb <- cli::cli_status("{cli::symbol$arrow_right} Calculating cost matrix")
    CostLoc <- Gmean_simultanee(dat, lmin = lmin,sameVar = sameSigma)
    cli::cli_alert_success("Cost matrix calculated")
    cli::cli_status_update(id = sb,
                      "{cli::symbol$arrow_right} Dynamic Programming")
    res.DynProg <- wrap_dynprog_cpp(CostLoc, Kmax)
    cli::cli_alert_success(
      "Optimal segmentation calculated \\
      for all number of segments <= {Kmax}")
    # CostLoc <- segTraj::Gmean_simultanee(dat, lmin = lmin)
    # res.DynProg <- segTraj::DynProg(CostLoc, Kmax)
    cli::cli_status_update(
      id = sb,
      "{cli::symbol$arrow_right} Calculating segment statistics")
    outputs <- lapply(1:Kmax,function(k){
      # print(k)
      out <- stat_segm(x, 
                       diag.var = diag.var, order.var = order.var,
                       param = res.DynProg, 
                       seg.type = 'segmentation', nseg=k)
      names(out) <- c("segments","states")
      return(out)
    })
    
    names(outputs) <- paste(1:Kmax, "segments")
    output_lavielle <- chooseseg_lavielle(res.DynProg$J.est, ...)
    # stationarity = test_stationarity(dat,outputs,Kmax)
    # seg_var = test_var(dat,outputs,Kmax)
    # seg_mean = test_mean(dat,outputs,Kmax)
    cli::cli_status_clear(id = sb)
    cli::cli_alert_success("Best segmentation estimated with \\
                      {output_lavielle$Kopt} segments, \\
                      according to Lavielle's criterium")
    cli::cli_text(cli::col_grey(
    'Other number of segment may be selected 
      by looking for likelihood breaks with plot_likelihood()'))
    cli::cli_text(cli::col_grey(
    'Results of the segmentation may be explored with plot() and segmap()'))
    segmented <- 
      list("data" = x,
           "seg.type" = "segmentation",
           "outputs" = outputs,
           "likelihood" = 
             data.frame(nseg=1:Kmax,
                        likelihood=-res.DynProg$J.est),
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
