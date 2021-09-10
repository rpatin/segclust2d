#' Check for deprecated 'type' and 'coord.names' argument
#'
#' Check whether argument 'type' and 'coord.names' 
#' were provided and communicate adequately if need be.
#'
#' @param ... additional parameters transmitted from \code{\link{segmentation}}
#' or \code{\link{segclust}}
#' @return  a NULL object


argcheck_type_coord <-
  function(...){
    if (hasArg(type)) {
      type <- list(...)$type
      cli::cli_alert_danger("Argument {cli::col_red('type')} \\
                     is deprecated and should not be used")
      if(type == "home-range"){
        if(!hasArg(coord.names)){
          coord.names <- c("x","y")
        } else {
          coord.names <- list(...)$coord.names
          cli::cli_alert_danger("Argument {cli::col_red('coord.names')} \\
                     is deprecated and should not be used")
        }
        cli::cli_alert("Please use instead \\
                    {.field seg.var = {deparse(coord.names)}} and \\
                    {.field scale.variable = FALSE}")
      } else if(type == "behaviour"){
        cli::cli_alert("Please use instead \\
                    {.field seg.var} argument and \\
                    {.field scale.variable = TRUE}")
        
      }
      stop("argument 'type' is deprecated and should not be used")
    } else {
      # check coord.names
      if(hasArg(coord.names)){
        coord.names <- list(...)$coord.names
        cli::cli_alert_danger("Argument {cli::col_red('coord.names')} \\
                     is deprecated and should not be used")
        cli::cli_alert("Please use instead \\
                    {.field seg.var = {deparse(coord.names)}}")
        stop("argument 'coord.names' is deprecated and should not be used")
      }
    }
  }


#' Check for argument 'seg.var'
#'
#' Check whether argument 'seg.var' was adequately provided.
#' If provided, additionnaly check for its length (1 or 2) and
#' for the existence of corresponding column names in x
#' If unprovided, use default value (segmentation only) and tests
#' if column names exists.
#' 
#' @param x data used for segmentation. Supported: data.frame, Move object,
#' ltraj object
#' @param seg.var for behavioral segmentation: names of the variables used for
#'   segmentation (either one or two names).
#' @param is_segclust TRUE if function is called from \code{\link{segclust}} ; 
#' FALSE otherwise, if function is called from \code{\link{segmentation}}.
#' @param ... additional parameters transmitted from \code{\link{segmentation}}
#' or \code{\link{segclust}}
#' @return  a list with a data.frame and a vector with two character strings

argcheck_seg.var <- function(x, seg.var, is_segclust){
  
  if(class(x) == "data.frame"){
    x.df <- x
  } else if (class(x) == "Move") {
    x.df <- x@data
  } else if (class(x) == "ltraj") {
    x.infolocs <- adehabitatLT::infolocs(ltrajshift)[[1]]
    x.df <- cbind(x.infolocs,
                  dplyr::select(x[[1]],
                                - dplyr::any_of(colnames(x.infolocs)))
    )
    
  } else if (class(x) == "sftraj") {
    
  }
  
  # is argument missing ?
  if(missing(seg.var)){
    if(is_segclust){
      cli::cli_alert_danger(
        "segclust requires argument {cli::col_red('seg.var')}.")
      cli::cli_alert("Please provide up to two variables names \\
        in {cli::col_red('seg.var')} ")
      stop("seg.var argument required")
    } else {
      cli::cli_alert_warning(
        "Argument {cli::col_yellow('seg.var')} missing")
      
      if(class(x) == "data.frame"){
        cli::cli_text(cli::col_grey(
          "taking default value \\
        {cli::col_yellow('seg.var = c(\"x\",\"y\")')}")) 
        seg.var <- c("x","y")
      } else if (class(x) == "Move") {
        
        cli::cli_text(cli::col_grey(
          "taking coordinates as default value for a move object"))  
        if(!requireNamespace("sp", quietly = TRUE)){
          cli::cli_alert_danger(
            "sp package not found.  Please run install.packages('sp')")
          stop("sp package required for calling
             segmentation on coordinates of a Move object.")
        } else {
          x_coord <- sp::coordinates(x)
          seg.var <- colnames(x_coord)
          x.df[,seg.var[1]] <- x_coord[,1]
          x.df[,seg.var[2]] <- x_coord[,2]
        }
        
      } else if (class(x) == "ltraj") {
        cli::cli_text(cli::col_grey(
          "taking coordinates as default value for a ltraj object"))
        seg.var <- c("x","y")
      } else if (class(x) == "sftraj") {
        
      }
      
    }
  } 
  # is argument length valid ?
  if( length(seg.var) == 1 ){
    cli::cli_alert_warning(
      "Argument {cli::col_yellow('seg.var')} \\
        has only one variable  {cli::col_yellow({deparse(seg.var)})}. \\
        Changing {cli::col_yellow('seg.var')} \\
        to {cli::col_yellow({deparse(rep(seg.var,2))})}")  
    seg.var <- rep(seg.var,2)
  } else if ( length(seg.var) > 2 ) {
    cli::cli_alert_danger(
      "Segmentation is not supported for more than two variables.")
    cli::cli_alert("Please provide only two variables in {cli::col_red('seg.var')} ")
    stop("seg.var should have length 1 or 2")
  } 
  # do variable exists?
  check_seg.var <- seg.var %in% colnames(x.df) 
  if(!all(check_seg.var)){
    cli::cli_alert_danger(
      "column names {cli::col_red(seg.var[!check_seg.var])} \\
      not found in the data provided.")
    stop("column names provided in seg.var not found")
  }
  cli::cli_alert_success("Segmentation with \\
                        {.field seg.var = {deparse(seg.var)}}") 
  return(list("x.df" = x.df, 
              "seg.var" = seg.var))
}


#' Check for argument 'lmin'
#'
#' Check whether argument 'lmin' was provided and is adequate
#' before subsampling
#' 
#' @param lmin minimum length of segments.
#' @param is_segclust TRUE if function is called from \code{\link{segclust}} ; 
#' FALSE otherwise, if function is called from \code{\link{segmentation}}.
#' @return  a NULL object

argcheck_lmin <- function(lmin, is_segclust){
  # is argument missing ?
  if(missing(lmin)){
    cli::cli_alert_danger(
      "{ifelse(is_segclust,'segclust()','segmentation()')} \\
        requires argument lmin.")
    cli::cli_alert("Please provide lmin, the minimum length of segments")
    stop("lmin argument required")
  } 
  if(lmin < 5){
    cli::cli_alert_danger(
      "User provided {cli::col_red('lmin = ', lmin )}
        {cli::col_red('lmin')} should be {cli::col_red('>= 5')} \\
        to avoid variance estimation instability.")
    stop("lmin < 5 not allowed")
  }
  cli::cli_alert_success("Using {cli::col_green('lmin = ', lmin )}") 
  
  return(invisible(NULL))
}



#' Check for argument 'Kmax'
#'
#' Check whether argument 'Kmax' was provided and is adequate
#' before subsampling. Propose adequate value if Kmax is not
#' provided.
#' 
#' @param lmin minimum length of segments.
#' @param Kmax maximum number of segments.
#' @param datalength length of data provided
#' @return an integer

argcheck_Kmax <- function(Kmax, lmin, datalength){
  # is argument missing ?
  if(missing(Kmax)){
    Kmax <- floor(0.75 * datalength/lmin)
    cli::cli_alert_warning(
      "Argument {cli::col_yellow('Kmax')} missing")
    cli::cli_text(cli::col_grey("Taking default value \\
                  Kmax = floor(0.75 * n/lmin) = {Kmax}."))
    cli::cli_alert(cli::col_grey(
      "Think about reducing Kmax if running time is high"))
  } else {
    cli::cli_alert_success("Using {cli::col_green('Kmax = ', Kmax )}") 
  }
  
  return(Kmax)
}

#' Check for argument 'ncluster'
#'
#' Check whether argument 'ncluster' was provided and is adequate
#' 
#' @param ncluster number of cluster into which segments should be grouped. Can
#'   be a vector if one want to test several number of clusters.
#' @param Kmax maximum number of segments.
#' @return  a NULL object

argcheck_ncluster <- function(ncluster, Kmax){
  # is argument missing ?
  if(missing(ncluster)){
    cli::cli_alert_danger(
      "segclust() requires argument {cli::col_red('ncluster')}")
    cli::cli_alert(
      "Please provide ncluster, a vector (of length >= 1) \\
      with the different number of cluster to be tested")
    stop("ncluster argument required")
  } 
  if(max(ncluster) > Kmax){
    cli::cli_alert_danger(
      "User tried {cli::col_red('ncluster = ', max(ncluster) )} with \\
      {cli::col_red('Kmax = ', Kmax )}
        {cli::col_red('ncluster')} should be {cli::col_red('<= Kmax')}.")
    stop("ncluster > Kmax not allowed")
  }
  cli::cli_alert_success("Using {cli::col_green('ncluster = ', 
                         deparse(ncluster) )}") 
  
  return(invisible(NULL))
}

#' Check for argument 'scale.variable'
#'
#' Check whether argument 'scale.variable' was provided.
#' If not, propose default value.
#' 
#' @param scale.variable minimum length of segments.
#' @param is_segclust TRUE if function is called from \code{\link{segclust}} ; 
#' FALSE otherwise, if function is called from \code{\link{segmentation}}.
#' @return a boolean

argcheck_scale.variable <- 
  function(scale.variable, is_segclust){
    # is argument missing ?
    if(missing(scale.variable)){
      cli::cli_alert_warning(
        "Argument {cli::col_yellow('scale.variable')} missing")
      if(is_segclust){
        cli::cli_text(cli::col_grey(
          "Taking default value scale.variable = TRUE for segclust()."))
        scale.variable <- TRUE
      } else {
        cli::cli_text(cli::col_grey(
          "Taking default value scale.variable = FALSE for segmentation()."))
        scale.variable <- FALSE
      }
    } else {
      cli::cli_alert_success("Using {cli::col_green('scale.variable = ', scale.variable)}") 
    }
    return(scale.variable)
  }

#' Check for argument 'diag.var'
#'
#' Check whether argument 'diag.var' was provided.
#' If not, propose default value.
#' 
#' @param diag.var names of the variables on which
#'   statistics are calculated.
#' @param seg.var for behavioral segmentation: names of the variables used for
#'   segmentation (either one or two names).

#' @return a vector of character string

argcheck_diag.var <- 
  function(diag.var, seg.var){
    # is argument missing ?
    if(missing(diag.var)){
      cli::cli_alert_info(
        "Argument {cli::col_cyan('diag.var')} was not provided")
      cli::cli_text(cli::col_grey(
        "Taking default seg.var as diagnostic variables diag.var."))
      cli::cli_text(cli::col_grey(
        "Setting diag.var = {deparse(seg.var)}"))
      diag.var <- seg.var
    } else {
      cli::cli_alert_success("Using diagnostic variables {cli::col_green('diag.var = ', {deparse(diag.var)})})}") 
    }
    return(diag.var)
  }

#' Check for argument 'order.var'
#'
#' Check whether argument 'order.var' was provided.
#' If not, propose default value.
#' 
#' @param order.var names of the variable with which states are ordered.
#' @param diag.var names of the variables on which
#'   statistics are calculated.
#' @return a vector of character string

argcheck_order.var <- 
  function(order.var, diag.var){
    # is argument missing ?
    if(missing(order.var)){
      cli::cli_alert_info(
        "Argument {cli::col_cyan('order.var')} was not provided")
      cli::cli_text(cli::col_grey(
        "Taking default diag.var[1] as ordering variable order.var."))
      cli::cli_text(cli::col_grey(
        "Setting order.var = {deparse(diag.var[1])}"))
      order.var <- diag.var[1]
    } else {
      cli::cli_alert_success("Using ordering variable {cli::col_green('diag.var = ', {deparse(order.var)})}") 
    }
    return(order.var)
  }


#' Check for argument 'order'
#'
#' Check whether argument 'order' was provided for plot.segmentation
#' and segmap. If not, propose default value.
#' 
#' @param order TRUE or FALSE depending on whether cluster be ordered
#' @param seg.type types of the segmentation.
#' @return a boolean

argcheck_ordering <- 
  function(order, seg.type, order.var){
    if (missing(order)) {
      cli::cli_alert_info(
        "Argument {cli::col_cyan('order')} missing.")
      order <- ifelse(seg.type == "segmentation", FALSE, TRUE)
      if(order){
        cli::cli_text(
          cli::col_grey(
            "Ordering cluster with variable 
          {order.var} for segmentation/clustering.
          To disable, use order = FALSE"))
      } else {
        cli::cli_text(
          cli::col_grey(
            "Ordering of segment disabled for segmentation.
          To allow ordering of segment, use order = TRUE"))
      }
    } else {
      if(order){
        cli::cli_alert_success(
          "Ordering cluster with variable 
          {order.var}. To disable, use order = FALSE")
      } else {
        cli::cli_alert_success(
          "Ordering of segment disabled. To allow ordering of segment, use order = TRUE")
      }
    }
    return(order)
  }


#' Check for argument 'ncluster' and 'nseg'
#'
#' Check whether argument 'ncluster' and 'nseg' were provided.
#'  If not, propose default value based on BIC.
#' 
#' @param ncluster number of cluster
#' @param nseg number of segment
#' @param ncluster.BIC optimal number of cluster selected by BIC
#' @param Kopt.BIC optimal number of segment selected by BIC for 
#' each number of cluster
#' @return a list with two integers

argcheck_segclust <- 
  function(ncluster, 
           nseg, 
           ncluster.BIC,
           Kopt.BIC){
    if (missing(ncluster)) {
      ncluster <- ncluster.BIC
      cli::cli_alert_warning(
        "Argument ncluster was not provided. Selecting values with BIC")
        
      if(!missing(nseg)){
        cli::cli_alert_danger(
          "Argument nseg can only be provided if ncluster is provided as well. \\
          Overwriting nseg with BIC-selected value")
      }
      nseg <- Kopt.BIC[ncluster]
      cli::cli_alert_info(
        "BIC-selected number of class : ncluster = {ncluster}
         BIC-selected number of segment : nseg = {nseg}")
      
    } else if (missing(nseg)) {
      cli::cli_alert_success(
        "Argument ncluster was provided but not nseg. \\
        Selecting nseg value with BIC")
      
      nseg <- Kopt.BIC[ncluster]
      cli::cli_alert_success(
        "User-selected number of class : ncluster = {ncluster}
         BIC-selected number of segment : nseg = {nseg}")
      
    } else {
      cli::cli_alert_info(
        "User-selected number of class : ncluster = {ncluster}
         User-selected number of segment : nseg = {nseg}")
    }
    return(list("nseg" = nseg,
                "ncluster" = ncluster))
  }


#' Check for argument 'nseg'
#'
#' Check whether argument  'nseg' was provided.
#'  If not, propose default value based on Lavielle's criterium
#' 
#' @param nseg number of segment
#' @param Kopt.lavielle optimal number of segment selected with 
#' Lavielle's criterium
#' @return an integer

argcheck_segmentation <- 
  function(nseg, Kopt.lavielle){
  if (missing(nseg)) {
      cli::cli_alert_warning(
        "Argument nseg was not provided. \\
        Selecting nseg value with Lavielle's criterium")
      nseg <- Kopt.lavielle
      cli::cli_alert_info(
        "Lavielle-selected number of segment : nseg = {nseg}")
    } else {
      cli::cli_alert_success(
        "User-selected number of segment : nseg = {nseg}")
    }
    return(nseg)
  }