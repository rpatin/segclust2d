
#' Internal function for subsampling
#'
#' if subsample = FALSE do nothing.
#' 
#' else if subsample_by is missing, subsample only 
#' if nrow(x) > subsample_over,
#' then it subsample with the minimum needed to get a
#' data.frame smaller than subsample_over
#' 
#' if subsample_by is provided, use it to subsample.
#' 
#' @param x data.frame to be subsampled
#' @param is_segclust TRUE or FALSE whether the function was called from
#' `segclust()` or `segmentation()`
#' @param subsample if FALSE disable subsampling
#' @param subsample_over maximum number of row accepted
#' @param subsample_by subsampling parameters
#' @return a data.frame

apply_subsampling <-
  function(x, is_segclust, subsample, subsample_over, subsample_by){
    
    if(!methods::hasArg(subsample)){
      cli::cli_alert_warning(
        "Subsampling automatically activated. \\
        To disable it, provide {cli::col_yellow('subsample = FALSE')}")
      subsample <- TRUE
    } else {
      if(subsample){
        cli::cli_alert_success("Subsampling was activated with \\
                           {cli::col_green('subsample = TRUE')}")
      }
    }
    x_nrow <- nrow(x)
    x_ind <-  seq_len(x_nrow)
    x$x_ind <- x_ind
    
    if(subsample){
      if(!methods::hasArg(subsample_by)){
        if(!methods::hasArg(subsample_over)){
          if(is_segclust){
            subsample_over <- 1000
          } else {
            subsample_over <- 10000
          }
          cli::cli_alert_info(
            "Argument {cli::col_cyan('subsample_over')} was not provided")
          cli::cli_text(cli::col_grey(
            "Taking default value for 
          {ifelse(is_segclust,'segclust()','segmentation()')}"))
          cli::cli_text(cli::col_grey(
            "Setting subsample_over = {subsample_over}"))
        } else {
          cli::cli_alert_success(
            "Using {cli::col_green('subsample_over = ', subsample_over )}") 
        }
        if( x_nrow > subsample_over){
          subsample_by <- ceiling(x_nrow/subsample_over)
          cli::cli_alert_success(
            "nrow(x) > subsample_over, \\
          {cli::col_green('subsampling by ',subsample_by)}")
          keep <- (x_ind %% subsample_by) == 1
          subsample_ind <- ifelse(keep,(x_ind %/% subsample_by)+1,NA)
          x$subsample_ind <- subsample_ind
        } else {
          subsample_by <- 1
          x$subsample_ind <- x_ind
          cli::cli_alert_success(
            "nrow(x) < subsample_over, no subsample needed")
        }
      } else {
        cli::cli_alert_success(
          "Using {cli::col_green('subsample_by = ', subsample_by )}") 
        if(subsample_by != 1){
          keep <- (x_ind %% subsample_by) == 1
          subsample_ind <- ifelse(keep,(x_ind %/% subsample_by)+1,NA)
          x$subsample_ind <- subsample_ind
          cli::cli_alert_success(
            "{cli::col_green('subsampling by ',subsample_by)}")
        } else {
          cli::cli_alert_success(
            "subsample_by set to {subsample_by}. No subsampling required.")
        }
      } # end if(missing(subsample_by))
    } else {
      x$subsample_ind <- x_ind
      subsample_by <- 1
      cli::cli_alert_success("Subsampling was deactivated with \\
                           {cli::col_green('subsample = FALSE')}")
    } # end if(subsample)
    attr(x,'subsample_by') <- subsample_by
    return(x)
  }

#' Internal function for subsampling
#'
#' merge subsampled data.frame df with fulldata to add segmentation information
#' on the full data.frame
#' @param df subsampled data.frame with additional information on segmentation
#' @param fulldata full data.frame
#' @param colname column name
subsample_rename <- function(df, fulldata, colname){
  translate_ind  <- with(fulldata,data.frame(x_ind,subsample_ind))
  translate_ind <- with(translate_ind, translate_ind[!is.na(subsample_ind),])
  var_join <- c("subsample_ind" = paste(colname) )
  df <- dplyr::right_join(translate_ind,df,var_join)
  df[,colname] <- df$x_ind
  df$x_ind <- NULL
  df
}
