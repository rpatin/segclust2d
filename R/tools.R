
#' Internal function for subsampling
subsample <- function(x,subsample_over){
  x_nrow <- nrow(x)
  x_ind <-  1:x_nrow
  x$x_ind <- x_ind

  if( x_nrow > subsample_over){
    subsample_by <- ceiling(x_nrow/subsample_over)
    warning(paste("Number of data > ",subsample_over,", subsampling by ",subsample_by,". Adjusting also lmin.",sep=""))
    keep <- (x_ind %% subsample_by) == 1
    subsample_ind <- ifelse(keep,(x_ind %/% subsample_by)+1,NA)
    x$subsample_ind <- subsample_ind
  } else {
    subsample_by <- 1
    x$subsample_ind <- x_ind
  }
  list(x = x, by = subsample_by)
}

#' Internal function for subsampling
subsample_rename <- function(df, fulldata, colname){
  translate_ind <- dplyr::select(fulldata,x_ind,subsample_ind) %>% dplyr::filter(!is.na(subsample_ind))
  eval_str <- paste("df <- left_join(df,translate_ind, c(\"",colname,"\" = \"subsample_ind\"))",sep="")
  eval(parse(text = eval_str))
  df[,colname] <- df$x_ind
  df$x_ind <- NULL
  df
}
