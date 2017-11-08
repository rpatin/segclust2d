#' Prepare depmixS4 output for proper plots
#'
#' \code{prepare_depmixS4}
#'
#' @param data data
#' @param depmixS4.model depmixS4.model
#' @param diag.var diag.var
#' @param order.var order.var
#' @export

prepare_depmixS4 <- function(data, depmixS4.model = NULL, diag.var, order.var = diag.var[1]){

  outputsshift <- stat_segm_depmixS4(data = data, depmixS4.model = depmixS4.model, diag.var = diag.var, order.var = order.var)
  names(outputsshift) <- c("segments","states")

  segmented <- list("data" = data,
                    "type" = "behavior",
                    "seg.type" = "depmixS4",
                    "outputs" = outputsshift,
                    "depmixS4.model" = depmixS4.model,
                    "Diagnostic variables" = diag.var,
                    "Order variable" = order.var)
  class(segmented) <- "segmentation"
  return(segmented)
}

#' Get segment statistic for HMM model
#'
#' \code{stat_segm_depmixS4}
#'
#' @inheritParams prepare_depmixS4

stat_segm_depmixS4 <- function(data, depmixS4.model = NULL, diag.var, order.var = NULL){
  df.segm <- prep_segm_depmixS4(depmixS4.model)
  data$indice <- 1:nrow(data)
  df.states <- calc_stat_states(data,df.segm,diag.var,order.var)
  return(list(df.segm,df.states))
}

#' Internal function for HMM
#'
#' \code{prep_segm_depmixS4}
#' @inheritParams prepare_depmixS4
#' @export

prep_segm_depmixS4 <- function(depmixS4.model){
  df.segm = depmixS4::posterior(depmixS4.model)
  tmp = data.frame()
  j = 1
  i = 1
  while( i < nrow(df.segm)){
    i<- i+1
    if( df.segm[i,"state"] != df.segm[j,"state"] ){
      tmp <- rbind(tmp, data.frame("state" = df.segm[j,"state"],begin = j, end = i ))
      j = i
    }
  }
  tmp <- rbind(tmp, data.frame("state" = df.segm[j,"state"],begin = j, end = i))

  # df.segm$begin <- 1:(nrow(df.segm)-1)
  # df.segm$end <- 2:(nrow(df.segm))
  return(tmp)
}
