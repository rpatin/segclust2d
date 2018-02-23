#' Prepare shiftfit output for proper comparison plots
#'
#' \code{prepare_shiftfit}
#' @param data data
#' @param shiftfit.model shiftfit.model
#' @param diag.var diag.var
#' @param order.var order.var
#' @export
#' @examples 
#' \dontrun{
#' data(simulshift)
#' # 1. subsample to a reasonable size
#' subdata <- simulshift[seq(1,30000,by = 100),]
#' # 2. use algorithm from marcher package
#' MWN.fit <- with(subdata, marcher::estimate_shift(T=indice, X=x, Y=y,n.clust = 3))
#' # 3. convert output
#' MWN.segm <- prepare_shiftfit(subdata,MWN.fit,diag.var = c("x","y"))
#' # 4. use segclust2d functions
#' plot(MWN.segm)
#' plot(MWN.segm,stationarity = TRUE)
#' segmap(MWN.segm)
#' }
prepare_shiftfit <- function(data, shiftfit.model = NULL, diag.var, order.var = diag.var[1]){
  outputsshift <- stat_segm_shiftfit(data = data, shiftfit.model = shiftfit.model, diag.var = diag.var, order.var = order.var)
  names(outputsshift) <- c("segments","states")

  segmented <- list("data" = data,
                    "type" = "home-range",
                    "seg.type" = "shiftfit",
                    "outputs" = outputsshift,
                    "shiftfit.model" = shiftfit.model,
                    "Diagnostic variables" = diag.var,
                    "Order variable" = order.var)
  class(segmented) <- "segmentation"
  return(segmented)
}

#' Get segment statistic for shiftfit model
#'
#' \code{stat_segm_shiftfit}
#' @inheritParams prepare_shiftfit

stat_segm_shiftfit <- function(data, shiftfit.model = NULL, diag.var, order.var = NULL){
  df.segm <- prep_segm_shiftfit(data,shiftfit.model)
  data$indice <- 1:nrow(data)
  df.states <- calc_stat_states(data,df.segm,diag.var,order.var)
  return(list(df.segm,df.states))
}

#' Internal function for HMM
#'
#' \code{prep_segm_shiftfit}
#' @inheritParams prepare_shiftfit
#' @export

prep_segm_shiftfit <- function(data,shiftfit.model){
  ncluster <- shiftfit.model$n.clust
  df.segm = NULL
  for(i in 1:ncluster){
    if(i == 1){
      begin = 1
    } else {
      begin = round(shiftfit.model$p.hat[paste("t",i-1,sep="")])
      begin <- which.min(abs(shiftfit.model$T - begin))
    }

    if(i == ncluster){
      end = nrow(data)
    } else {
      end = round(shiftfit.model$p.hat[paste("t",i,sep="")]-1)
      end <- which.min(abs(shiftfit.model$T - end))
    }

    df.segm = rbind(df.segm,data.frame("begin"=begin,
                                       "end"=end,
                                       "state"=i))
  }
  return(df.segm)
}
