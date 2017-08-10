#' Prepare HMM output for proper comparison plots
#'
#' \code{prepare_HMM}
#' @export

prepare_HMM <- function(data, hmm.model = NULL, diag.var, order.var = diag.var[1]){
  outputsHMM <- stat_segm_HMM(data = data, hmm.model = hmm.model, diag.var = diag.var, order.var = order.var)
  names(outputsHMM) <- c("segments","states")

  segmented <- list("data" = data,
                    "type" = "behavior",
                    "seg.type" = "HMM",
                    "outputs" = outputsHMM,
                    "model.hmm" = mod.hmm,
                    "Diagnostic variables" = diag.var,
                    "Order variable" = order.var)
  class(segmented) <- "segmentation"
  return(segmented)
}

#' Get segment statistic for HMM model
#'
#' \code{stat_segm_HMM}

stat_segm_HMM <- function(data, hmm.model = NULL, diag.var, order.var = NULL){
  df.segm <- prep_segm_HMM(data,hmm.model)
  data$indice <- 1:nrow(data)
  df.states <- calc_stat_states(data,df.segm,diag.var,order.var)
  return(list(df.segm,df.states))
}

#' Internal function for HMM
#'
#' \code{prep_segm_HMM}
#' @export

prep_segm_HMM <- function(data,hmm.model){
  cluster <- moveHMM::viterbi(hmm.model)
  proba_states <- as.data.frame(moveHMM::stateProbs(hmm.model))
  nstates <- dim(proba_states)[2]
  colnames(proba_states) <- paste("state",1:nstates,sep="")
  df.segm <- data.frame("begin"=1:(nrow(data)-1),
                        "end"=2:(nrow(data)),
                        "state"=cluster[1:(nrow(data)-1)])
  df.segm <- cbind(df.segm,proba_states[1:(nrow(data)-1),])
  return(df.segm)
}
