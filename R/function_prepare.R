#' Calculate statistics on a given segmentation
#'
#' \code{stat_segm} calculates statistics of a given segmentation : mean and
#' variance of the different states. it also creates standard objects for plot.
#' Actually implementation exists for \code{moveHMM} and \code{segTraj} methods.
#' @param data the data.frame with the different variable
#' @param diag.var names of the variables on which statistics are calculated
#' @param order.var names of the variable with which states are ordered
#' @param type type of model either 'hmm' or 'picard'
#' @param hmm.model the moveHMM::fitHMM output
#' @param picard.param the param output of the segmentation
#' @param picard.type either 'hybrid' or 'dynprog'
#' @param picard.nseg number of segment chosen
#' @return  a list which first element is a data.frame with states of the
#'   different segments and which second element is a data.frame with mean and
#'   variance of the different states
#'
#' @examples
#' stat_segm(data,diag.var=c("dist","angle"),order.var='dist',type='hmm',hmm.model=mod1.hmm)
#' @export
#'
#
# test <- stat_segm(data=subdf2,diag.var,order.var,model.type='picard',picard.param=param,picard.type='hybrid')
# attributes(data$scaled_speed)<- NULL
# attributes(data$scaled_angle)<- NULL
# plot_states(outputs,diag.var)

stat_segm <- function(data, diag.var, order.var = NULL, model.type = 'hmm', hmm.model = NULL, picard.param = NULL, picard.type = NULL, picard.nseg){
  if(model.type == 'hmm'){
    df.segm <- prep_segm_HMM(data,hmm.model)
  } else if (model.type == 'picard') {
    df.segm <- prep_segm_picard(data,picard.param,picard.type,picard.nseg)
  } else {
    stop("formal argument type not recognized")
  }
  data$indice <- 1:nrow(data)
  df.states <- calc_stat_states(data,df.segm,diag.var,order.var)
  return(list(df.segm,df.states))
}

#' Find segment and states for a HMM model
#'
#' \code{prep_segm_HMM} find the different segment and states of a given HMM
#' model
#' @param data the data.frame with the different variable
#' @param hmm.model the moveHMM::fitHMM output
#' @return  a data.frame with states of the different segments
#'
#' @examples
#' prep_segm_HMM(data,hmm.model=mod1.hmm)

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

#' Find segment and states for a Picard model
#'
#' \code{prep_segm_picard} find the different segment and states of a given HMM
#' model
#' @param data the data.frame with the different variable
#' @param diag.var names of the variables on which statistics are calculated
#' @param order.var names of the variable with which states are ordered
#' @param type type of model either 'hmm' or 'picard'
#' @param hmm.model the moveHMM::fitHMM output
#' @param picard.param the param output of the segmentation
#' @param picard.type either 'hybrid' or 'dynprog'
#' @param picard.nseg number of segment chosen
#' @return  a data.frame with states of the different segments
#'
#' @examples
#' prep_segm_picard(data,picard.param,picard.type='hybrid',picard.nseg=NULL)

# attributes(subdf2$scaled_speed)<- NULL
# attributes(subdf2$scaled_angle)<- NULL
# outputs <- segtools::stat_segm(subdf2,diag.var,order.var,model.type='picard',picard.param=param,picard.type='hybrid')

prep_segm_picard <- function(data,picard.param,picard.type='hybrid',picard.nseg=NULL){

  if(picard.type=="hybrid"){
    df.segm <- as.data.frame(picard.param$rupt)
    colnames(df.segm) <- c("begin","end")
    df.segm$state <- picard.param$cluster
    tmp.tau <- as.data.frame(picard.param$tau)
    nstates <- dim(tmp.tau)[2]
    colnames(tmp.tau) <- paste("state",1:nstates,sep="")
    df.segm <- cbind(df.segm,tmp.tau)
    return(df.segm)
  } else {
    rupt = picard.param$t.est[picard.nseg,1:picard.nseg]
    df.segm <- data.frame(begin=c(1,rupt[1:(picard.nseg-1)]),end=rupt,state=1:picard.nseg)
    return(df.segm)
  }
}


#' Calculate state statistics
#'
#' \code{calc_stat_states} calculates statistics of a given segmentation : mean
#' and variance of the different states.
#' @param data the data.frame with the different variable
#' @param diag.var names of the variables on which statistics are calculated
#' @param order.var names of the variable with which states are ordered
#' @param df.segm output of prep_segm function
#' @return  a data.frame with mean and variance of the different states
#'
#' @examples
#' calc_stat_states(data,diag.var=c("dist","angle"),order.var='dist',type='hmm',hmm.model=mod1.hmm)
#' @importFrom magrittr "%>%"
#' @export

calc_stat_states <- function(data,df.segm,diag.var,order.var=NULL)
{
  data <- dplyr::mutate(data,state= df.segm[findInterval(indice,df.segm$begin,rightmost.closed = F,left.open = F),"state"])

  eval_str <- paste("dplyr::group_by(data,state) %>% dplyr::summarise(prop=n()/nrow(data),",paste("mu.",diag.var," = mean(",diag.var,",na.rm=T)",collapse=",",sep=""),",",paste("sd.",diag.var," = sd(",diag.var,",na.rm=T)",collapse=",",sep=""),") %>% as.data.frame()",sep="")
  df.states <- eval(parse(text=eval_str))
  df.states$state_ordered  <- order(df.states[,paste("mu",order.var,sep=".")])
  return(df.states)
}


#' Find mean and standard deviation of segments
#'
#' \code{find_mu_sd} calculates statistics of a given segmentation : mean
#' and variance of the different states.
#' @param df.states a list of data.frame
#' @param diag.var names of the variables on which statistics are calculated
#' @param order.var names of the variable with which states are ordered
#' @param df.segm output of prep_segm function
#' @return  a data.frame with mean and variance of the different states
#'
#' @examples
#' calc_stat_states(data,diag.var=c("dist","angle"),order.var='dist',type='hmm',hmm.model=mod1.hmm)

find_mu_sd <- function(df.states,diag.var){
  if(is.null(df.states$model)) df.states$model <- 'model'
  eval_str <-  paste("reshape2::melt(df.states,measure.var = c(", paste("\"mu.",diag.var,"\"",collapse=",",sep=""),"))",sep="")
    mu.melt <- eval(parse(text=eval_str))
    mu.melt$variable <- plyr::laply(strsplit(as.character(mu.melt$variable),split=".",fixed=T),function(x){paste(x[-1],collapse = ".")})
    mu.melt <- dplyr::rename(mu.melt,mu=value)
    mu.melt <- dplyr::select(mu.melt,state,state_ordered,variable,mu,prop,model)

    eval_str <-  paste("reshape2::melt(df.states,measure.var = c(", paste("\"sd.",diag.var,"\"",collapse=",",sep=""),"))",sep="")
    sd.melt <- eval(parse(text=eval_str))
    sd.melt$variable <- plyr::laply(strsplit(as.character(sd.melt$variable),split=".",fixed=T),function(x){paste(x[-1],collapse = ".")})
    sd.melt <- dplyr::rename(sd.melt,sd=value)
    sd.melt <- dplyr::select(sd.melt,state,state_ordered,variable,sd,prop,model)

    mu.melt <- dplyr::left_join(mu.melt,sd.melt,by = c("state", "prop", "state_ordered", "variable","model"))
  return(mu.melt)
}


