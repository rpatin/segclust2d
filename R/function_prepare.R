#' Calculate statistics on a given segmentation
#'
#' \code{stat_segm} calculates statistics of a given segmentation : mean and
#' variance of the different states. it also creates standard objects for plot.
#' Actually implementation exists for \code{moveHMM} and \code{segTraj} methods.
#' @param data the data.frame with the different variable
#' @param diag.var names of the variables on which statistics are calculated
#' @param order.var names of the variable with which states are ordered
#' @param seg.type either 'hybrid' or 'dynprog'
#' @param nseg number of segment chosen
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

stat_segm <- function(data, diag.var, order.var = NULL, param = NULL, seg.type = NULL, nseg){
  df.segm <- prep_segm(data,param,nseg = nseg, seg.type = seg.type)
  data$indice <- 1:nrow(data)
  df.states <- calc_stat_states(data,df.segm,diag.var,order.var)
  return(list(df.segm,df.states))
}


#' Find segment and states for a Picard model
#'
#' \code{prep_segm_picard} find the different segment and states of a given HMM
#' model
#' @param data the data.frame with the different variable
#' @param diag.var names of the variables on which statistics are calculated
#' @param order.var names of the variable with which states are ordered
#' @param param the param output of the segmentation
#' @param seg.type either 'hybrid' or 'dynprog'
#' @param nseg number of segment chosen
#' @return a data.frame with states of the different segments
#'
#' @examples
#' prep_segm_picard(data,picard.param,picard.type='hybrid',picard.nseg=NULL)

# attributes(subdf2$scaled_speed)<- NULL
# attributes(subdf2$scaled_angle)<- NULL
# outputs <- segtools::stat_segm(subdf2,diag.var,order.var,model.type='picard',picard.param=param,picard.type='hybrid')

prep_segm <- function(data,param,seg.type=NULL,nseg=NULL){

  if(seg.type=="segclust"){
    df.segm <- as.data.frame(param$rupt)
    colnames(df.segm) <- c("begin","end")
    df.segm$state <- param$cluster
    tmp.tau <- as.data.frame(param$tau)
    nstates <- dim(tmp.tau)[2]
    colnames(tmp.tau) <- paste("state",1:nstates,sep="")
    df.segm <- cbind(df.segm,tmp.tau)
    return(df.segm)
  } else {
    rupt = param$t.est[nseg,1:nseg]
    if(nseg == 1) {
      df.segm <- data.frame(begin=c(1),end=rupt[1],state=1)
    } else {
      df.segm <- data.frame(begin=c(1,rupt[1:(nseg-1)]),end=rupt,state=1:nseg)
    }
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
  df.states$state_ordered  <- rank(df.states[,paste("mu",order.var[1],sep=".")])
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
#' @export

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


#' Calculate BIC
#'
#' \code{BIC} calculates BIC given log-likelihood, number of segment and number of class
#' @param likelihood log-likelihood
#' @param ncluster number of cluster
#' @param nseg number of segment
#' @return a data.frame with BIC, number of cluster and number of segment
#'
#' @examples
#' calc_stat_states(data,diag.var=c("dist","angle"),order.var='dist',type='hmm',hmm.model=mod1.hmm)
#' @export

calc_BIC <- function(likelihood,ncluster,nseg,n){
  BIC = likelihood - 0.5*(5*ncluster-1)*log(2*n) - 0.5 * nseg * log(2*n)
  return(data.frame(BIC=BIC,ncluster=ncluster,nseg=nseg))
}
