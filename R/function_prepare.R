#' Calculate statistics on a given segmentation
#'
#' \code{stat_segm} calculates statistics of a given segmentation : mean and
#' variance of the different states. it also creates standard objects for plot.
#' @param data the data.frame with the different variable
#' @param diag.var names of the variables on which statistics are calculated
#' @param order.var names of the variable with which states are ordered
#' @param seg.type either 'hybrid' or 'dynprog'
#' @param nseg number of segment chosen
#' @param param parameters of ouptput segmentation
#' @return  a list which first element is a data.frame with states of the
#'   different segments and which second element is a data.frame with mean and
#'   variance of the different states
#'
#' @examples
#' \dontrun{
#' #res.segclust is a result of a segmentation-clustering algorithm
#' param <- res.segclust$param[["3 class"]]
#' nseg = 10
#' out <- stat_segm(data, diag.var = c("dist","angle"),
#'  order.var = "dist", param = param, nseg=nseg, seg.type = "segclust")
#' 
#' }
#' @export
#'

stat_segm <- function(data, diag.var, order.var = NULL, param = NULL, seg.type = NULL, nseg){
  subdata <- data[!is.na(data$subsample_ind),]
  df.segm <- prep_segm(subdata,param,nseg = nseg, seg.type = seg.type)

  subdata$indice <- 1:nrow(subdata)
  df.states <- calc_stat_states(subdata,df.segm,diag.var,order.var)
  df.segm <- subsample_rename(df.segm,data,"begin")
  df.segm <- subsample_rename(df.segm,data,"end")

  return(list(df.segm,df.states))
}


#' Find segment and states for a Picard model
#'
#' \code{prep_segm} find the different segment and states of a given HMM
#' model
#' @param data the data.frame with the different variable
#' @param param the param output of the segmentation
#' @param seg.type either 'hybrid' or 'dynprog'
#' @param nseg number of segment chosen
#' @return a data.frame with states of the different segments
#'

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
      df.segm <- data.frame(begin=c(1,rupt[1:(nseg-1)]+1),end=rupt,state=1:nseg)
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
#' \dontrun{calc_stat_states(data, diag.var = c("dist","angle"),
#' order.var='dist', type='hmm',hmm.model=mod1.hmm)}
#' @importFrom magrittr "%>%"
#' @export

calc_stat_states <- function(data,df.segm,diag.var,order.var=NULL)
{
  data$state <- df.segm[findInterval(data$indice,df.segm$begin,rightmost.closed = F,left.open = F),"state"]

  
  tmp.gp <- dplyr::group_by(data,state)
  eval_str <- paste("dplyr::summarise(tmp.gp, prop=dplyr::n()/nrow(data),",
                    paste("mu.",diag.var," = mean(",diag.var,",na.rm=T)",collapse=",",sep=""),
                    ",",
                    paste("sd.",diag.var," = stats::sd(",diag.var,",na.rm=T)",collapse
                          =",",sep=""), ")")
  tmp.states <- eval(parse(text=eval_str))
  df.states <- as.data.frame(tmp.states)
  df.states$state_ordered  <- rank(df.states[,paste("mu",order.var[1],sep=".")])
  return(df.states)
}


#' Find mean and standard deviation of segments
#'
#' \code{find_mu_sd} calculates statistics of a given segmentation : mean
#' and variance of the different states.
#' @param df.states a list of data.frame
#' @param diag.var names of the variables on which statistics are calculated
#' @return  a data.frame with mean and variance of the different states
#'
#' @export

find_mu_sd <- function(df.states,diag.var){
  if(is.null(df.states$model)) df.states$model <- 'model'
  var_measure <- c(paste("mu.",diag.var,sep=""))
  mu.melt <-  reshape2::melt(df.states,measure.var = var_measure)
    mu.melt$variable <- plyr::laply(strsplit(as.character(mu.melt$variable),split=".",fixed=T),function(x){paste(x[-1],collapse = ".")})
    mu.melt$mu <- mu.melt$value
    mu.melt$value <- NULL
    mu.melt <- data.frame("state" = mu.melt$state,
                          "state_ordered" = mu.melt$state_ordered,
                          "variable" = mu.melt$variable,
                          "mu" = mu.melt$mu,
                          "prop" = mu.melt$prop,
                          "model" = mu.melt$model)

    var_measure <- c(paste("sd.",diag.var,sep=""))
    sd.melt <-  reshape2::melt(df.states,measure.var = var_measure)
    sd.melt$variable <- plyr::laply(strsplit(as.character(sd.melt$variable),split=".",fixed=T),function(x){paste(x[-1],collapse = ".")})
    sd.melt$sd <- sd.melt$value
    sd.melt$value <- NULL

    # sd.melt <- with(sd.melt(data.frame(state,state_ordered,variable,sd,prop,model)))
    sd.melt <- data.frame("state" = sd.melt$state,
                          "state_ordered" = sd.melt$state_ordered,
                          "variable" = sd.melt$variable,
                          "sd" = sd.melt$sd,
                          "prop" = sd.melt$prop,
                          "model" = sd.melt$model)

    mu.melt <- dplyr::left_join(mu.melt,sd.melt,by = c("state", "prop", "state_ordered", "variable","model"))
  return(mu.melt)
}


#' Calculate BIC
#'
#' \code{BIC} calculates BIC given log-likelihood, number of segment and number of class
#' @param likelihood log-likelihood
#' @param ncluster number of cluster
#' @param nseg number of segment
#' @param n number of observations
#' @return a data.frame with BIC, number of cluster and number of segment
#'
#' @export

calc_BIC <- function(likelihood,ncluster,nseg,n){
  BIC = likelihood - 0.5*(5*ncluster-1)*log(2*n) - 0.5 * nseg * log(2*n)
  return(data.frame(BIC=BIC,ncluster=ncluster,nseg=nseg))
}

#' Check for repetition in the series
#'
#' \code{check_repetition} checks whether the series have identical or near-identical repetition larger than lmin.
#' if that is the case, throw an error, the algorithm cannot yet handle these repetition,
#' because variance on the segment would be null. 
#' @param x the bivariate series to be tested
#' @param lmin minimum length of segment
#' @param rounding whether or not series are rounded
#' @param magnitude number of magnitude of standard deviation below which values are rounded. i.e if magnitude = 3, difference smaller than one thousandth of the standard deviation are rounded to the same value.
#' @return a boolean, TRUE if there is any repetition larger or equal to lmin.
#'
#' @export
#' @examples 
#' set.seed(42) 
#' dat <- rbind(base::sample(seq(1,10),  size= 100, replace = TRUE),
#'              base::sample(seq(1,10),  size= 100, replace = TRUE))
#' check_repetition(dat, lmin = 3)
#' check_repetition(dat, lmin = 5)             

check_repetition <- function(x,lmin, rounding = FALSE, magnitude = 3){
    if(rounding){
      sd_x1 <- stats::sd(x[1,])
      magn1 <- - base::floor(log10(sd_x1)) +magnitude
      x1 <- base::round(x[1,], digits = magn1)
      sd_x2 <- stats::sd(x[2,])
      magn2 <- - base::floor(log10(sd_x2)) +magnitude
      x2 <- base::round(x[2,], digits = magn2)
      rep_1 <- rle(x1)
      rep_2 <- rle(x2)
      if( any(rep_1$length >= lmin) || any(rep_2$length >= lmin)){
        return(TRUE)
      } else {
        return(FALSE)
      } 
    } else {
      rep_1 <- rle(x[1,])
      rep_2 <- rle(x[2,])
      if( any(rep_1$length >= lmin) || any(rep_2$length >= lmin)){
        return(TRUE)
      } else {
        return(FALSE)
      } 
    }
   
}

#' Relabel states of a segmentation/clustering output
#' 
#' \code{relabel_states} relabel the states of a segmentation/clustering output.
#' This allows merging different states into the same if for instance several of
#' the model states represent the same behavioural states.
#' @param mode.segclust segclust output
#' @param newlabel a vector with the new names ordered, corresponding to 
#'   state_ordered
#' @param ncluster the number of cluster for which you want relabeling
#' @param nseg the number of segment for which you want relabeling
#' @param order boolean, whether this changes the ordered states or not. FALSE 
#'   value obsolete for now
#' @return a segmentation object with state names changed for the segmentation
#'   specified by ncluster and nseg
#'   
#' @export

relabel_states <- function(mode.segclust, newlabel, ncluster, nseg, order = TRUE){
  tmp <- mode.segclust$outputs[[paste0(ncluster," class - ",nseg," segments")]]  
  tmp$states$state_ordered <- newlabel[tmp$states$state_ordered]
  mode.segclust$outputs[[paste0(ncluster," class - ",nseg," segments")]]  <- tmp
  mode.segclust
}
