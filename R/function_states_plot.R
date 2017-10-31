#' Plot states statistics
#'
#' \code{plot_states} plot states statistics
#' @param diag.var names of the variables on which statistics are calculated
#' @param position_width width between different model if several models are compared
#' @param output outputs of the segmentation  or segclust algorithm for one number of segment
#' @param order should cluster be ordered
#' @return a graph
#'
#' @examples
#' stat_segm(data,diag.var=c("dist","angle"),order.var='dist',type='hmm',hmm.model=mod1.hmm)
#' @export

plot_states <- function(outputs,diag.var, position_width=0.3,order = F){
  state_variable <- ifelse(order,"factor(state_ordered)","factor(state)")

  df.states <- outputs[[2]]
  if(is.null(df.states$model)) df.states$model <- 'model'
  mu.list <- find_mu_sd(df.states,diag.var)

  if(class(mu.list)== 'data.frame'){
    g <- ggplot2::ggplot(data=mu.list, ggplot2::aes_string(x=state_variable,y="mu"))+
      ggplot2::geom_point()+
      ggplot2::facet_wrap(~variable,scale="free")+
      ggplot2::geom_errorbar(ggplot2::aes(ymin=mu-sd,ymax=mu+sd),width=0.1)+
      ggplot2::xlab("State")+
      ggplot2::ylab("Distribution (mean +/- sd)")
  } else {
    mu.merged <- do.call("rbind",mu.list)

    g <- ggplot2::ggplot(data=mu.merged, ggplot2::aes_string(x=state_variable,y="mu",col="model"))+
      ggplot2::geom_point(position = ggplot2::position_dodge(width=position_width))+
      ggplot2::facet_wrap(~variable,scale="free")+
      ggplot2::geom_errorbar(ggplot2::aes(ymin=mu-sd,ymax=mu+sd),width=0.1,position=ggplot2::position_dodge(width=position_width))+
      xlab("State")+
      ylab("Distribution (mean +/- sd)")

  }
  return(g)
}
