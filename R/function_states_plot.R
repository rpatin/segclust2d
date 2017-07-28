#' Plot states statistics
#'
#' \code{plot_states} plot states statistics
#' @param data the data.frame with the different variable
#' @param diag.var names of the variables on which statistics are calculated
#' @param position_width width between different model if several models are compared
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
    g <- ggplot(data=mu.list,aes_string(x=state_variable,y="mu"))+geom_point()+facet_wrap(~variable,scale="free")+geom_errorbar(aes(ymin=mu-sd,ymax=mu+sd),width=0.1)+xlab("State")+ylab("Distribution (mean +/- sd)")
  } else {
    mu.merged <- do.call("rbind",mu.list)

    g <- ggplot(data=mu.merged,aes_string(x=state_variable,y="mu",col="model"))+geom_point(position=position_dodge(width=position_width))+facet_wrap(~variable,scale="free")+geom_errorbar(aes(ymin=mu-sd,ymax=mu+sd),width=0.1,position=position_dodge(width=position_width))+xlab("State")+ylab("Distribution (mean +/- sd)")

  }
  return(g)
}
