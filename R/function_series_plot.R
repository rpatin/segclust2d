#' Plot segmentation on time-serie
#'
#' \code{plot_segm} plot segmented time serie.
#' @param data the data.frame with the different variable
#' @param diag.var names of the variables on which statistics are calculated
#' @param position_width width between different model if several models are compared
#' @return a graph
#'
#' @examples
#' stat_segm(data,diag.var=c("dist","angle"),order.var='dist',type='hmm',hmm.model=mod1.hmm)
#' @importFrom magrittr "%>%"
#' @export
# plot_segm(subdf2,outputs,separate=T,interactive=T,diag.var=diag.var,x_col='date')

plot_segm <- function(data,output,separate=T,interactive=F,diag.var,x_col="expectTime",html=F){
  # if(class(df.states) != "list"){
  #   df.states <- list(df.states)
  # }
  # if(class(df.segm.list) != "list"){
  #   df.segm.list <- list(df.segm.list)
  # }
  #
  # if(length(df.segm.list) != length(df.states)){ message("Error, df.segm.list and df.states of different list size"); stop("Different list size")}
  data$indice <- 1:nrow(data)
  df.states <- output[[2]]
  df.segm <- output[[1]]

  prepMu <- find_mu_sd(df.states,diag.var)
  # 1 seule segmentation
  if(class(prepMu)=='data.frame'){
    segmentation <- dplyr::left_join(df.segm,prepMu, by = c("state"))
    data.melt <- reshape2::melt(data,measure.vars = diag.var)
    segmentation$begin_date <- data[segmentation$begin,x_col]
    segmentation$end_date <- data[segmentation$end,x_col]
    # separate = T
    if(separate & !(interactive)){
      g <-   ggplot2::ggplot(data.melt)+ggplot2::geom_line(ggplot2::aes_string(x=x_col,y="value"))+
        ggplot2::facet_wrap(~variable,ncol=1,scales="free_y")+
        ggplot2::geom_rect(data=segmentation,ggplot2::aes(xmin=begin_date,xmax=end_date,ymin=mu-sd,ymax=mu+sd,fill=factor(state_ordered)),alpha=0.2)+
        ggplot2::geom_segment(data=segmentation,ggplot2::aes(x=begin_date,xend=end_date,y=mu,yend=mu,col=factor(state_ordered)))
      return(g)
    }
    if(separate & interactive){
      data <- dplyr::mutate(data,state= df.segm[findInterval(indice,df.segm$begin,rightmost.closed = F,left.open = F),"state"])
      data_seg <- dplyr::left_join(data,df.states, by = "state")
      a <- NULL
      for(var in diag.var){
        muvar = paste("mu",var,sep=".")
        data.xts <- xts::xts(data_seg[,c(var,muvar)],order.by=data_seg[,x_col])
        if(html){
          b <- htmltools::tagList(dygraphs::dygraph(data.xts, main = var, group = "Lionsteak",height=300) %>%
                                    dygraphs::dySeries(var, drawPoints = FALSE, color = "grey") %>%
                                    dygraphs::dySeries(muvar, stepPlot = FALSE, fillGraph = FALSE, color = "red",strokeWidth=2))
          if(var == dplyr::last(diag.var)){
            b <- b %>% dygraphs::dyRangeSelector()
          }
          if(is.null(a)){
            a <- b
          } else {
            a <- c(a,b)
          }
        } else {
          dygraphs::dygraph(data.xts, main = var, group = "Lionsteak") %>%
            dygraphs::dySeries(var, drawPoints = FALSE, color = "grey") %>%
            dygraphs::dySeries(muvar, stepPlot = FALSE, fillGraph = FALSE, color = "red",strokeWidth=2)
        }
      }
      return(a)
    }
    if(!separate){
      separate = F
      # to be done
      message("Error: Separate=F not coded yet"); stop("Ask the lazy programmer")
    }
  } #end length ==1
}
