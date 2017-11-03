#' Plot segmentation on time-serie
#'
#' \code{plot_segm} plot segmented time serie.
#' @param data the data.frame with the different variable
#' @param diag.var names of the variables on which statistics are calculated
#' @param output outputs of the segmentation  or segclust algorithm for one number of segment
#' @param interactive should graph be interactive through leaflet ?
#' @param html should the graph be incorporated in a markdown file through htmltools::tagList()
#' @param order should cluster be ordered
#' @param x_col column name for time
#' @return a graph
#'
#' @examples
#' \dontrun{stat_segm(data,diag.var=c("dist","angle"), order.var='dist',
#' type='hmm',hmm.model=mod1.hmm)}
#' @importFrom magrittr "%>%"
#' @export
# plot_segm(subdf2,outputs,separate=T,interactive=T,diag.var=diag.var,x_col='date')

plot_segm <- function(data,output,interactive=F,diag.var,x_col="expectTime",html=F,order=F){
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
  state_color <- ifelse(order,"state_ordered","state")
  prepMu <- find_mu_sd(df.states,diag.var)
  # 1 seule segmentation
  segmentation <- dplyr::left_join(df.segm,prepMu, by = c("state"))
  data.melt <- reshape2::melt(data,measure.vars = diag.var)
  segmentation$begin_date <- data[segmentation$begin,x_col]
  segmentation$end_date <- data[segmentation$end,x_col]
  # separate = T
  if(!(interactive)){
    # g <-   ggplot2::ggplot(data.melt)+ggplot2::geom_line(ggplot2::aes_string(x=x_col,y="value"))+
    #   ggplot2::facet_wrap(~variable,ncol=1,scales="free_y")+
    #   ggplot2::geom_rect(data=segmentation,ggplot2::aes_string(xmin="begin_date",xmax="end_date",ymin="mu-sd",ymax="mu+sd",fill=paste("factor(",state_color,")",sep="")),alpha=0.2)+
    #   ggplot2::geom_segment(data=segmentation,ggplot2::aes_string(x="begin_date",xend="end_date",y="mu",yend="mu",col=paste("factor(",state_color,")")))


    g <-   ggplot2::ggplot(data.melt)+
      ggplot2::geom_line(ggplot2::aes_string(x=x_col,y="value"))+
      ggplot2::facet_wrap(~variable,ncol=1,scales="free_y")+
      ggplot2::geom_rect(data=segmentation,ggplot2::aes_string(xmin="begin_date",xmax="end_date",ymin="mu-sd",ymax="mu+sd",fill=paste("factor(",state_color,")",sep="")),alpha=0.2)+
      ggplot2::geom_segment(data=segmentation,ggplot2::aes_string(x="begin_date",xend="end_date",y="mu",yend="mu",col=paste("factor(",state_color,")")))

    df.label <- data.frame()

    return(g)
  } else {
    data$state= df.segm[findInterval(data$indice,df.segm$begin,rightmost.closed = F,left.open = F),"state"]
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
}
