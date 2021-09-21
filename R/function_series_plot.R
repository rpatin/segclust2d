#' Plot segmentation on time-serie
#'
#' \code{plot_segm} plot segmented time serie.
#' @param data the data.frame with the different variable
#' @param diag.var names of the variables on which statistics are calculated
#' @param output outputs of the segmentation  or segclust algorithm for one
#'   number of segment
#' @param interactive should graph be interactive through leaflet ?
#' @param html should the graph be incorporated in a markdown file through
#'   htmltools::tagList()
#' @param order should cluster be ordered
#' @param stationarity if TRUE, cut each segment in three and plot each part
#'   with its own mean to assess stationarity of each segment
#' @param x_col column name for time
#' @return a graph
#'
#' @importFrom magrittr "%>%"
#' @export
#' @examples 
#' \dontrun{
#' #res.segclust is the results of the segmentation-clustering algorithm
#' ncluster = 3
#' nseg = 10
#'  g <- plot_segm(data = res.segclust$data, output =
#'   res.segclust$outputs[[paste(ncluster,"class -",nseg, "segments")]], 
#'    diag.var = x$`Diagnostic variables`,x_col = 'dateTime)
#' #res.seg is the results of the segmentation-only algorithm
#' nseg = 10
#'  g <- plot_segm(data = res.segclust$data, 
#'  output = res.segclust$outputs[[paste(nseg, "segments")]], 
#'   diag.var = x$`Diagnostic variables`,x_col = 'dateTime)
#'  
#' }

plot_segm <- 
  function(data, output,
           interactive=FALSE, 
           diag.var, x_col="expectTime",
           html=FALSE, order=FALSE,
           stationarity = FALSE){

  data$indice <- seq_len(nrow(data))
  df.states <- output[[2]]
  df.segm <- output[[1]]
  state_color <- ifelse(order,"state_ordered","state")
  prepMu <- find_mu_sd(df.states,diag.var)
  # 1 seule segmentation
  chosen_segmentation <- dplyr::left_join(df.segm,prepMu, by = c("state"))
  data.melt <- reshape2::melt(data,measure.vars = diag.var)
  chosen_segmentation$begin_date <- data[chosen_segmentation$begin,x_col]
  chosen_segmentation$end_date <- data[chosen_segmentation$end,x_col]
  if(stationarity){
    df_stat <- NULL
    for(seg in seq_len(nrow(df.segm))){
      if(df.segm[seg,'end']-df.segm[seg,'begin'] < 3){
        df_stat <- rbind(df_stat,df.segm[seg,])
      } else {
        tmp <- rbind(df.segm[seg,],df.segm[seg,],df.segm[seg,])
        begin <- df.segm[seg,'begin']
        end <- df.segm[seg,'end']
        tmp$end[1] <- begin + round((end-begin) / 3) - 1
        tmp$begin[2] <- begin + round((end-begin) / 3)
        tmp$end[2] <- begin + round((end-begin) * 2 / 3) -1
        tmp$begin[3] <- begin + round((end-begin) * 2 / 3)
        df_stat <- rbind(df_stat,tmp)
      }
    }
    df_stat$state <- seq_len(nrow(df_stat))
    df_stat_states <- 
      calc_stat_states(data, df_stat,
                       diag.var = diag.var, 
                       order.var = diag.var[1])
    prepMu_stat <- find_mu_sd(df_stat_states,diag.var)
    segmentation_stat <- 
      dplyr::left_join(df_stat,prepMu_stat, by = c("state"))
    segmentation_stat$begin_date <-
      data[segmentation_stat$begin,x_col]
    segmentation_stat$end_date <- 
      data[segmentation_stat$end,x_col]
  }
  # separate = T
  if(!(interactive)){
    # g <-
    # ggplot2::ggplot(data.melt)+
    # ggplot2::geom_line(ggplot2::aes_string(x=x_col,y="value"))+
    # ggplot2::facet_wrap(~variable,ncol=1,scales="free_y")+
    # ggplot2::geom_rect(data=segmentation,
    # ggplot2::aes_string(xmin="begin_date",
    # xmax="end_date",ymin="mu-sd",ymax="mu+sd",
    # fill=paste("factor(",state_color,")",sep="")),alpha=0.2)+
    # ggplot2::geom_segment(data=segmentation,
    # ggplot2::aes_string(x="begin_date",
    # xend="end_date",y="mu",yend="mu",
    # col=paste("factor(",state_color,")")))


    g <-   ggplot2::ggplot(data.melt)+
      ggplot2::geom_line(ggplot2::aes_string(x=x_col,y="value"),col='grey60')+
      ggplot2::facet_wrap(~variable,ncol=1,scales="free_y")+
      ggplot2::geom_rect(
        data=chosen_segmentation,
        ggplot2::aes_string(xmin="begin_date",
                            xmax="end_date",
                            ymin="mu-sd",
                            ymax="mu+sd",
                            fill=paste0("factor(",state_color,")")),
        alpha=0.2)+
      ggplot2::geom_segment(
        data=chosen_segmentation,
        ggplot2::aes_string(
          x="begin_date",
          xend="end_date",
          y="mu",
          yend="mu",
          col=paste0("factor(",state_color,")")))+
      ggplot2::theme_bw()
    
    if(stationarity){
      g <- g +
        ggplot2::geom_segment(
          data=segmentation_stat,
          ggplot2::aes_string(x="begin_date",
                              xend="end_date",
                              y="mu",
                              yend="mu"),
          col="black",size = .4)
    }
    df.label <- data.frame()

    return(g)
  } else {
    data$state <- df.segm[
      findInterval(data$indice,
                   df.segm$begin,
                   rightmost.closed = FALSE,
                   left.open = FALSE),
      "state"]
    data_seg <- dplyr::left_join(data,df.states, by = "state")
    a <- NULL
    for(var in diag.var){
      muvar  <- paste("mu",var,sep=".")
      data.xts <- xts::xts(data_seg[,c(var,muvar)],
                           order.by=data_seg[,x_col])
      if(html){
        b <- htmltools::tagList(
          dygraphs::dygraph(data.xts, main = var, 
                            group = "Lionsteak",height=300) %>%
            dygraphs::dySeries(var, drawPoints = FALSE,
                               color = "grey") %>%
            dygraphs::dySeries(muvar, stepPlot = FALSE, 
                               fillGraph = FALSE,
                               color = "red",strokeWidth=2))
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
          dygraphs::dySeries(muvar, stepPlot = FALSE,
                             fillGraph = FALSE, color = "red",strokeWidth=2)
      }
    }
    return(a)
  }
}
