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

plot_segm <- function(data,output,separate=T,interactive=F,diag.var,x_col="expectTime",html=F,order=F, stationarity = NULL, mean = NULL, var = NULL){
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
  if(class(prepMu)=='data.frame'){
    segmentation <- dplyr::left_join(df.segm,prepMu, by = c("state"))
    data.melt <- reshape2::melt(data,measure.vars = diag.var)
    segmentation$begin_date <- data[segmentation$begin,x_col]
    segmentation$end_date <- data[segmentation$end,x_col]
    # separate = T
    if(separate & !(interactive)){
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

      if(!is.null(stationarity)){
        colnames(stationarity) <- c("seg",diag.var,"stat_tot")
        stationarity.melt <- reshape2::melt(stationarity,measure.var = diag.var) %>% dplyr::filter(value == "stationary") %>% mutate(variable = as.factor(variable))

        if(nrow(stationarity.melt) >0 ){
          df.label <- dplyr::select(stationarity.melt, seg, variable)
          df.label$label.stat = "s"
        }
      }

      if(!is.null(mean)){
        samemean.x <- which(mean$x == "Same Mean",arr.ind = T) %>% data.frame() %>% mutate(label.mean = paste("m",row,col,sep="")) %>% reshape2::melt(measure.var = c("row","col"),value.name="seg") %>% mutate(variable = "x")
        samemean.y <- which(mean$y == "Same Mean",arr.ind = T) %>% data.frame() %>% mutate(label.mean = paste("m",row,col,sep="")) %>% reshape2::melt(measure.var = c("row","col"),value.name="seg") %>% mutate(variable = "y")
        samemean <- rbind(samemean.x,samemean.y) %>% dplyr::group_by(seg,variable) %>% dplyr::summarise(label.mean = paste(label.mean,collapse = "-")) %>% mutate(variable = as.factor(variable))
        if(nrow(samemean) > 0 ){
          if(nrow(df.label) > 0){
            df.label <- dplyr::full_join(df.label, samemean, by = c("variable","seg"))
          } else {
            df.label <- dplyr::select(samemean, seg, variable, label.mean)
            df.label$label.stat = ""
          }
        } else {
          df.label$label.mean = ""
        }
      } else {
        df.label$label.mean = ""
      }

      if(!is.null(var)){
        samevar.x <- which(var$x == "Same Variance",arr.ind = T) %>% data.frame() %>% mutate(label.var = paste("v",row,col,sep="")) %>% reshape2::melt(measure.var = c("row","col"),value.name="seg") %>% mutate(variable = "x")
        samevar.y <- which(var$y == "Same Variance",arr.ind = T) %>% data.frame() %>% mutate(label.var = paste("v",row,col,sep="")) %>% reshape2::melt(measure.var = c("row","col"),value.name="seg") %>% mutate(variable = "y")
        samevar <- rbind(samevar.x,samevar.y) %>% dplyr::group_by(seg,variable) %>% dplyr::summarise(label.var = paste(label.var,collapse = "-")) %>% mutate(variable = as.factor(variable))
        if(nrow(samevar) > 0 ){
          if(nrow(df.label) > 0){
            df.label <- dplyr::full_join(df.label, samevar, by = c("variable","seg"))
          } else {
            df.label <- dplyr::select(samevar, seg, variable, label.var)
            df.label$label.stat = ""
            df.label$label.mean = ""
          }
        } else {
          df.label$label.var = ""
        }
      } else {
        df.label$label.var = ""
      }

      if(nrow(df.label) > 0){
        df.label$label.stat[is.na(df.label$label.stat)] <- ''
        df.label$label.mean[is.na(df.label$label.mean)] <- ''
        df.label$label.var[is.na(df.label$label.var)] <- ''
        df.label$label <- mapply(function(x,y,z){paste(x,y,z,collapse=" ")},df.label$label.stat,df.label$label.mean,df.label$label.var)
        range_x <- diff(range(data[,diag.var[1]]))/10
        range_y <- diff(range(data[,diag.var[2]]))/10
        segmentation <- mutate(segmentation,variable = factor(variable))
        df.label.plot <- dplyr::left_join(df.label,segmentation,by = c("variable"="variable", "seg"="state")) %>% mutate(mean_date = (begin_date+end_date)/2)
        df.label.plot.x <- dplyr::filter(df.label.plot,variable == diag.var[1])
        df.label.plot.y <- dplyr::filter(df.label.plot,variable == diag.var[2])
        g <- g + geom_text(data = df.label.plot.x, aes(label = label,x = mean_date, y = mu+sd, col = factor(seg)),nudge_x = 0, nudge_y = range_x)
        g <- g + geom_text(data = df.label.plot.y, aes(label = label,x = mean_date, y = mu+sd, col = factor(seg)),nudge_x = 0, nudge_y = range_y)

      }


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
