#' Plot map of segmentation
#'
#' \code{plot_segm} plot segmented traject on a map.
#' @param data the data.frame with the different variable
#' @param diag.var names of the variables on which statistics are calculated
#' @param position_width width between different model if several models are compared
#' @return a graph
#'
#' @examples
#' stat_segm(data,diag.var=c("dist","angle"),order.var='dist',type='hmm',hmm.model=mod1.hmm)
#' @importFrom magrittr "%>%"
#' @export

map_segm <- function(data,output,interactive=F,x_col="expectTime",html=F,scale=100,UTMstring="+proj=utm +zone=35 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0",width=400,height=400,order=NULL,pointsize = 1, linesize = 0.5){
  # print("test")
  df.segm <- dplyr::left_join(output[[1]],output[[2]],by="state")
  data$indice <- 1:nrow(data)

  data <- dplyr::mutate(data,state = df.segm[findInterval(indice,df.segm$begin,rightmost.closed = F,left.open = F),ifelse(order,"state_ordered","state")])

  data <- dplyr::mutate(data,x=x/scale,y=y/scale)
  if(!interactive){
    g <- ggplot2::ggplot(data,ggplot2::aes(x=x,y=y))+
      ggplot2::geom_path(size = linesize)+
      ggplot2::geom_point(ggplot2::aes(col=factor(state)),size = pointsize)
    return(g)
  } else {
    # coordinates(data) <- c("x","y")

    LinesList <- list()
    for(i in 1:nrow(df.segm))
    {
      begin <- df.segm$begin[i]
      end <- df.segm$end[i]
      X <- data$x[begin:end]
      Y <- data$y[begin:end]
      LinesList[[i]]<- sp::Lines(sp::Line(cbind(X,Y)),ID=i)
    }

    ConnectionList <- list()
    for(i in 1:(nrow(df.segm)-1))
    {
      begin <- df.segm$end[i]
      X <- data$x[begin:(begin+1)]
      Y <- data$y[begin:(begin+1)]
      ConnectionList[[i]]<- sp::Lines(sp::Line(cbind(X,Y)),ID=i)
    }

    list_colors = RColorBrewer::brewer.pal(n = length(unique(df.segm$state)),name='Set1')

    dfSPLDF <- sp::SpatialLinesDataFrame(sp::SpatialLines(LinesList,proj4string = sp::CRS(UTMstring)),data=data.frame("ID"=seq(1,nrow(df.segm)),"state"=df.segm$state))
    dfSPLDF@data <- dplyr::mutate(dfSPLDF@data,colors=list_colors[state])

    # dfConnection <- SpatialLinesDataFrame(SpatialLines(ConnectionList,proj4string=CRS(UTMstring)),data=data.frame("ID"=seq(nrow(dfSeg)+1,2*nrow(dfSeg)-1),"behaviour"=c("transition"),"state"=4))
    dfConnection <- sp::SpatialLinesDataFrame(sp::SpatialLines(ConnectionList,proj4string=sp::CRS(UTMstring)),data=data.frame("ID"=seq(nrow(df.segm)+1,2*nrow(df.segm)-1),"behaviour"=c("transition"),"state"=4))

    dfSPPDF <- sp::SpatialPointsDataFrame(cbind(data$x,data$y),data = data.frame("ID"=seq(1,nrow(data)),"state"=data$state,"date"=data$indice),proj4string = sp::CRS(UTMstring))
    dfSPPDF@data <- dplyr::mutate(dfSPPDF@data,colors=list_colors[state])



    a <- leaflet::leaflet(width=width,height=height)  %>% leaflet::addPolylines(data=dfConnection,color = c("grey"),dashArray = '5,5')  %>%  leaflet::addPolylines(data=dfSPLDF,color =dfSPLDF@data$color)  %>% leaflet::addCircles(data=dfSPPDF,radius=1,color = dfSPPDF@data$colors,popup=paste("point ",as.character(dfSPPDF@data$date)," ; state ",as.character(dfSPPDF@data$state),sep=""))
    return(a)
  }
}
