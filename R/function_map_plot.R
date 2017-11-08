#'
#' \code{plot_segm} plot segmented traject on a map.
#' @param data the data.frame with the different variable
#' @param output outputs of the segmentation  or segclust algorithm for one
#'   number of segment
#' @param interactive should graph be interactive with leaflet ?
#' @param html should the graph be incorporated in a markdown file through
#'   htmltools::tagList()
#' @param order should cluster be ordered
#' @param pointsize size of points
#' @param height height
#' @param width width
#' @param linesize size of lines
#' @param scale for dividing coordinates to have compatibility with leaflet
#' @param UTMstring projection of the coordinates
#' @param coord.names names of coordinates
#' @return a graph
#'
#' @examples
#' \dontrun{stat_segm(data, diag.var=c("dist","angle"),
#' order.var='dist',type='hmm',hmm.model=mod1.hmm)}
#' @importFrom magrittr "%>%"
#' @export

map_segm <- function(data,output,interactive=F,html=F, scale=1,
                     UTMstring="+proj=longlat +datum=WGS84 +no_defs",
                     width=400,height=400,order=NULL,pointsize = 1, linesize = 0.5, coord.names = c("x","y")){
  # print("test")
  df.segm <- dplyr::left_join(output[[1]],output[[2]],by="state")
  data$indice <- 1:nrow(data)

  data$state = df.segm[findInterval(data$indice,df.segm$begin,rightmost.closed = F,left.open = F),ifelse(order,"state_ordered","state")]

  data$x <- data$x/scale
  data$y <- data$y/scale
  if(!interactive){
    g <- ggplot2::ggplot(data,ggplot2::aes_string(x=coord.names[1],y=coord.names[2]))+
      ggplot2::geom_path(size = linesize)+
      ggplot2::geom_point(ggplot2::aes_string(col="factor(state)"),size = pointsize)
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
    dfSPLDF@data$colors <- list_colors[dfSPLDF@data$state]

    # dfConnection <- SpatialLinesDataFrame(SpatialLines(ConnectionList,proj4string=CRS(UTMstring)),data=data.frame("ID"=seq(nrow(dfSeg)+1,2*nrow(dfSeg)-1),"behaviour"=c("transition"),"state"=4))
    dfConnection <- sp::SpatialLinesDataFrame(sp::SpatialLines(ConnectionList,proj4string=sp::CRS(UTMstring)),data=data.frame("ID"=seq(nrow(df.segm)+1,2*nrow(df.segm)-1),"behaviour"=c("transition"),"state"=4))

    dfSPPDF <- sp::SpatialPointsDataFrame(cbind(data$x,data$y),data = data.frame("ID"=seq(1,nrow(data)),"state"=data$state,"date"=data$indice),proj4string = sp::CRS(UTMstring))
    dfSPPDF@data$colors=list_colors[dfSPPDF@data$state]



    a <- leaflet::leaflet(width=width,height=height)  %>% leaflet::addPolylines(data=dfConnection,color = c("grey"),dashArray = '5,5')  %>%  leaflet::addPolylines(data=dfSPLDF,color =dfSPLDF@data$color)  %>% leaflet::addCircles(data=dfSPPDF,radius=1,color = dfSPPDF@data$colors,popup=paste("point ",as.character(dfSPPDF@data$date)," ; state ",as.character(dfSPPDF@data$state),sep=""))
    return(a)
  }
}
