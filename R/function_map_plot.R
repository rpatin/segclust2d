#'
#' \code{plot_segm} plot segmented movement data on a map.
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
#' @param ... additional arguments
#' @return a ggplot object
#'
#' @importFrom magrittr "%>%"
#' @export
#' @examples 
#' \dontrun{
#' #res.seg is a result of the segmentation-only algorithm : 
#' nseg = 10
#' outputs = res.seg$outputs[[paste(nseg, "segments")]]
#' map <- map_segm(data=res.seg$data,output=outputs)
#' #res.segclust is a result of the segmentation-clusturing algorithm : 
#' nseg = 10; ncluster = 3
#' outputs = res.segclust$outputs[[paste(ncluster,"class -",nseg, "segments")]]
#' map <- map_segm(data=res.seg$data,output=outputs)
#' }
#' 


map_segm <- function(data,output,interactive=FALSE,html=FALSE, scale=1,
                     UTMstring="+proj=longlat +datum=WGS84 +no_defs",
                     width=400,height=400,order=NULL,
                     pointsize = 1, linesize = 0.5, coord.names = c("x","y"),
                     ...){
  # print("test")
  df.segm <- dplyr::left_join(output[[1]],output[[2]],by="state")
  data$indice <- seq_len(nrow(data))
  
  data$state  <- df.segm[
    findInterval(data$indice, df.segm$begin,
                 rightmost.closed = FALSE, left.open = FALSE),
    ifelse(order,"state_ordered","state")
  ]
  
  data[,coord.names[1]] <- data[,coord.names[1]]/scale
  data[,coord.names[2]] <- data[,coord.names[2]]/scale
  if(!interactive){
    g <- ggplot2::ggplot(
      data,
      ggplot2::aes_string(x=coord.names[1],y=coord.names[2]))+
      ggplot2::geom_path(size = linesize)+
      ggplot2::geom_point(
        ggplot2::aes_string(col="factor(state)"),
        size = pointsize)
    return(g)
  } else {
    # coordinates(data) <- c("x","y")
    
    LinesList <- list()
    for(i in seq_len(nrow(df.segm))) {
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
    
    list_colors <-
      RColorBrewer::brewer.pal(
        n = length(unique(df.segm$state)),
        name='Set1'
        )
    
    dfSPLDF <- 
      sp::SpatialLinesDataFrame(
        sp::SpatialLines(LinesList,
                         proj4string = sp::CRS(UTMstring)),
        data=data.frame(
          "ID" = seq(1,nrow(df.segm)),
          "state"=df.segm$state
          )
        )
    dfSPLDF@data$colors <- list_colors[dfSPLDF@data$state]
    
    # dfConnection <-
    # SpatialLinesDataFrame(
    #   SpatialLines(ConnectionList,
    #                proj4string=CRS(UTMstring)),
    #   data=data.frame(
    #     "ID" = seq(nrow(dfSeg)+1,2*nrow(dfSeg)-1),
    #     "behaviour" = c("transition"),
    #     "state"=4)
    #   )
    dfConnection <-
      sp::SpatialLinesDataFrame(
        sp::SpatialLines(ConnectionList,
                         proj4string = sp::CRS(UTMstring)),
        data=data.frame(
          "ID" = seq(nrow(df.segm)+1,2*nrow(df.segm)-1),
          "behaviour" = c("transition"),
          "state"=4
          )
      )
    
    dfSPPDF <-
      sp::SpatialPointsDataFrame(
        cbind(data$x,data$y),
        data = data.frame(
          "ID"=seq(1,nrow(data)),
          "state"=data$state,
          "date"=data$indice
          ),
        proj4string = sp::CRS(UTMstring)
        )
    dfSPPDF@data$colors <- list_colors[dfSPPDF@data$state]
    
    
    
    a <- leaflet::leaflet(width=width,height=height)  %>% 
      leaflet::addPolylines(
        data=dfConnection,
        color = c("grey"),
        dashArray = '5,5')  %>% 
      leaflet::addPolylines(
        data=dfSPLDF,
        color =dfSPLDF@data$color)  %>% 
      leaflet::addCircles(
        data=dfSPPDF,radius=1,
        color = dfSPPDF@data$colors,
        popup=paste0("point ", 
                     as.character(dfSPPDF@data$date),
                     " ; state ",
                     as.character(dfSPPDF@data$state)))
    return(a)
  }
}





#' \code{segmap_list} create maps with a list of object of \code{segmentation}
#' class
#' @param x_list list of segmentation objects for different individuals or path
#' @param ncluster_list list of number of cluster to be selected for each
#'   individual. If empty, the function takes the default one
#' @param nseg_list list of number of segment to be selected for each
#'   individual. If empty, the function takes the default one
#' @param pointsize size of points
#' @param linesize size of lines
#' @param coord.names names of coordinates
#' @return a ggplot2 graph
#'
#' @export
segmap_list <- 
  function(x_list, 
           ncluster_list = NULL,
           nseg_list = NULL,
           pointsize = 1,
           linesize = 0.5,
           coord.names = c("x","y")){
  
  g <- ggplot2::ggplot()
  for(i in seq_along(x_list)){
    x <- x_list[[i]]
    if(is.null(ncluster_list)){
      ncluster <- x$ncluster.BIC
    } else {
      ncluster <- ncluster_list[i]
    }
    if(is.null(nseg_list)){
      nseg <- x$Kopt.BIC[ncluster]
    } else {
      nseg <- nseg_list[i]
    }
    outputs <- 
      x$outputs[[
        paste(ncluster,"class -",nseg, "segments")
        ]]
    data <- x$data
    df.segm <- dplyr::left_join(outputs[[1]],outputs[[2]],by="state")
    data$indice <- seq_len(nrow(data))
    data$state  <- 
      df.segm[
        findInterval(data$indice,
                     df.segm$begin,
                     rightmost.closed = FALSE,
                     left.open = FALSE),
        "state_ordered"]
    
    g <- g + 
      ggplot2::geom_path(
        data = data,
        ggplot2::aes_string(
          x=coord.names[1],
          y=coord.names[2]), 
        size = linesize)+
      ggplot2::geom_point(
        data = data,
        ggplot2::aes_string(
          x=coord.names[1],
          y=coord.names[2],
          col="factor(state)"),
        size = pointsize)
    
  }
  return(g)
}
