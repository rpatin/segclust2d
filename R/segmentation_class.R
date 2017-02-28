
#' Print object of \code{segmentation} class
#'
#' @rdname segmentation
#' @export


print.segmentation <- function(x,max.level = 1){
  str(x,max.level = max.level)
}


#' Plot object of \code{segmentation} class
#' Wrapper for \code{\link{plot_segm}}
#' @rdname segmentation
#' @export

#  x <- test.explo
plot.segmentation <- function(x,nseg=NULL, separate=T, interactive=F, xcol="indice", html = F) {
  if(x$type == 'picard'){
    if(is.null(nseg)) stop("nseg must be chosen for plotting picard segmentation")
    g <- plot_segm(data = x$data, output = x$outputs[[paste(nseg, "segments")]], separate = T, interactive=interactive, diag.var = x$`Diagnostic variables`,x_col = xcol, html = html)
    return(g)
  } else if (x$type == 'HMM'){
    g <- plot_segm(data = x$data, output = x$outputs, separate = T, interactive=interactive, diag.var = x$`Diagnostic variables`,x_col = xcol, html = html)
    return(g)
  }
}




#' Plot likelihood estimates of a \code{segmentation} object
#' Works only for picard segmentation.
#' @rdname segmentation
#' @export

likelihood.segmentation <- function(x) {
  if(x$type != 'picard') stop("likelihood only pertinent for picard segmentation")
  g <- ggplot2::ggplot(x$likelihood,ggplot2::aes(x=nseg,y=likelihood))+ggplot2::geom_point()+ggplot2::geom_line()+ggplot2::xlab("Number of segments")+ggplot2::ylab("log-Likelihood")
  return(g)
}


#' Plot state distribution of a \code{segmentation} object
#' @rdname segmentation
#' @export

stateplot <- function(x,nseg = NULL){
  if(x$type == 'picard'){
    if(is.null(nseg)) stop("nseg must be chosen for plotting picard states statistics")
    g <- plot_states(x$outputs[[paste(nseg, "segments")]],x$`Diagnostic variables`)
    return(g)
  } else if (x$type == 'HMM'){
    g <- plot_states(x$outputs,x$`Diagnostic variables`)
    return(g)
  }

}

#' Return data.frame with states statistics a \code{segmentation} object
#' @rdname segmentation
#' @export

states <- function(x,nseg = NULL){
  if(x$type == 'picard'){
    if(is.null(nseg)) stop("nseg must be chosen for getting states statistics")
    return(x$outputs[[paste(nseg, "segments")]]$states)
  } else if (x$type == 'HMM'){
    return(x$outputs$states)
  }
}

#' Return data.frame with segment information of a \code{segmentation} object
#' @rdname segmentation
#' @export

segment <- function(x,nseg = NULL){
  if(x$type == 'picard'){
    if(is.null(nseg)) stop("nseg must be chosen for getting segment statistics")
    return(x$outputs[[paste(nseg, "segments")]]$segments)
  } else if (x$type == 'HMM'){
    return(x$outputs$segments)
  }
}

#' Return data.frame with original data and state information of a \code{segmentation} object
#' @rdname segmentation
#' @export

augment.segmentation<- function(x,nseg = NULL,colname_state = "state"){
  if(any(colnames(x$data) == colname_state)) stop(paste(colname_state,"already exists as column names of the data.frame. Cannot erase"))

  if(x$type == 'picard'){
    if(is.null(nseg)) stop("nseg must be chosen for getting segment statistics")
    df.segm  <- segment(x,nseg=nseg)
    x$data$indice <- 1:nrow(x$data)
    evalstr  <- paste("data <- dplyr::mutate(x$data,",colname_state,"= df.segm[findInterval(indice,df.segm$begin,rightmost.closed = F,left.open = F),\"state\"])",sep="")
    eval(parse(text=evalstr))
    return(data)
    } else if (x$type == 'HMM'){
      cluster <- moveHMM::viterbi(x$model.hmm)
      data[,colname_state] <- cluster
  }
  return(data)
}



#' Create maps with object of \code{segmentation} class
#'
#' @rdname segmentation
#' @export

segmap <-  function(x,interactive=F,nseg=NULL,x_col="expectTime",html=F,scale=100,UTMstring="+proj=utm +zone=35 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0",width=400,height=400){
  if(x$type == 'picard'){
    if(is.null(nseg)) stop("nseg must be chosen for getting segment statistics")
    map <- map_segm(data=x$data,output=x$outputs[[paste(nseg, "segments")]],interactive = interactive, x_col = x_col, html = html, scale=scale, UTMstring = UTMstring,width=width,height=height)
  } else if (x$type == 'HMM'){
    map <- map_segm(data = x$data, output = x$outputs, interactive = interactive, x_col = x_col, html = html, scale=scale, UTMstring = UTMstring,width=width,height=height)
  }
  return(map)
}

