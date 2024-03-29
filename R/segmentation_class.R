#' segmentation class description
#'
#' @param x a \code{segmentation} object generated by
#'   \code{\link{segmentation}}
#' @param nseg number of segment chosen
#' @param ncluster number of classes chosen
#' @param ... additional arguments
#' @name segmentation-class
NULL

#' \code{print.segmentation} prints object of \code{segmentation} class
#' @param max.level argument to be passed to utils::str()
#' @rdname segmentation-class
#' @export


print.segmentation <- function(x, max.level = 1, ...) {
  utils::str(x, max.level = max.level)
}


#' \code{plot.segmentation} plot object of \code{segmentation} class -
#' wrapper for \code{\link{plot_segm}}
#' @param interactive whether plot are 
#' interactive (dygraphs/leaflet) or not (ggplot2)
#' @param xcol column for x axis. can be POSIXct
#' @param html whether htmltools::tagList should be applied on the returned
#'   object object for integrating in html pages
#' @rdname segmentation-class
#' @inheritParams plot_segm
#' @export
#'
#' @examples
#' \dontrun{
#' plot(res.segclust)
#' plot(res.segclust, nseg = 10, ncluster = 3)
#' }
#'
plot.segmentation <-
  function(x, nseg, 
           ncluster, interactive = FALSE,
           xcol = "indice", order, ...) {
    
    order <- argcheck_ordering(order, x$seg.type,
                               x$`Order variable`)
    
    if (x$seg.type == "segclust") {
      selected_param <- argcheck_segclust(ncluster, 
                                          nseg, 
                                          x$ncluster.BIC,
                                          x$Kopt.BIC)
      ncluster <- selected_param$ncluster
      nseg <- selected_param$nseg
      outputs <- x$outputs[[paste(ncluster, "class -", nseg, "segments")]]
      g <-
        plot_segm(data = x$data,
                  output = outputs, 
                  interactive = interactive, 
                  diag.var = x$`Diagnostic variables`, 
                  x_col = xcol, order = order, ...)
    } else if (x$seg.type == "segmentation") {
      nseg <- argcheck_segmentation(nseg, 
                                    x$Kopt.lavielle)
      
      outputs <- x$outputs[[paste(nseg, "segments")]]
      
      g <- plot_segm(
        data = x$data,
        output = outputs,
        interactive = interactive,
        diag.var = x$`Diagnostic variables`,
        x_col = xcol,
        order = order, ...
      )
    } else if (x$seg.type == "HMM" | 
               x$seg.type == "shiftfit" |
               x$seg.type == "depmixS4") {
      g <- plot_segm(data = x$data, output = x$outputs,
                     diag.var = x$`Diagnostic variables`,
                     x_col = xcol, order = order, ...)
    }
    return(g)
  }




#' \code{likelihood.segmentation} deprecated function for plotting likelihood
#' estimates of \code{segmentation} object. Now use \link{plot_likelihood}.
#' @rdname segmentation-class
#' @export
#'

likelihood.segmentation <- function(x, ...) {
  .Deprecated("plot_likelihood")
  plot_likelihood(x)
}

#' \code{plot_likelihood} plot likelihood estimates of a \code{segmentation}
#' object - works only for picard segmentation.
#' @rdname segmentation-class
#' @export
#' @examples
#' \dontrun{
#' plot_likelihood(res.seg)
#' }
#'
plot_likelihood <- function(x) {
  if (x$seg.type == "segclust") {
    likedat <- x$likelihood
    # nseg.bic <- x$Kopt.lavielle
    # tmpdf =  filter(li("nseg"=nseg.bic, "likelihood" =
    # x$likelihood$likelihood[which(x$likelihood$nseg == nseg.bic)])
    g <- ggplot2::ggplot(likedat[
      (likedat$ncluster != 0) & 
        is.finite(likedat$likelihood), ], 
      ggplot2::aes_string(x = "nseg",
                          y = "likelihood",
                          col = "factor(ncluster)")) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::xlab("Number of segments") +
      ggplot2::ylab("log-Likelihood") +
      # ggplot2::geom_point(data = tmpdf,ggplot2::aes(x=nseg,y=likelihood,col =
      # fact),size = 3)+
      ggplot2::scale_color_discrete(name = "Number of \nCluster")
  } else if (x$seg.type == "segmentation") {
    nseg.lav <- x$Kopt.lavielle
    tmpdf <-
      data.frame("nseg" = nseg.lav,
                 "likelihood" = x$likelihood$likelihood[
                   which(x$likelihood$nseg == nseg.lav)
                 ])
    nudgeY <- (max(x$likelihood$likelihood, na.rm = TRUE) -
                 min(x$likelihood$likelihood, na.rm = TRUE)) / 20
    nudgeX <- (max(x$likelihood$nseg, na.rm = TRUE) -
                 min(x$likelihood$nseg, na.rm = TRUE)) / 6
    g <-
      ggplot2::ggplot(
        x$likelihood,
        ggplot2::aes_string(x = "nseg", y = "likelihood")
      ) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::xlab("Number of segments") +
      ggplot2::ylab("log-Likelihood") +
      ggplot2::scale_color_discrete(name = "Number of \nCluster") +
      ggplot2::geom_point(data = tmpdf, 
                          ggplot2::aes_string(x = "nseg", y = "likelihood"),
                          size = 3) +
      ggplot2::geom_text(data = tmpdf, 
                         ggplot2::aes_string(x = "nseg", y = "likelihood"),
                         label = "Lavielle-selected optimum", 
                         nudge_x = nudgeX, nudge_y = -nudgeY, size = 3) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks())
  }
  return(g)
}

#' \code{get_likelihood} returns 
#' likelihood estimates of a \code{segmentation} object.
#' Deprecated, now use \link{logLik.segmentation}.
#' @rdname segmentation-class
#' @export

get_likelihood <- function(x) {
  .Deprecated("logLik.segmentation")
  return(x$likelihood)
}


#' \code{logLik.segmentation} returns 
#' log-likelihood estimates of a \code{segmentation} object
#' @rdname segmentation-class
#' @export
#' @examples
#' \dontrun{
#' logLik(res.seg)
#' }
#'
logLik.segmentation <- function(object, ...) {
  return(object$likelihood)
}




#' \code{plot_BIC} plot BIC estimates of a \code{segmentation} object
#' - works only for segclust algorithm.
#' @rdname segmentation-class
#' @export
#' @examples
#' \dontrun{
#' plot_BIC(res.segclust)
#' }
#'
plot_BIC <- function(x) {
  if (x$seg.type == "segclust") {
    likedat <- x$BIC
    ncluster.BIC <- x$ncluster.BIC
    Kopt.BIC <- x$Kopt.BIC[ncluster.BIC]
    ClusterOpt <-
      data.frame("ncluster" = ncluster.BIC,
                 "nseg" = Kopt.BIC, 
                 "BIC" = likedat[
                   (likedat$ncluster == ncluster.BIC) &
                     (likedat$nseg == Kopt.BIC), 
                 ]$BIC
      )
    
    nudgeY <-
      (max(likedat$BIC[is.finite(likedat$BIC)], na.rm = TRUE) - 
         min(likedat$BIC[is.finite(likedat$BIC)], na.rm = TRUE)) / 20
    nudgeX <- (max(likedat$nseg, na.rm = TRUE) -
                 min(likedat$nseg, na.rm = TRUE)) / 6
    
    ncluster.BIC <- 2:max(likedat$ncluster)
    Kopt.BIC <- x$Kopt.BIC[-1]
    SegOpt <- data.frame(ncluster = ncluster.BIC, nseg = Kopt.BIC)
    SegOpt <- dplyr::left_join(SegOpt, likedat, by = c("ncluster", "nseg"))
    g <- ggplot2::ggplot(
      likedat[is.finite(likedat$BIC), ], 
      ggplot2::aes_string(x = "nseg",
                          y = "BIC", 
                          col = "factor(ncluster)")) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::xlab("Number of segments") +
      ggplot2::ylab("BIC-based penalized log-Likelihood") +
      ggplot2::geom_point(data = SegOpt, shape = 15, size = 2) +
      ggplot2::geom_point(data = ClusterOpt, shape = 19, size = 3.5) +
      ggplot2::geom_text(data = ClusterOpt, size = 3,
                         label = "selected optimum",
                         nudge_x = -nudgeX, nudge_y = nudgeY) +
      ggplot2::scale_color_discrete(name = "Number of \nClusters") +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks())
  } else if (x$seg.type == "segmentation") {
    stop("no BIC estimates for segmentation only algorithm")
  }
  return(g)
}


#' \code{BIC} returns BIC-based penalized log-likelihood estimates of a
#' \code{segmentation} object when segmentation/clustering has been run.
#' @param object a segmentation-class object, created by segclust.
#' @rdname segmentation-class
#' @importFrom stats BIC
#' @export
#' @examples
#' \dontrun{
#' plot_BIC(res.segclust)
#' }
#'
BIC.segmentation <- function(object, ...) {
  return(object$BIC)
}

#' \code{stateplot} plot state distribution of a \code{segmentation} object
#' @rdname segmentation-class
#' @export
#' @examples
#' \dontrun{
#' stateplot(res.segclust)
#' stateplot(res.seg)
#' }
stateplot <- function(x, nseg, ncluster, order) {
  order <- argcheck_ordering(order, x$seg.type,
                             x$`Order variable`)
  
  if (x$seg.type == "segclust") {
    selected_param <- argcheck_segclust(ncluster, 
                                        nseg, 
                                        x$ncluster.BIC,
                                        x$Kopt.BIC)
    ncluster <- selected_param$ncluster
    nseg <- selected_param$nseg
    outputs <- x$outputs[[paste(ncluster, "class -", nseg, "segments")]]
    g <- plot_states(
      outputs,
      x$`Diagnostic variables`,
      order = order)
  } else if (x$seg.type == "segmentation") {
    nseg <- argcheck_segmentation(nseg, 
                                  x$Kopt.lavielle)
    
    outputs <- x$outputs[[paste(nseg, "segments")]]
    
    g <- plot_states(
      outputs, 
      x$`Diagnostic variables`, order = order)
  } else if (x$seg.type == "HMM" |
             x$seg.type == "shiftfit" |
             x$seg.type == "depmixS4") {
    g <- plot_states(x$outputs,
                     x$`Diagnostic variables`, 
                     order = order)
  }
  return(g)
}

#' \code{states} return data.frame with states statistics a \code{segmentation}
#' object
#' @rdname segmentation-class
#' @export
#' @examples
#' \dontrun{
#' states(res.segclust)
#' states(res.seg)
#' }
#'
states <- function(x, nseg, ncluster) {
  if (x$seg.type == "segclust") {
    selected_param <- argcheck_segclust(ncluster, 
                                        nseg, 
                                        x$ncluster.BIC,
                                        x$Kopt.BIC)
    ncluster <- selected_param$ncluster
    nseg <- selected_param$nseg
    outputs <- x$outputs[[paste(ncluster, "class -", nseg, "segments")]]
    return(outputs)
  } else if (x$seg.type == "segmentation") {
    nseg <- argcheck_segmentation(nseg, 
                                  x$Kopt.lavielle)
    
    outputs <- x$outputs[[paste(nseg, "segments")]]
    return(outputs)
    
  } else if (x$seg.type == "HMM" | 
             x$seg.type == "shiftfit" |
             x$seg.type == "depmixS4") {
    return(x$outputs$states)
  }
}

#' \code{segment} return data.frame with segment information of a
#' \code{segmentation} object
#' @rdname segmentation-class
#' @export
#' @examples
#' \dontrun{
#' segment(res.segclust)
#' segment(res.segclust, ncluster = 3, nseg = 30)
#' segment(res.seg)
#' segment(res.seg, nseg = 4)
#' }
segment <- function(x, nseg, ncluster) {
  if (x$seg.type == "segclust") {
    selected_param <- argcheck_segclust(ncluster, 
                                        nseg, 
                                        x$ncluster.BIC,
                                        x$Kopt.BIC)
    ncluster <- selected_param$ncluster
    nseg <- selected_param$nseg
    outputs <- x$outputs[[paste(ncluster, "class -", nseg, "segments")]]
    statesdf <- outputs$states
    segmentdf <- outputs$segments
    totdf <- dplyr::left_join(segmentdf, statesdf, by = "state")
    return(totdf)
  } else if (x$seg.type == "segmentation") {
    nseg <- argcheck_segmentation(nseg, 
                                  x$Kopt.lavielle)
    
    outputs <- x$outputs[[paste(nseg, "segments")]]
    statesdf <- outputs$states
    segmentdf <-outputs$segments
    totdf <- dplyr::left_join(segmentdf, statesdf, by = "state")
    return(totdf)
  } else if (x$seg.type == "HMM" |
             x$seg.type == "shiftfit" |
             x$seg.type == "depmixS4") {
    statesdf <- x$outputs$states
    segmentdf <- x$outputs$segments
    totdf <- dplyr::left_join(segmentdf, statesdf, by = "state")
    return(totdf)
  }
}

#' \code{augment.segmentation} return data.frame with original data and state
#' information of a \code{segmentation} object
#' @param colname_state column name for the added state column
#' @rdname segmentation-class
#' @export
#' @examples
#' \dontrun{
#' augment(res.segclust)
#' augment(res.segclust, ncluster = 3, nseg = 30)
#' augment(res.seg)
#' augment(res.seg, nseg = 4)
#' }
augment.segmentation <- 
  function(x, 
           nseg, ncluster,
           colname_state = "state", ...) {
    if (any(colnames(x$data) == colname_state)) {
      cli::cli_alert_danger(
        "{colname_state} already exists as column names of the \\
        original data.frame. Please provide argument colname_state \\
        with the name for the new state column.")
      stop("cannot erase already existing column name")
    }
    if (x$seg.type == "segclust") {
        selected_param <- argcheck_segclust(ncluster, 
                                            nseg, 
                                            x$ncluster.BIC,
                                            x$Kopt.BIC)
        ncluster <- selected_param$ncluster
        nseg <- selected_param$nseg
        outputs <- x$outputs[[paste(ncluster, "class -", nseg, "segments")]]
        statesdf <- outputs$states
        df.segm <- suppressMessages(segment(x, nseg = nseg))
      } else if (x$seg.type == "segmentation") {
        nseg <- argcheck_segmentation(nseg, 
                                      x$Kopt.lavielle)
        
        outputs <- x$outputs[[paste(nseg, "segments")]]
        statesdf <- outputs$states
        df.segm <- suppressMessages(segment(x, nseg = nseg))
      } else if (x$seg.type == "HMM" |
               x$seg.type == "shiftfit" |
               x$seg.type == "depmixS4") {
      statesdf <- states(x)
      df.segm <- suppressMessages(segment(x))
    }
    x$data$indice <- seq_len(nrow(x$data))
    data <- x$data
    data[, colname_state] <- 
      df.segm[findInterval(data$indice,
                           df.segm$begin,
                           rightmost.closed = FALSE, 
                           left.open = FALSE),
              "state"]
    totdf <- dplyr::left_join(data, statesdf, by = "state")
    
    return(totdf)
  }



#' \code{segmap} create maps with object of \code{segmentation} class
#'   (interpreting latitude/longitude)
#' @rdname segmentation-class
#' @inheritParams map_segm
#' @export
#'
#' @examples
#' \dontrun{
#' segmap(res.segclust, coord.names = c("x", "y"))
#' segmap(res.segclust, ncluster = 3, nseg = 30)
#' segmap(res.seg)
#' segmap(res.seg, nseg = 4)
#' }
segmap <- function(
  x, interactive = FALSE,
  nseg, ncluster, 
  html = FALSE,
  scale = 1, width = 400, height = 400, 
  order,
  pointsize = 1, linesize = 0.5, ...) {
  
  order <- argcheck_ordering(order, x$seg.type,
                             x$`Order variable`)
  
  if (x$seg.type == "segclust") {
    selected_param <- argcheck_segclust(ncluster, 
                                        nseg, 
                                        x$ncluster.BIC,
                                        x$Kopt.BIC)
    ncluster <- selected_param$ncluster
    nseg <- selected_param$nseg
    outputs <- x$outputs[[paste(ncluster, "class -", nseg, "segments")]]
  }  else if (x$seg.type == "segmentation") {
    nseg <- argcheck_segmentation(nseg, 
                                  x$Kopt.lavielle)
    
    outputs <- x$outputs[[paste(nseg, "segments")]]
  } else if (x$seg.type == "HMM" |
             x$seg.type == "shiftfit" | 
             x$seg.type == "depmixS4") {
    outputs <- x$outputs
  }
  map <- map_segm(data = x$data, output = outputs, 
                  interactive = interactive, html = html,
                  scale = scale, width = width, height = height, 
                  order = order,
                  pointsize = pointsize, linesize = linesize, ...)
  
  return(map)
}
