#' Covariate Calculations
#'
#' Add several covariates to movement observations
#' \code{add_covariates} add several covariates to a data frame with movement
#' information. It adds : distance between location, spatial angle, speed,
#' smoothed speed, persistence and rotation velocity (calculated with spatial
#' angle).
#'
#' @param x movement data
#' @param ... additional arguments
#' @export
#'

add_covariates <- function(x, ...) {
  UseMethod("add_covariates", x)
}


#' add_covariates method for Move object
#'
#'
#' @rdname add_covariates
#' @examples
#' \dontrun{add_covariates(move_object, coord.names = c("x","y"), smoothed = T)}
#' @export

add_covariates.Move <- function(x, coord.names = c("x","y"), ...){
  if(!requireNamespace("move", quietly = TRUE))
    stop("move package required 
         for calling segclust on a Move object.")
  if(!requireNamespace("sp", quietly = TRUE))
    stop("sp package required for calling 
         segclust (home-range) on a Move object.")
  x.dat <- x@data
  if(! all(coord.names %in% colnames(x.dat))){
    dat <- data.frame(sp::coordinates(x))
    colnames(dat) <- coord.names
    x.dat <- cbind(dat,x.dat)
  }
  add_covariates.data.frame(x.dat, coord.names = coord.names, ...)
}

#' add_covariates method for ltraj object
#'
#' @export
#' @rdname add_covariates

add_covariates.ltraj <- function(x, coord.names = c("x","y"), ...){
  if(!requireNamespace("adehabitatLT", quietly = TRUE)){
    stop("adehabitatLT package required 
         for calling segclust on a ltraj object.")
  }
  
  
  x.dat <- adehabitatLT::infolocs(x)[[1]]
  
  if(! all(coord.names %in% colnames(x.dat))){
    tmp <- x[[1]]
    if(any(is.na(tmp$x))) stop("Please filter NA from ltraj object")
    dat <- data.frame(cbind(tmp$x,tmp$y))
    colnames(dat) <- coord.names
    x.dat <- cbind(dat,x.dat)
  }
  add_covariates.data.frame(x.dat, coord.names = coord.names, ...)
}


#' add_covariates method for data.frame
#'
#' @param coord.names names of coordinates column in \code{x}
#' @param timecol names of POSIXct time column
#' @param smoothed whether speed are smoothed or not
#' @param units units for time calculation. Default "hour"
#' @param radius for spatial angle calculations
#' @return data.frame with additional covariates
#' @rdname add_covariates
#' @export
#' @examples 
#' \dontrun{
#' data(simulmode)
#' simple_data <- simulmode[,c("dateTime","x","y")]
#' full_data   <- add_covariates(simple_data, coord.names = c("x","y"),
#'  timecol = "dateTime",smoothed = TRUE, units ="min")
#' }


add_covariates.data.frame <-
  function(x, coord.names = c("x","y"), 
           smoothed = FALSE, timecol = "dateTime",
           units = "hour", radius = NULL, ...){
    if(any(is.na(x[,timecol]))) stop("time should not contain NA")
    if(any(is.na(x[,coord.names[1]]))) stop("x should not contain NA")
    if(any(is.na(x[,coord.names[2]]))) stop("y should not contain NA")
    if(missing(units)) {
      message("Using hours as default time unit. 
            You can change this with argument units.")
    }
    x_dist <- calc_dist(x, coord.names = coord.names, smoothed = FALSE)
    x_dist_smoothed <-
      zoo::rollapply(x_dist, 2, mean,
                     by.column = FALSE, 
                     fill = NA, align = "right")
    n <- nrow(x)

    tmptime  <- as.numeric(
      difftime(x[2:n,timecol, drop = TRUE],
               x[1:(n-1),timecol, drop = TRUE], 
               units = units)
    )
    x_speed <- x_dist/c(tmptime,NA)
    x_speed_smoothed <-
      zoo::rollapply(x_speed, 2, mean,
                     by.column = FALSE,
                     fill = NA, align = "right")
    x_ang_spa <- spatial_angle(x, coord.names = coord.names, radius = radius)
    x_vit_p <- x_speed*cos(x_ang_spa)
    x_vit_r <- x_speed*sin(x_ang_spa)
    x$angular_speed <- angular_speed(x, coord.names = coord.names)/c(tmptime,NA)
    x$angle <- angular_speed(x, coord.names = coord.names)
    x$dist <- x_dist
    x$dist_smoothed <- x_dist_smoothed
    x$speed <- x_speed
    x$speed_smoothed <- x_speed_smoothed
    x$spatial_angle <- x_ang_spa
    x$vit_p <- x_vit_p
    x$vit_r <- x_vit_r
    return(x)
  }

#' Calculate distance between locations
#'
#' \code{calc_dist} calculate distance between locations, taking a dataframe as
#' input. Distance can also be smoothed over the two steps before and after the
#' each point.
#' @param x data.frame with locations
#' @param coord.names names of coordinates column in \code{x}
#' @param smoothed whether distance are smoothed or not
#' @return vector of distance
#'
#' @examples
#' \dontrun{calc_dist(df,coord.names = c("x","y"), smoothed = T)}
#' @export
#' @author Remi Patin


calc_dist <- function(x, coord.names = c("x","y"), smoothed = FALSE){
  tmp <- 
    zoo::rollapply(
      cbind(x[,coord.names[1], drop = TRUE],
            x[,coord.names[2], drop = TRUE]), 2,
      function(x){
        x.tmp <- stats::dist(x)
        return(as.numeric(x.tmp))
      },
      by.column = FALSE,
      fill = NA
    )
  # tmp <- as.numeric(tmp)
  if( smoothed ){
    tmp <-
      zoo::rollapply(tmp, 2, mean, 
                     by.column = FALSE, 
                     fill = NA, align = "right")
  }
  return(tmp)
}

#' Calculate speed along a path
#'
#' \code{calc_dist} calculate speed between locations, taking a dataframe as
#' input. Speed can also be smoothed over the two steps before and after the
#' each point.
#' @param x data.frame with locations
#' @param coord.names names of coordinates column in \code{x}
#' @param timecol names of POSIXct time column
#' @param smoothed whether speed are smoothed or not
#' @param units units for time calculation. Default "hour"
#' @return vector of distance
#'
#' @examples
#' \dontrun{calc_speed(df,coord.names = c("x","y"), timecol = "dateTime",
#' smoothed = T)}
#' @export
#' @author Remi Patin


calc_speed <- 
  function(x,
           coord.names = c("x","y"),
           timecol = "dateTime",
           smoothed = FALSE, units = "hour"
  ){
    tmpdist  <- calc_dist(x, coord.names = coord.names)
    n <- nrow(x)
    tmptime  <- as.numeric(
      difftime(x[2:n,timecol, drop = TRUE],
               x[1:(n - 1),timecol, drop = TRUE],
               units = units)
    )
    tmpspeed  <- tmpdist/c(tmptime,NA)
    if (smoothed) {
      tmpspeed <-
        zoo::rollapply(tmpspeed, 2, mean,
                       by.column = FALSE, 
                       fill = NA, align = "right")
    }
    return(tmpspeed)
  }


#' Calculate spatial angle along a path
#'
#' \code{spatial_angle} calculate spatial angle between locations, taking a
#' dataframe as input. Spatial angle is considered as the angle between the
#' focus point, the first location entering a given circle and the last location
#' inside.
#' @param x data.frame with locations
#' @param coord.names names of coordinates column in \code{x}
#' @param radius for angle calculation. Default is median of step length.
#'
#' @return vector of spatial angle.
#'
#' @export
#' @author Remi Patin, Simon Benhamou.
#' @examples 
#' \dontrun{
#' data(simulmode)
#' spatial_angle(simulmode)
#' }



spatial_angle <- function(x, coord.names = c("x","y"), radius = NULL){
  tmpdist  <- calc_dist(x, coord.names = coord.names)
  xx <- x[ , coord.names[1], drop = TRUE]
  yy <- x[ , coord.names[2], drop = TRUE]
  if(is.null(radius)) radius <- stats::median(tmpdist,na.rm=TRUE)
  radius2 <- radius^2
  ri2 <- 0.998*radius2
  re2 <- 1.002*radius2
  n <- nrow(x)
  angle_spa <- c(NA)
  arg1 <- c(NA)
  arg2 <- c(NA)
  for(i in 2:(n-1)){
    # message(i)
    flag <- TRUE
    j <- i; d2 <- 0; dxs <- 0; dys <- 0
    # forward. Find first point outside of circle
    while (d2 <= ri2 & j < n){
      j <- j+1
      dxp <- dxs
      dyp <- dys
      dxs <- xx[i] - xx[j]
      dys <- yy[i] - yy[j]
      d2 <- dxs^2+dys^2
    }
    if (d2>ri2){
      # interpolation if needed
      # (no interpolation if points fall between ri and re)
      if (d2<re2){
        xinter <- xx[j]
        yinter <- yy[j]
      } else {
        c <- (xx[j]-xx[j-1])/tmpdist[j-1] # cosinus
        s <- (yy[j]-yy[j-1])/tmpdist[j-1] # sinus
        rd <- ( dxp*c + dyp*s+ 
                  sqrt( radius2-(dyp*c-dxp*s)^2 ) 
        )/tmpdist[j-1] # beware the change in coordinate system 
        ard <- 1-rd
        xinter <- xx[j-1]*ard+xx[j]*rd
        yinter <- yy[j-1]*ard+yy[j]*rd
      }
      cs <- xinter-xx[i]
      ss <- yinter-yy[i]
    } else {
      flag <- FALSE
    }
    j <- i; d2 <- 0; dxs <- 0; dys <- 0
    # backward, find first point outside circle
    while (d2 <= ri2 & j > 1){
      j <- j-1
      dxp <- dxs
      dyp <- dys
      dxs <- xx[i] - xx[j]
      dys <- yy[i] - yy[j]
      d2 <- dxs^2+dys^2
    }
    if (d2>ri2){
      #interpolation
      if (d2<re2){
        xinter <- xx[j]
        yinter <- yy[j]
      } else {
        c <- (xx[j]-xx[j+1])/tmpdist[j] # cosinus
        s <- (yy[j]-yy[j+1])/tmpdist[j] # sinus
        rd <- (dxp*c+dyp*s+
                 sqrt(radius2-(dyp*c-dxp*s)^2)
        )/tmpdist[j] # attention au changement de repere.
        ard <- 1-rd
        xinter <- xx[j+1]*ard+xx[j]*rd
        yinter <- yy[j+1]*ard+yy[j]*rd
      }
      cp <- xx[i]-xinter
      sp <- yy[i]-yinter
    } else {
      flag <- FALSE
    }
    tmp <- NA
    if(flag){
      tmp <- atan2(ss*cp-cs*sp,cp*cs+sp*ss)
    }
    angle_spa <- c(angle_spa,tmp)
    
  }
  angle_spa <- c(angle_spa,NA)
  return(angle_spa)
}

#' Calculate angular speed along a path
#'
#' \code{angular_speed} calculate turning angle between locations, taking a
#' dataframe as input.
#' @param x data.frame with locations
#' @param coord.names names of coordinates column in \code{x}
#' @return vector of turning angle.
#'
#' @export
#' @author Remi Patin, Simon Benhamou.

angular_speed <- function(x, coord.names = c("x","y")){
  
  xx <- diff(x[,coord.names[1], drop = TRUE])
  yy <- diff(x[,coord.names[2], drop = TRUE])
  b<-sign(xx)
  b[b==0] <- 1  #corrects for the fact that sign(0) == 0
  bearings <- b*(yy<0)*pi+atan(xx/yy)
  current_angular_speed <- diff(bearings)
  current_angular_speed <- ifelse(current_angular_speed < -pi,
                                  current_angular_speed + 2*pi,
                                  ifelse(current_angular_speed > pi,
                                         current_angular_speed -2*pi, 
                                         current_angular_speed))
  angular_speed <- c(NA,current_angular_speed,NA)
}
