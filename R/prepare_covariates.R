#' Add several covariates to movement observations
#' \code{add_covariates} add several covariates to a data frame with movement
#' information. It adds : distance between location, spatial angle, speed,
#' smoothed speed, persistence and rotation velocity (calculated with spatial
#' angle).
#' @param x data.frame with locations
#' @param coord.names names of coordinates column in \code{x}
#' @param timecol names of POSIXct time column
#' @param smoothed whether speed are smoothed or not
#' @param units units for time calculation. Default "hour"
#' @param radius for spatial angle calculations
#' @return data.frame increased with
#'
#' @examples
#' calc_dist(df,coord.names = c("x","y"), smoothed = T)
#' @export
#' @author Remi Patin


add_covariates <- function(x, coord.names = c("x","y"), smoothed = F, timecol = "dateTime", units = "hour", radius = NULL){
  if(any(is.na(x[,timecol]))) stop("time should not contain NA")
  if(any(is.na(x[,coord.names[1]]))) stop("x should not contain NA")
  if(any(is.na(x[,coord.names[2]]))) stop("y should not contain NA")

  dist <- calc_dist(x, coord.names = coord.names, smoothed = F)
  dist_smoothed <- zoo::rollapply(dist, 2, mean, by.column = FALSE, fill = NA, align = "right")
  n <- nrow(x)
  tmptime =  as.numeric(difftime(x[2:n,timecol],x[1:(n-1),timecol], units = units))
  speed = dist/c(tmptime,NA)
  speed_smoothed <- zoo::rollapply(speed, 2, mean, by.column = FALSE, fill = NA, align = "right")
  ang_spa <- spatial_angle(x, coord.names = coord.names, radius = radius)
  vit_p <- speed*cos(ang_spa)
  vit_r <- speed*sin(ang_spa)

  x$dist <- dist
  x$dist_smoothed <- dist_smoothed
  x$speed <- speed
  x$speed_smoothed <- speed_smoothed
  x$spatial_angle <- ang_spa
  x$vit_p <- vit_p
  x$vit_r <- vit_r
  return(x)
}

#' Calculate distance between locations
#'
#' \code{calc_dist} calculate distance between locations, taking a dataframe as input. Distance can also be smoothed over the two steps before and after the each point.
#' @param x data.frame with locations
#' @param coord.names names of coordinates column in \code{x}
#' @param smoothed whether distance are smoothed or not
#' @return vector of distance
#'
#' @examples
#' calc_dist(df,coord.names = c("x","y"), smoothed = T)
#' @export
#' @author Remi Patin


calc_dist <- function(x, coord.names = c("x","y"), smoothed = F){
  tmp = zoo::rollapply(cbind(x[,coord.names[1]], x[,coord.names[2]]), 2, dist, by.column = FALSE, fill = NA)
  if( smoothed ){
    tmp <- zoo::rollapply(tmp, 2, mean, by.column = FALSE, fill = NA, align = "right")
  }
  return(tmp)
}

#' Calculate speed along a trajectory
#'
#' \code{calc_dist} calculate speed between locations, taking a dataframe as input. Speed can also be smoothed over the two steps before and after the each point.
#' @param x data.frame with locations
#' @param coord.names names of coordinates column in \code{x}
#' @param timecol names of POSIXct time column
#' @param smoothed whether speed are smoothed or not
#' @param units units for time calculation. Default "hour"
#' @return vector of distance
#'
#' @examples
#' calc_speed(df,coord.names = c("x","y"), timecol = "dateTime", smoothed = T)
#' @export
#' @author Remi Patin


calc_speed <- function(x, coord.names = c("x","y"), timecol = "dateTime", smoothed = F, units = "hour"){
  tmpdist = calc_dist(x, coord.names = coord.names)
  n <- nrow(x)
  tmptime =  as.numeric(difftime(x[2:n,timecol],x[1:(n-1),timecol], units = units))
  tmpspeed = tmpdist/c(tmptime,NA)
  if(smoothed){
    tmpspeed <- zoo::rollapply(tmpspeed, 2, mean, by.column = FALSE, fill = NA, align = "right")
  }
  return(tmpspeed)
}


#' Calculate spatial angle along a trajectory
#'
#' \code{spatial_angle} calculate spatial angle between locations, taking a dataframe as input. Spatial angle is considered as the angle between the focus point, the first location entering a given circle and the last location inside.
#' @param x data.frame with locations
#' @param coord.names names of coordinates column in \code{x}
#' @param radius for angle calculation. Default is median of step length.
#' @return vector of spatial angle.
#'
#' @examples
#' calc_speed(df,coord.names = c("x","y"), timecol = "dateTime", smoothed = T)
#' @export
#' @author Remi Patin, Simon Benhamou.


spatial_angle <- function(df, coord.names = c("x","y"), radius = NULL){
  tmpdist = calc_dist(df, coord.names = coord.names)
  x <- df[,coord.names[1]]
  y <- df[,coord.names[2]]
  if(is.null(radius)) radius <- median(tmpdist,na.rm=T)
  radius2 <- radius^2
  ri2 <- 0.998*radius2
  re2 <- 1.002*radius2
  n <- nrow(df)
  angle_spa <- c(NA)
  arg1 <- c(NA)
  arg2 <- c(NA)
  for(i in 2:(n-1)){
    # message(i)
    flag <- T
    j <- i; d2 <- 0; dxs <- 0; dys <- 0
    # forward. Find first point outside of circle
    while (d2 <= ri2 & j < n){
      j <- j+1
      dxp <- dxs
      dyp <- dys
      dxs <- x[i] - x[j]
      dys <- y[i] - y[j]
      d2 <- dxs^2+dys^2
    }
    if (d2>ri2){
      # interpolation if needed (no interpolation if points fall between ri and re)
      if (d2<re2){
        xinter <- x[j]
        yinter <- y[j]
      } else {
        c <- (x[j]-x[j-1])/tmpdist[j-1] # cosinus
        s <- (y[j]-y[j-1])/tmpdist[j-1] # sinus
        rd <- ( dxp*c + dyp*s+ sqrt( radius2-(dyp*c-dxp*s)^2 ) )/tmpdist[j-1] # attention au changement de repere.
        ard <- 1-rd
        xinter <- x[j-1]*ard+x[j]*rd
        yinter <- y[j-1]*ard+y[j]*rd
      }
      cs <- xinter-x[i]
      ss <- yinter-y[i]
    } else {
      flag = F
    }
    j <- i; d2 <- 0; dxs <- 0; dys <- 0
    # backward, find first point outside circle
    while (d2 <= ri2 & j > 1){
      j <- j-1
      dxp <- dxs
      dyp <- dys
      dxs <- x[i] - x[j]
      dys <- y[i] - y[j]
      d2 <- dxs^2+dys^2
    }
    if (d2>ri2){
      #interpolation
      if (d2<re2){
        xinter <- x[j]
        yinter <- y[j]
      } else {
        c <- (x[j]-x[j+1])/tmpdist[j] # cosinus
        s <- (y[j]-y[j+1])/tmpdist[j] # sinus
        rd <- (dxp*c+dyp*s+sqrt(radius2-(dyp*c-dxp*s)^2))/tmpdist[j] # attention au changement de repere.
        ard <- 1-rd
        xinter <- x[j+1]*ard+x[j]*rd
        yinter <- y[j+1]*ard+y[j]*rd
      }
      cp <- x[i]-xinter
      sp <- y[i]-yinter
    } else {
      flag = F
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
