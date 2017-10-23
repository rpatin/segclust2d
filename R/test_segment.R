#' #' test_stationarity checks the stationarity within each segment
#' #'
#' #' \code{test_stationarity}
#' #' @export
#'
#' test_stationarity <- function(dat,outputs,Kmax){
#'   segment_stationarity <- list()
#'   for (nseg in 1:Kmax){
#'     tmp = data.frame()
#'     out <- outputs[[paste(nseg,'segments')]]
#'     for (seg in 1:nseg){
#'       start <- out$segments[seg,'begin']
#'       stop <- out$segments[seg,'end']
#'       serie_x <- dat[1,start:stop]
#'       serie_y <- dat[2,start:stop]
#'       chunks <- chunk2(serie_x, 3)
#'       tmp_df <- data.frame(serie_x,serie_y,chunks)
#'       tmp.lm.x <- lm(data = tmp_df, formula = serie_x ~ factor(chunks))
#'       tmp.lm.y <- lm(data = tmp_df, formula = serie_y ~ factor(chunks))
#'       tmp.aov.x <- summary(aov(tmp.lm.x))[[1]]$`Pr(>F)`[1]
#'       tmp.aov.y <- summary(aov(tmp.lm.y))[[1]]$`Pr(>F)`[1]
#'       # t1 <- Box.test(serie_x,lag = min(length(serie_x),20),type = "Ljung-Box")$p.value <= 0.05
#'       # t2 <- Box.test(serie_y,lag = min(length(serie_x),20),type = "Ljung-Box")$p.value <= 0.05
#'       t1 <- tmp.aov.x > 0.05
#'       t2 <- tmp.aov.y > 0.05
#'       tmp = rbind(tmp, data.frame('seg' = seg,
#'                                   'stat1' = ifelse(t1,"stationary","non-stationary"),
#'                                   'stat2' = ifelse(t2,"stationary","non-stationary"),
#'                                   'stat_tot' = ifelse(t1 & t2,"stationary","non-stationary")))
#'     }
#'     segment_stationarity[[paste(nseg,'segments')]] <- tmp
#'   }
#'   return(segment_stationarity)
#' }
#'
#' #' test_var checks the variance between each pair of segment
#' #'
#' #' \code{test_stationarity}
#' #' @export
#'
#' test_var <- function(dat,outputs,Kmax){
#'   segment_var <- list()
#'   for (nseg in 2:Kmax){
#'     tmpx = matrix(NA,nrow = nseg, ncol = nseg)
#'     tmpy = matrix(NA,nrow = nseg, ncol = nseg)
#'     tmpboth = matrix(NA,nrow = nseg, ncol = nseg)
#'     out <- outputs[[paste(nseg,'segments')]]
#'     for (seg1 in 1:(nseg-1)){
#'       for (seg2 in (seg1+1):(nseg)){
#'         start1 <- out$segments[seg1,'begin']
#'         stop1 <- out$segments[seg1,'end']
#'         serie_x1 <- dat[1,start1:stop1]
#'         serie_y1 <- dat[2,start1:stop1]
#'
#'         start2 <- out$segments[seg2,'begin']
#'         stop2 <- out$segments[seg2,'end']
#'         serie_x2 <- dat[1,start2:stop2]
#'         serie_y2 <- dat[2,start2:stop2]
#'
#'         tx <- var.test(serie_x1, serie_x2, ratio = 1,
#'                        alternative = c("two.sided"))$p.value <= 0.05
#'         ty <- var.test(serie_y1, serie_y2, ratio = 1,
#'                        alternative = c("two.sided"))$p.value <= 0.05
#'         tboth <- ifelse(tx | ty,"Different Variances","Same Variance")
#'         tx = ifelse(tx,"Different Variances","Same Variance")
#'         ty = ifelse(ty,"Different Variances","Same Variance")
#'         tmpx[seg1,seg2] <- tx
#'         tmpy[seg1,seg2] <- ty
#'         tmpboth[seg1,seg2] <- tboth
#'       }
#'     }
#'     segment_var[[paste(nseg,'segments')]] <- list("x" = tmpx,
#'                                                   "y" = tmpy,
#'                                                   "both"= tmpboth)
#'   }
#'   return(segment_var)
#' }
#'
#'
#' #' test_mean checks the mean between each pair of segment
#' #'
#' #' \code{test_mean}
#' #' @export
#'
#' test_mean <- function(dat,outputs,Kmax){
#'   segment_mean <- list()
#'   for (nseg in 2:Kmax){
#'     tmpx = matrix(NA,nrow = nseg, ncol = nseg)
#'     tmpy = matrix(NA,nrow = nseg, ncol = nseg)
#'     tmpboth = matrix(NA,nrow = nseg, ncol = nseg)
#'     out <- outputs[[paste(nseg,'segments')]]
#'     for (seg1 in 1:(nseg-1)){
#'       for (seg2 in (seg1+1):(nseg)){
#'         start1 <- out$segments[seg1,'begin']
#'         stop1 <- out$segments[seg1,'end']
#'         serie_x1 <- dat[1,start1:stop1]
#'         serie_y1 <- dat[2,start1:stop1]
#'
#'         start2 <- out$segments[seg2,'begin']
#'         stop2 <- out$segments[seg2,'end']
#'         serie_x2 <- dat[1,start2:stop2]
#'         serie_y2 <- dat[2,start2:stop2]
#'
#'         tx <- t.test(serie_x1, serie_x2, var.equal = F,
#'                      alternative = c("two.sided"))$p.value <= 0.05
#'         ty <-  t.test(serie_y1, serie_y2, var.equal = F,
#'                       alternative = c("two.sided"))$p.value <= 0.05
#'         tboth <- ifelse(tx | ty,"Different Means","Same Mean")
#'         tx = ifelse(tx,"Different Means","Same Mean")
#'         ty = ifelse(ty,"Different Means","Same Mean")
#'         tmpx[seg1,seg2] <- tx
#'         tmpy[seg1,seg2] <- ty
#'         tmpboth[seg1,seg2] <- tboth
#'       }
#'     }
#'     segment_mean[[paste(nseg,'segments')]] <- list("x" = tmpx,
#'                                                    "y" = tmpy,
#'                                                   "both"= tmpboth)
#'   }
#'   return(segment_mean)
#' }
#'
#' #' chunk2 chunks a serie into several groups
#' #'
#' #' \code{chunk2}
#'
#'
#' chunk2 <- function(x,n){
#'   return(as.numeric(Hmisc::cut2(1:length(x),g=n)))
#' }
