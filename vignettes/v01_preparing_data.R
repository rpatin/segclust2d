## ----option chunk, echo = FALSE-----------------------------------------------
options(Encoding="UTF-8")
knitr::opts_chunk$set(
  fig.width = 8,
  fig.height = 5,
  collapse = TRUE,
  comment = "#>"
)

## ----library and data, fig.show='hold'----------------------------------------
library(segclust2d)
data(simulshift)
data(simulmode)
simulmode$abs_spatial_angle <- abs(simulmode$spatial_angle)
simulmode <- simulmode[!is.na(simulmode$abs_spatial_angle), ]

## ----ex df, eval = FALSE------------------------------------------------------
#  segmentation(x_data.frame, lmin = 5, Kmax = 25)

## ----ex Move, eval = FALSE----------------------------------------------------
#  segmentation(x_move, lmin = 5, Kmax = 25)

## ----ex ltraj, eval = FALSE---------------------------------------------------
#  segmentation(x_ltraj, lmin = 5, Kmax = 25)

## ----disabling subsampling, fig.show='hold', eval = FALSE---------------------
#  shiftseg <- segmentation(simulshift,
#                           Kmax = 30, lmin=5,
#                           seg.var = c("x","y"),
#                           subsample = FALSE)
#  
#  mode_segclust <- segclust(simulmode,
#                            Kmax = 30, lmin=5, ncluster = c(2,3,4),
#                            seg.var = c("speed","abs_spatial_angle"),
#                            subsample = FALSE)

## ----automatic subsampling, fig.show='hold', eval = FALSE---------------------
#  shiftseg <- segmentation(simulshift,
#                           Kmax = 30, lmin=5,
#                           seg.var = c("x","y"),
#                           subsample_over = 2000)
#  
#  mode_segclust <- segclust(simulmode,
#                            Kmax = 30, lmin=5, ncluster = c(2,3,4),
#                            seg.var = c("speed","abs_spatial_angle"),
#                            subsample_over = 500)

## ----manual subsampling, fig.show='hold', eval = FALSE------------------------
#  shiftseg <- segmentation(simulshift,
#                           Kmax = 30, lmin=5,
#                           seg.var = c("x","y"),
#                           subsample_by = 60)
#  
#  mode_segclust <- segclust(simulmode,
#                            Kmax = 30, lmin=5, ncluster = c(2,3,4),
#                            seg.var = c("speed","abs_spatial_angle"),
#                            subsample_by = 2)

## ----lmin and subsampling, echo = FALSE---------------------------------------
lmin <-  240
subsample_by  <-  60
  cli::cli_alert_success("Using {cli::col_green('lmin = ', lmin )}") 
lmin <- max(lmin/subsample_by,5)

      cli::cli_alert_success(
        "Adjusting lmin to subsampling. 
        {cli::col_grey('Dividing lmin by ',subsample_by,', with a minimum of 5')}")
      cli::cli_alert("After subsampling, {cli::col_green('lmin = ', lmin)}. 
                    {cli::col_grey('Corresponding to lmin = ',lmin*subsample_by,
                     ' on the original time scale')}")


## ----calculate covariate, fig.show='hold', eval = FALSE-----------------------
#  simple_data <- simulmode[,c("dateTime","x","y")]
#  full_data   <- add_covariates(simple_data,
#                                coord.names = c("x","y"),
#                                timecol = "dateTime",
#                                smoothed = TRUE,
#                                units ="min")
#  head(full_data)

## ----example repeated, eval = FALSE-------------------------------------------
#  df <- data.frame(x = rep(1,500), y = rep(2, 500))
#  segclust(df,
#           seg.var = c("x","y"),
#           lmin = 50, ncluster = 3 )

## ----example repeated message, echo = FALSE-----------------------------------
seg.var <- c("x","y")
      cli::cli_alert_danger(
        "Data have repetition of nearly-identical \\
        values longer than lmin. 
        {cli::col_grey('The algorithm cannot estimate variance \\
        for segment with repeated values. \\
        This is potentially caused by interpolation \\
         of missing values or rounding of values.')}
        {cli::symbol$arrow_right} Please check for repeated \\
        or very similar values of {seg.var}")

