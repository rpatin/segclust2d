#' Simulations of home-range shift
#'
#' A dataset containing a simulation of home-range shift
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{indice}{index of position}
#'   \item{x}{x coordinates}
#'   \item{y}{y coordinates}
#'   \item{dateTime}{arbitrary date in POSIXct format}
#' }

"simulshift"

# simulshift =   ACS <- read.table("../../These/Segmentation/SegTraj Draft/data/Raw/shift_sensibility/DV_shift_mix_10.txt",
#                                 col.names = c("indice", "x", "y"))
# 
# simulshift$dateTime <- as.POSIXct(simulshift$indice*60, origin = Jan1)
# 
# devtools::use_data(simulshift, overwrite = TRUE)

#' Simulations of behavioural mode
#'
#' A dataset containing a simulation of 3 different behavioural mode
#'
#' @format A data frame with 302 rows and 10 variables:
#' \describe{
#'   \item{indice}{index of position}
#'   \item{x}{x coordinates}
#'   \item{y}{y coordinates}
#'   \item{speed}{smoothed speed}
#'   \item{spatial_angle}{angle at constant step length}
#'   \item{dist}{raw speed}
#'   \item{angle}{angular speed}
#'   \item{vit_p}{persistence speed}
#'   \item{vit_r}{rotation speed}
#'   \item{vit_p_spa}{persistence speed calculated with spatial angles}
#'   \item{vit_r_spa}{rotation speed calculated with spatial angles}
#'   \item{dateTime}{arbitrary date in POSIXct format}
#' }

"simulmode"

# simulmode =   ACS <- read.table("../../These/Segmentation/SegTraj Draft/data/Raw/acs_sensibility/acs3_1_25.txt",
#                                 col.names = c("indice", "x", "y", "speed", "spatial_angle",
#                                               "dist", "angle", "vit_p", "vit_r", "vit_p_spa",
#                                               "vit_r_spa"),
#                                 na.strings = "9999.000", fill=TRUE)
# simulmode$dateTime <- as.POSIXct(simulmode$indice*60, origin = Jan1)
# devtools::use_data(simulmode, overwrite = TRUE)
