## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 8,
  fig.height = 5,
  collapse = TRUE,
  comment = "#>"
)

## ----library and data, fig.show='hold'----------------------------------------
library(segclust2d)
data(simulmode)
simulmode$abs_spatial_angle <- abs(simulmode$spatial_angle)
simulmode <- simulmode[!is.na(simulmode$abs_spatial_angle), ]

## ---- fig.show='hold', message = FALSE----------------------------------------
mode_segclust <- segclust(simulmode,
                          Kmax = 20, lmin=10, 
                          ncluster = c(2,3),
                          seg.var = c("speed","abs_spatial_angle"))

