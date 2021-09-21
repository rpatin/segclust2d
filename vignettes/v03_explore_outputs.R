## ---- echo = FALSE------------------------------------------------------------
options(Encoding="UTF-8")
knitr::opts_chunk$set(
  fig.width = 8,
  fig.height = 5,
  collapse = TRUE,
  comment = "#>"
)

## ----library and data, fig.show='hold'----------------------------------------
library(segclust2d)

## ----loading data and segclust, fig.show='hold', message = FALSE--------------
data(simulmode)
simulmode$abs_spatial_angle <- abs(simulmode$spatial_angle)
simulmode <- simulmode[!is.na(simulmode$abs_spatial_angle), ]
mode_segclust <- segclust(simulmode,
                          Kmax = 20, lmin=10, ncluster = c(2,3),
                          seg.var = c("speed","abs_spatial_angle"),
                          scale.variable = TRUE)

## ----plot.segmentation, fig.show='hold', message = FALSE----------------------
plot(mode_segclust, ncluster = 3)

## ----segmap, fig.show='hold', message = FALSE, fig.width=5, fig.height=5------
segmap(mode_segclust, ncluster = 3)

## ----stateplot, fig.show='hold', message = FALSE, fig.width=3, fig.height=3----
stateplot(mode_segclust, ncluster = 3)

## ----augment, fig.show='hold', message = FALSE, eval = FALSE------------------
#  augment(mode_segclust, ncluster = 3)

## ----segment, fig.show='hold', results = "hide", message = FALSE--------------
segment(mode_segclust, ncluster = 3)

## ----states, fig.show='hold', results = "hide", message = FALSE---------------
states(mode_segclust, ncluster = 3)

## ----simulshift, fig.show='hold', message = FALSE-----------------------------
data("simulshift")
shift_seg <- segmentation(simulshift,
                          seg.var = c("x","y"),
                          lmin = 240, Kmax = 25,
                          subsample_by = 60)

## ----logLik, eval = FALSE-----------------------------------------------------
#  logLik(shift_seg)
#  

## ----plot_likelihood, fig.show='hold', message = FALSE------------------------
plot_likelihood(shift_seg)

## ----BIC, fig.show='hold', results = "hide", message = FALSE------------------
BIC(mode_segclust)

## ----plot_BIC, fig.show='hold', message = FALSE-------------------------------
plot_BIC(mode_segclust)

