## ----option chunk, echo = FALSE-----------------------------------------------
options(Encoding="UTF-8")
knitr::opts_chunk$set(
  fig.width = 8,
  fig.height = 5,
  collapse = TRUE,
  comment = "#>"
)

## ----load segclust2d----------------------------------------------------------
library(segclust2d)

## ----load simulshift, fig.show='hold'-----------------------------------------
data(simulshift)

## ----plot simulshift, echo = FALSE, fig.width = 5,  fig.height = 4, fig.cap = "`simulshift`: simulation of movement within three successive home-range. Data shown after subsampling by 100."----
library(ggplot2)
tmpdf <- simulshift[seq(1,30000, by = 100),]
tmpdf$class <- factor(rep(c(1,2,3), each = 102)[1:300])
ggplot(tmpdf)+
  geom_path(aes(x = x, y = y))+
  geom_point(aes(x= x, y = y, col = factor(class)))+
  scale_color_discrete("home-range")+
  theme(legend.position= "top")

## ----load simulmode, fig.show='hold'------------------------------------------
data(simulmode)
simulmode$abs_spatial_angle <- abs(simulmode$spatial_angle)
simulmode <- simulmode[!is.na(simulmode$abs_spatial_angle), ]

## ----map simulmode, echo = FALSE, fig.width = 5,  fig.height = 4, fig.cap = "`simulmode`: simulation of movement with three different behavioural modes."----

library(ggplot2)
tmpdf <- simulmode
tmpdf$class <- factor(rep(c(1,2,3), each = 20, 5))
ggplot(tmpdf)+
  geom_path(aes(x = x, y = y))+
  geom_point(aes(x= x, y = y, col = factor(class)))+
  scale_color_discrete("behavioural mode")+
  theme(legend.position= "top")

## ----wrong type argument call, eval = FALSE-----------------------------------
#  df.seg <- segmentation(simulshift,
#                         type = "home-range",
#                         lmin = 300, Kmax = 10,
#                         subsample_by = 60)

## ----wrong type argument message, echo=FALSE, results='asis'------------------
        cli::cli_alert_danger("Argument {cli::col_red('type')} \\
                     is deprecated and should not be used")
        cli::cli_alert_danger("Argument {cli::col_red('coord.names')} \\
                     is deprecated and should not be used")
        coord.names <- c("x","y")
        cli::cli_alert("Please use instead \\
                    {.field seg.var = {deparse(coord.names)}} and \\
                    {.field scale.variable = FALSE}")

## ----minimal segmentation call, eval = FALSE, message = FALSE-----------------
#  shift_seg <- segmentation(simulshift,
#                            lmin = 240,
#                            subsample_by = 60)

## ----Kmax segmentation call, eval = TRUE, message = FALSE---------------------
shift_seg <- segmentation(simulshift,
                          lmin = 240, Kmax = 25,
                          subsample_by = 60)

## ----segmentation warning lmin*Kmax, echo=FALSE, results='asis'---------------
Kmax = 25
      cli::cli_alert_warning(
        "Adjusting Kmax so that lmin*Kmax < nrow(x). Now, \\
      {cli::col_yellow('Kmax = ', Kmax)}")
    

## ----segmentation error lmin*Kmax, error = TRUE, echo=FALSE, results='asis'----
        cli::cli_alert_danger(
          "lmin*Kmax > nrow(x) and Kmax cannot be adjusted. \\
          Please provide lower values for lmin")
   stop("lmin*Kmax > nrow(x)")

## ----providing seg.var, eval = FALSE, message = FALSE-------------------------
#  shift_seg <- segmentation(simulshift,
#                            seg.var = c("x","y"),
#                            lmin = 240, Kmax = 25,
#                            subsample_by = 60)

## ----segmentation summary, echo = FALSE, results = 'asis'---------------------
   cli::cli_alert_success("Best segmentation estimated with \\
                     {shift_seg$Kopt.lavielle} segments, \\
                     according to Lavielle's criterium")
   cli::cli_text(cli::col_grey(
   'Other number of segments may be selected 
     by looking for likelihood breaks with plot_likelihood()'))
   cli::cli_text(cli::col_grey(
   'Results of the segmentation may be explored with plot() and segmap()'))


## ----plot_likelihood example, fig.show = 'hold', fig.width = 4, fig.height = 3----
plot_likelihood(shift_seg)

## ----plot_likelihood mode_seg, fig.show = 'hold', fig.width = 4, fig.height = 3, message = FALSE----

mode_seg <- segmentation(simulmode,
                          lmin = 10, Kmax = 20,
                         seg.var = c("speed","abs_spatial_angle"),
                          scale.variable = TRUE)

plot_likelihood(mode_seg)

## ----simulmode segclust default, fig.show='hold', message = FALSE-------------
mode_segclust <- segclust(simulmode,
                          Kmax = 20, lmin=10, 
                          ncluster = c(2,3),
                          seg.var = c("speed","abs_spatial_angle"))

## ----simulmode segclust scale variable, fig.show='hold', eval = FALSE---------
#  mode_segclust <- segclust(simulmode,
#                            Kmax = 20, lmin=10,
#                            ncluster = c(2,3),
#                            seg.var = c("speed","abs_spatial_angle"),
#                            scale.variable = TRUE)

## ----segclust summary, echo = FALSE, results = 'asis'-------------------------
  cli::cli_alert_success(
    "Best segmentation/clustering estimated with \\
    {mode_segclust$ncluster.BIC} clusters and \\
    {mode_segclust$Kopt.BIC[mode_segclust$ncluster.BIC]} segments according to BIC")
  cli::cli_text(cli::col_grey(
    '{cli::symbol$arrow_right} Number of clusters should preferentially be selected 
    according to biological knowledge. Exploring the BIC plot with plot_BIC()
    can also provide advice to select the number of clusters.'))
  cli::cli_text(cli::col_grey(
    '{cli::symbol$arrow_right} Once number of clusters is selected, \\
    the number of segments can be selected according to BIC.'))
  
  cli::cli_text(cli::col_grey(
    '{cli::symbol$arrow_right} Results of the segmentation/clustering
    may further be explored with plot() and segmap()'))
  

## ----plot_BIC, fig.show='hold', message = FALSE-------------------------------
plot_BIC(mode_segclust)

## ----simulmode plot BIC cluster 2-5, fig.show='hold', message = FALSE---------
mode_segclust <- segclust(simulmode,
                          Kmax = 20, lmin=10, 
                          ncluster = 2:5,
                          seg.var = c("speed","abs_spatial_angle"),
                          scale.variable = TRUE)

plot_BIC(mode_segclust)

## ----3 clusters, fig.show='hold', message = FALSE-----------------------------
plot(mode_segclust, ncluster = 3)

## ----4 clusters, fig.show='hold', message = FALSE-----------------------------
plot(mode_segclust, ncluster = 4)

