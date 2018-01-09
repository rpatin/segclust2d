# segclust2d

`segclust2d` provides R code for two methods of segmentation and joint
segmentation/clustering of bivariate time-series. It was originally intended for
ecological segmentation (home-range and behavioural modes) but can be easily
applied on other type of time-series. The package also provides tools for
analysing outputs from R packages `moveHMM` and `marcher`.

## Installation

You can install `segclust2d` from github with:

``` r
# install.packages("devtools")
devtools::install_github("rpatin/segclust2d")
```

## Examples

The algorithm can perform a [segmentation](#segmentation) of the time-serie into
homogeneous segments. A typical case is the identification of home-range
behaviour. It can also perform an integrated classification of those segments
into clusters of homogeneous behaviour through a
[segmentation/clustering](#segmentation-clustering) algorithm. This is generally used to
identify behavioural modes. Input data can be a `data.frame` (shown in the first examples), a `Move` object or a `ltraj` object (from package `adehabitatLT`), both shown in section [Other data types](#other-data-types)

### Segmentation

``` r
library(segclust2d)
data(simulshift)
```

`simulshift` is an example dataset containing a simulation of home-range behaviour with two shifts. It is a data.frame with two columns for coordinates : x and y. We can now run a simple segmentation with this dataset to find the different home-ranges. By default `segmentation()` function will take `type` argument as 'home-range' and choose default column x and y as coordinates but one can specify something else using argument `seg.var`.

For a segmentation one has to specify arguments `lmin`, the minimum length of a segment and `Kmax`, the maximum number of segments. By default `Kmax` will be set to `floor(n/lmin)`, with `n` the number of observations. However this can considerably slow the calculations so do not hesitate to reduce it to a reasonable value.

``` r
shift_seg <- segmentation(simulshift, lmin = 5, Kmax = 25, type = "home-range")
```

Segmentation is performed through a Dynamic Programming algorithm that finds the best segmentation given a number of segment. For each number of segment, the optimal segmentation is associated with a likelihood value. By default, the algorithm choose the number of segment given a criterium developped by Marc Lavielle based on the value of the second derivative of the penalized likelihood. This criterium use a threshold value of `S = 0.75`, but a different threshold can be specified. By default segmentation is done on coordinates `c("x","y")` but one can specify other names through arguments `coord.names`.

`segmentation()` returns an object of `segmentation-class` for which several methods are available (see section [segmentation-class](#segmentation-class)). The most important one is plot.segmentation, that shows the segmented time-series. 

``` r
plot(shift_seg)
```
By default,  `plot.segmentation` shows the best segmentation, but one can specify a given number of segments (inside the range `1:Kmax`). See [segmentation-class](#plot.segmentation) for additional informations.

``` r
plot(shift_seg, nseg = 10)
```

The second important method is `plot_likelihood` that shows the log-likelihood of the best segmentation versus the number of segments and highlights the one chosen with Lavielle's criterium.


``` r
plot_likelihood(shift_seg, nseg = 10)
```

### Segmentation-Clustering

``` r
data(simulmode)
simulmode <- simulmode[!is.na(simulmode$spatial_angle), ]
```

`simulmode` is an example dataset containing a movement simulation with three different movement mode. It is a data.frame with 11 columns, with coordinates and several covariates. Be careful to check your dataset for missing value.

We can now run a joint segmentation/clustering on this dataset to identify the different behavioural modes. By default `segclust()` function will take `type` argument as 'behavior' but one has to specify two variables for segmentation using argument `seg.var`. 

For a joint segmentation/clustering one has to specify arguments `lmin`, the minimum length of a segment and `Kmax`, the maximum number of segments, and `ncluster` a vector of number of class. By default `Kmax` will be set to `floor(n/lmin)`, with `n` the number of observations. However this can considerably slow the calculations so do not hesitate to reduce it to a reasonable value.

By default `segclust()` will standardized variables, but one can change this by setting `scale.variable = F`. 

``` r
mode_segclust <- segclust(simulmode, Kmax = 30, lmin=5, ncluster = c(2,3,4), type = "behavior", seg.var = c("speed","spatial_angle"))

```

`segclust()` also an object of `segmentation-class` for which the same methods are available (see section [segmentation-class](#segmentation-class)). The most important one is again plot.segmentation, that shows the segmented time-series. 

``` r
plot(mode_segclust)
```

By default for a segmentation/clustering,  `plot.segmentation` shows the best segmentation, maximizing BIC-penalized likelihood, but one can specify a given number of cluster/or segment See [segmentation-class](#plot.segmentation) for additional informations.

``` r
plot(mode_segclust, ncluster = 3)
plot(mode_segclust, ncluster = 3, nseg = 7)

```
One can also inspect the BIC-penalized log-likelihood through functions `plot_BIC()`.

``` r
plot_BIC(mode_segclust)
```

## Other data types

We have shown examples for using data.frames but one can also segment data from `ltraj` and `Move` object that contains a single individual.

### Concerning segmentation

For a simple segmentation, the algorithm will assume a home-range segmentation and use coordinates directly.

``` r
segmentation(ltraj_object, lmin = 5, Kmax = 25)
segmentation(Move_object, lmin = 5, Kmax = 25)
```
### Concerning segclust

For a segmentation/clustering, one has to provide the variables used for segmentation

``` r
segmentation(ltraj_object, lmin = 5, Kmax = 25, ncluster = c(2,3), seg.var = c("speed","spatial_angle"))
segmentation(Move_object, lmin = 5, Kmax = 25, ncluster = c(2,3), seg.var = c("speed","spatial_angle"))
```

Of course the variable names provided must exist as column in `Move_object@data` and `adehabitatLT::infolocs(ltraj_object[1])`.

## segmentation-class

Both functions `segmentation()` and `segclust()` returns a `segmentation-class` object for which several methods are available.

### extract informations

#### augment  

`augment.segmentation()` is a method for `broom::augment`. It returns an augmented data.frame with outputs of the model - here, the attribution to segment or cluster

``` r
augment(shift_seg)
augment(mode_segclust)
```
By default `augment.segmentation` will use data for the best segmentation (maximum of BIC for `segclust()` and Lavielle's criterium for `segmentation()`) but one can ask for a specific segmentation

``` r
augment(shift_seg, nseg = 10) # segmentation()
augment(mode_segclust, ncluster = 2) # segclust()
augment(mode_segclust, ncluster = 2, nseg = 5) # segclust()
```

#### segment

`segment()` allows retrieving informations on the different segment of a given segmentation. Each segment is associated with the mean and standard deviation for each variable, the state (equivalent to the segment number for `segmentation`) and the state ordered given a variable - by default the first variable given by `seg.var`. One can specify the variable for ordering states through the `order.var` of `segmentation()` and `segclust()`.

``` r
segment(shift_seg)
segment(shift_seg, nseg = 3)
segment(mode_segclust)
segment(mode_segclust, nclust = 3, nseg = 8)
```

#### states 

`states()` return information on the different states of the segmentation. For `segmentation()` it is quite similar to `segment()`. For `segclust`, however it gives the different cluster found and the statistics associated.

``` r
states(shift_seg)
states(shift_seg, nseg = 3)
states(mode_segclust)
states(mode_segclust, nclust = 3, nseg = 8)
```

#### log-Likelihood - logLik 

`logLik.segmentation()` return information on the log-likelihood of the different segmentations possible. It returns a data.frame with the number of segment, the log-likelihood and eventually the number of cluster.

``` r
logLik(shift_seg)
logLik(mode_segclust)
```

#### BIC (segclust)

`BIC.segmentation()` return information on the BIC-penalized log-likelihood of the different segmentations possible. It returns a data.frame with the number of segment, the BIC-penalized log-likelihood and the number of cluster. For `segclust()` only.

``` r
BIC(mode_segclust)
```

### Graphical outputs

`segmentation-class` also provides methods for plotting results of segmentations. All plot methods use ggplot2 library.

#### plot.segmentation

`plot.segmentation()` can be used to plot the output of a segmentation as a serie-plot. A specific segmentation can be chosen with `nseg` and `ncluster` arguments. If the original data had a specific x-axis, like a `POSIXct` time column, this can be specified using argument `xcol`. By default, data are plotted by their number. If you want clusters or segments to be ordered according to one of the variables, this can be specified using argument `order`. By default segmentation/clustering output are plotted using ordered states.

``` r
plot(shift_seg)
plot(mode_segclust, ncluster = 3, nseg = 10, xcol = "indice", order = T)
```

#### segmap

`segmap()` plot the results of the segmentation as a map. This can be done only if data have a geographic meaning. Coordinate names are by default "x" and "y" but this can be provided through argument `coord.names`.

``` r
segmap(shift_seg, nseg = 10)
segmap(mode_segclust, ncluster = 3, nseg = 10)
```

#### stateplot

`stateplot()` show statistics for each state or segment.
``` r
stateplot(shift_seg, nseg = 10)
stateplot(mode_segclust, ncluster = 3, nseg = 10)
```


#### plot_likelihood

`plot_likelihood()` plot the log-likelihood of the segmentation for all the tested number of segments and clusters.
``` r
plot_likelihood(shift_seg)
plot_likelihood(mode_segclust)
```

#### plot_BIC

`plot_likelihood()` plot the BIC-penalized log-likelihood of the segmentation for all the tested number of segments and clusters.
``` r
plot_BIC(mode_segclust)
```

## Covariate calculations

