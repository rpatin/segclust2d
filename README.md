# segtools

`segtools` provides R code for two methods of segmentation and joint
segmentation/clustering of bivariate time-series. It was originally intended for
ecological segmentation (home-range and behavioural modes) but can be easily
applied on other type of time-series. The package also provides tools for
analysing outputs from R packages `moveHMM` and `marcher`.

## Installation

You can install `segtools` from github with:

``` r
# install.packages("devtools")
devtools::install_github("rpatin/segtools")
```

## Examples

The algorithm can perform a [segmentation](#Segmentation) of the time-serie into
homogeneous segments. A typical case is the identification of home-range
behaviour. It can also perform an integrated classification of those segments
into clusters of homogeneous behaviour through a
[segmentation/clustering](#Segmentation-Clustering) algorithm. This is generally used to
identify behavioural modes. Input data can be a `data.frame` (shown in the first examples), a `Move` object or a `ltraj` object (from package `adehabitatLT`), both shown in section [Other data types](#Other-data-types)

### Segmentation

``` r
library(segtools)
data(simulshift)

head(simulshift)

```

`simulshift` is an example dataset containing a simulation of home-range behaviour with two shifts. It is a data.frame with two columns for coordinates : x and y. We can now run a simple segmentation with this dataset to find the different home-ranges. By default `segmentation()` function will take `type` argument as 'home-range' and choose default column x and y as coordinates but one can specify something else using argument `seg.var`.

For a segmentation one has to specify arguments `lmin`, the minimum length of a segment and `Kmax`, the maximum number of segments. By default `Kmax` will be set to `floor(n/lmin)`, with `n` the number of observations. However this can considerably slow the calculations so do not hesitate to reduce it to a reasonable value.

``` r
shift_seg <- segmentation(simulshift, lmin = 5, Kmax = 25, type = "home-range")
```

Segmentation is performed through a Dynamic Programming algorithm that finds the best segmentation given a number of segment. For each number of segment, the optimal segmentation is associated with a likelihood value. By default, the algorithm choose the number of segment given a criterium developped by Marc Lavielle based on the value of the second derivative of the penalized likelihood. This criterium use a threshold value of `S = 0.75`, but a different threshold can be specified.

`segmentation()` returns an object of `segmentation-class` for which several methods are available (see section [segmentation-class](#segmentation-class))



``` r
plot(shift_seg)
```



``` r
plot(shift_seg)
```

### Segmentation-Clustering



``` r
plot(shift_seg)
```

### Other data types

### segmentation-class

#### extract informations

##### augment  

##### states 

##### log-Likelihood - logLik 

##### get_BIC

#### Graphical outputs

##### plot.segmentation


``` r
plot(shift_seg)
```
##### segmap


``` r
segmap(shift_seg)
```
##### stateplot

##### plot_likelihood

##### plot_BIC


[segmentation](#Segmentation) 
