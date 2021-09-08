# Introduction

`segclust2d` provides R code for two methods of segmentation and joint
segmentation/clustering of bivariate time-series. It was originally intended for
ecological segmentation (home-range and behavioural modes) but can be easily
applied on other type of time-series. The package also provides tools for
analysing outputs from R packages `moveHMM` and `marcher`.

# Website

Full documentation for segclust2d is available on its website:
[https://rpatin.github.io/segclust2d/](https://rpatin.github.io/segclust2d/)

Three differents topic are discussed there, and are also available as vignettes in the R package:
- First, preparation of data to be analysed with segmentation/clustering algorithm, including covariate calculations, subsampling and guidelines to prepare movement data.
- Second, a guide to run the segmentation and segmentation/clustering method, including advise on setting parameters
- Finally, an overview of the possible outputs that can be generated with the package.

# Installation

For the CRAN version : 
``` r
install.packages("segclust2d")
```

If you want the newest version, you can install `segclust2d` from github with:

``` r
devtools::install_github("rpatin/segclust2d")
```
