segclust2d: bivariate segmentation with optional clustering for R
================

[![](https://www.r-pkg.org/badges/version/segclust2d?color=orange)](https://cran.r-project.org/package=segclust2d)
[![](http://cranlogs.r-pkg.org/badges/grand-total/segclust2d?color=yellow)](https://cran.r-project.org/package=segclust2d)
[![](https://img.shields.io/badge/devel%20version-0.3.1-blue.svg)](https://github.com/rpatin/segclust2d)
[![](https://img.shields.io/github/last-commit/rpatin/segclust2d.svg)](https://github.com/rpatin/segclust2d/commits/master)

# Introduction

`segclust2d` provides R code for a segmentation method that can be used
on all bivariate time-series. The segmentation method can additionally
be associated with a clustering algorithm. It was originally intended
for ecological segmentation (home-range and behavioural modes) but can
be easily applied on other type of time-series. The package also
provides tools for analysing outputs from R packages `moveHMM` and
`marcher`.

# Website

Full documentation for segclust2d is available on this website:
<https://rpatin.github.io/segclust2d/>

Three topics are discussed there, and are also available as vignettes in
the R package:

-   [First, preparation of data to be analyzed with
    segmentation/clustering algorithm, including covariate calculations,
    subsampling and guidelines to prepare movement
    data.](https://rpatin.github.io/segclust2d/articles/v01_preparing_data.html)
-   [Second, a guide to run the segmentation and segmentation/clustering
    method, including advise on setting
    parameters.](https://rpatin.github.io/segclust2d/articles/v02_run_segclust2d.html)
-   [Finally, an overview of the possible outputs that can be generated
    with the
    package.](https://rpatin.github.io/segclust2d/articles/v03_explore_outputs.html)

# Installation

For the
[![](https://www.r-pkg.org/badges/version/segclust2d?color=orange)](https://cran.r-project.org/package=segclust2d)
version :

``` r
install.packages("segclust2d")
```

If you want the newest
[![](https://img.shields.io/badge/devel%20version-0.3.1-blue.svg)](https://github.com/rpatin/segclust2d),
you can install `segclust2d` from github with:

``` r
devtools::install_github("rpatin/segclust2d")
```
