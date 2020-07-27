
<!-- README.md is generated from README.Rmd. Please edit that file -->

# shapegrid

Creates a grid of points within an ellipse. Code taken from
<https://people.sc.fsu.edu/~jburkardt/cpp_src/ellipse_grid/ellipse_grid.html>.

Install package dependencies:

``` r
pkgs = c("Rcpp", "remotes" )
install.packages(pkgs)
```

Then install `shapegrid` package from github:

``` r
remotes::install_github("daffp/shapegrid")
```

Run some examples to see that it is working

``` r
library(shapegrid)
p = shapegrid::ellipse_grid(10, c(1,10), c(1,5))
plot(p, pch=16, cex=0.5)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />
