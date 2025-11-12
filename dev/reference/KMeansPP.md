# k-means++ clustering

k-means++ clustering (Arthur and Vassilvitskii 2007) improves the speed
and accuracy of standard [`kmeans`](https://rdrr.io/r/stats/kmeans.html)
clustering (Hartigan and Wong 1979) by preferring initial cluster
centres that are far from others. A scalable version of the algorithm
has been proposed for larger data sets (Bahmani et al. 2012) , but is
not implemented here.

## Usage

``` r
KMeansPP(x, k = 2, nstart = 10, ...)
```

## Arguments

- x:

  Numeric matrix of data, or an object that can be coerced to such a
  matrix (such as a numeric vector or a data frame with all numeric
  columns).

- k:

  Integer specifying the number of clusters, *k*.

- nstart:

  Positive integer specifying how many random sets should be chosen

- ...:

  additional arguments passed to
  [`kmeans`](https://rdrr.io/r/stats/kmeans.html)

## References

Arthur D, Vassilvitskii S (2007). “K-Means++: The Advantages of Careful
Seeding.” In *Proceedings of the Eighteenth Annual ACM-SIAM Symposium on
Discrete Algorithms*, SODA '07, 1027–1035.  
  
Bahmani B, Moseley B, Vattani A, Kumar R, Vassilvitskii S (2012).
“Scalable K-Means++.” *arXiv*.
[doi:10.48550/arXiv.1203.6402](https://doi.org/10.48550/arXiv.1203.6402)
, 1203.6402.  
  
Hartigan JA, Wong MA (1979). “Algorithm AS 136: a *K*-means clustering
algorithm.” *Journal of the Royal Statistical Society. Series C (Applied
Statistics)*, **28**(1), 100–108.
[doi:10.2307/2346830](https://doi.org/10.2307/2346830) .

## See also

[`kmeans`](https://rdrr.io/r/stats/kmeans.html)

Other cluster functions:
[`cluster-statistics`](https://ms609.github.io/TreeDist/dev/reference/cluster-statistics.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
# Generate random points
set.seed(1)
x <- cbind(c(rnorm(10, -5), rnorm(5, 1), rnorm(10, 6)),
           c(rnorm(5, 0), rnorm(15, 4), rnorm(5, 0)))

# Conventional k-means may perform poorly
klusters <- kmeans(x, cent = 5)
plot(x, col = klusters$cluster, pch = rep(15:19, each = 5))


# Here, k-means++ recovers a better clustering
plusters <- KMeansPP(x, k = 5)
plot(x, col = plusters$cluster, pch = rep(15:19, each = 5))
```
