# Variation of information for all split pairings

Calculate the variation of clustering information (Meila 2007) for each
possible pairing of non-trivial splits on *n* leaves (Smith 2020) ,
tabulating the number of pairings with each similarity.

## Usage

``` r
AllSplitPairings(n)
```

## Arguments

- n:

  Integer specifying the number of leaves in a tree.

## Value

`AllSplitPairings()` returns a named vector. The name of each element
corresponds to a certain variation of information, in bits; the value of
each element specifies the number of pairings of non-trivial splits that
give rise to that variation of information. Split `AB|CD` is treated as
distinct from `CD|AB`. If pairing `AB|CD`=`CD|AB` is considered
equivalent to `CD|AB`=`CD|AB` (etc), then values should be divided by
four.

## References

Meila M (2007). “Comparing clusterings—an information based distance.”
*Journal of Multivariate Analysis*, **98**(5), 873–895.
[doi:10.1016/j.jmva.2006.11.013](https://doi.org/10.1016/j.jmva.2006.11.013)
.  
  
Smith MR (2020). “Information theoretic Generalized Robinson-Foulds
metrics for comparing phylogenetic trees.” *Bioinformatics*, **36**(20),
5007–5013.
[doi:10.1093/bioinformatics/btaa614](https://doi.org/10.1093/bioinformatics/btaa614)
.

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
AllSplitPairings(6)
#>                0                1 1.33333333333333 1.74837083261218 
#>              100              480              360              480 
#> 1.83659166810898 1.91829583405449 
#>              360              720 
# Treat equivalent splits as identical by dividing by four:
AllSplitPairings(6) / 4L
#>                0                1 1.33333333333333 1.74837083261218 
#>               25              120               90              120 
#> 1.83659166810898 1.91829583405449 
#>               90              180 
```
