# SPR approximation via Oliveira Martins *et al.* (2008)

SPR approximation via Oliveira Martins *et al.* (2008)

## Usage

``` r
.SPRPairDeO(tree1, tree2, check = TRUE)
```

## Examples

``` r
# de Oliveira Martins et al 2008, fig. 7
tree1 <- ape::read.tree(text = "((1, 2), ((a, b), (c, d)), (3, (4, (5, (6, 7)))));")
tree2 <- ape::read.tree(text = "((1, 2), 3, (4, (5, (((a, b), (c, d)), (6, 7)))));")
oPar <- par(mfrow =c(2, 1), mar = rep(0, 4))
plot(tree1)
plot(tree2)

par(oPar)
SPRDist(tree1, tree2, method = "deO")
#> [1] 1
```
