# Solve linear assignment problem using LAPJV

Use the algorithm of Jonker and Volgenant (1987) to solve the Linear Sum
Assignment Problem (LSAP).

## Usage

``` r
LAPJV(x)
```

## Arguments

- x:

  Matrix of costs.

## Value

`LAPJV()` returns a list with two entries: `score`, the score of the
optimal matching; and `matching`, the columns matched to each row of the
matrix in turn.

## Details

The Linear Assignment Problem seeks to match each row of a matrix with a
column, such that the cost of the matching is minimized.

The Jonker & Volgenant approach is a faster alternative to the Hungarian
algorithm (Munkres 1957) , which is implemented in `clue::solve_LSAP()`.

Note: the JV algorithm expects integers. In order to apply the function
to a non-integer *n*, as in the tree distance calculations in this
package, each *n* is multiplied by the largest available integer before
applying the JV algorithm. If two values of *n* exhibit a trivial
difference – e.g. due to floating point errors – then this can lead to
interminable run times. (If numbers of the magnitude of billions differ
only in their last significant digit, then the JV algorithm may undergo
billions of iterations.) To avoid this, integers over 2^22 that differ
by a value of 8 or less are treated as equal.

## References

Jonker R, Volgenant A (1987). “A shortest augmenting path algorithm for
dense and sparse linear assignment problems.” *Computing*, **38**,
325–340. [doi:10.1007/BF02278710](https://doi.org/10.1007/BF02278710)
.  
  
Munkres J (1957). “Algorithms for the assignment and transportation
problems.” *Journal of the Society for Industrial and Applied
Mathematics*, **5**(1), 32–38.
[doi:10.1137/0105003](https://doi.org/10.1137/0105003) .

## See also

Implementations of the Hungarian algorithm exist in adagio,
RcppHungarian, and clue and lpSolve; for larger matrices, these are
substantially slower. (See discussion at [Stack
Overflow](https://stackoverflow.com/questions/72806265/).)

The JV algorithm is implemented for square matrices in the Bioconductor
package
[`GraphAlignment::LinearAssignment()`](https://www.bioconductor.org/packages/release/bioc/html/GraphAlignment.html).

## Author

[C++
code](https://github.com/yongyanghz/LAPJV-algorithm-c/blob/master/src/lap.cpp)
by Roy Jonker, MagicLogic Optimization Inc. <roy_jonker@magiclogic.com>,
with contributions from Yong Yang <yongyanglink@gmail.com>, after [Yi
Cao](https://uk.mathworks.com/matlabcentral/profile/authors/69713-yi-cao)

## Examples

``` r
problem <- matrix(c(7, 9, 8, 9, 9,
                    2, 8, 5, 7, 9,
                    1, 6, 6, 9, 9,
                    3, 6, 2, 2, 9), 4, 5, byrow = TRUE)

LAPJV(problem)
#> $score
#> [1] 17
#> 
#> $matching
#> [1] 2 3 1 4
#> 
```
