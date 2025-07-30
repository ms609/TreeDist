ub <- microbenchmark::microbenchmark
set.seed(1)
UnifMat <- function(n) matrix(runif(n * n), n, n)
test40 <- UnifMat(40)
test400 <- UnifMat(400)
test2000 <- UnifMat(2000)

ub(LAPJV(test40), LAPJV(test400), times = 100)
# With stl vectors in caller only
#            expr    min     lq     mean  median      uq     max neval
#  LAPJV(test40)   24.8   26.4   32.784   28.75   37.50    68.2   100
# LAPJV(test400) 3356.5 3486.1 3992.628 3572.05 3730.15 39626.6   100
# 
#   LAPJV(test40)   24.7   26.25   29.720   28.90   31.30   50.1   100
#  LAPJV(test400) 3360.7 3393.65 3499.819 3455.85 3542.15 3852.9   100

#   LAPJV(test40)   24.7   26.30   30.111   27.90   32.05   44.0   100
#  LAPJV(test400) 3283.0 3391.35 3526.312 3451.55 3630.90 4359.1   100
 
# stl vectors throughout (again after install)

#           expr    min      lq     mean  median     uq    max neval
#  LAPJV(test40)   25.0   26.65   30.214   29.25   31.3   93.8   100
# LAPJV(test400) 3281.2 3304.60 3436.435 3442.10 3490.9 4078.0   100

# LAPJV(test40)   24.4   26.45   32.038   29.70   32.6  154.6   100
# LAPJV(test400) 3290.8 3384.45 3544.428 3451.15 3577.0 5605.4   100

# no stl vectors

#             expr    min        lq     mean      median      uq     max neval
# LAPJV(test40)    24.001   25.8515   29.31508   28.851   31.101   43.101   100
# LAPJV(test400) 3228.601 3321.4505 3401.51697 3383.201 3438.751 3843.301   100
# LAPJV(test40)    23.800   26.001   30.07794   29.1005   30.651   81.901   100
# LAPJV(test400) 3300.702 3332.251 3500.24799 3393.8510 3528.001 5721.001   100

# Using flat template

#           expr      min       lq       mean    median        uq      max neval
#  LAPJV(test40)   23.702   25.151   30.58508   27.4515   31.8515   65.901   100
# LAPJV(test400) 3537.501 3640.601 3757.95597 3718.6505 3839.1500 4362.601   100

ub(LAPJV(test2000), times = 25)
# With stl vectors in caller only
# Unit: milliseconds
#             expr      min      lq     mean   median       uq      max neval
#  LAPJV(test2000) 113.3613 116.661 123.9962 121.1821 126.8234 190.1466    25
#             expr      min       lq     mean   median       uq      max neval
#  LAPJV(test2000) 114.1162 116.0226 121.3515 117.5171 125.5254 157.3014    25
# 
# stl vectors throughout
# Unit: milliseconds
#             expr      min       lq     mean   median       uq      max neval
#  LAPJV(test2000) 113.9332 115.7949 118.2589 116.5993 117.7188 144.8049    25
#  LAPJV(test2000) 112.4989 117.2052 120.0313 118.7492 120.4724 159.197    25
#  LAPJV(test2000) 112.1221 113.8703 116.5592 114.8536 118.8638 123.6153    25
#  LAPJV(test2000) 114.2271 116.5065 121.3377 118.9613 120.2638 187.0768    25
#  LAPJV(test2000) 114.514  117.4176 123.4273 118.4348 121.0561 180.2217    25

# Main, with no stl vectors

# LAPJV(test2000) 115.1755 117.8086 126.4238 122.9364 125.3566 200.8201    25
# LAPJV(test2000) 110.4296 112.9714 115.2006 114.5085 117.2046 123.2914    25
# LAPJV(test2000) 109.5763 112.2592 121.7592 114.8587 117.5283 203.9729    25
# LAPJV(test2000) 111.3518 113.1852 114.5551 114.7754 115.3454 120.0273    25

