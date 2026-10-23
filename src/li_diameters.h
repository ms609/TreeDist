#ifndef _TREEDIST_LI_DIAMETERS_H
#define _TREEDIST_LI_DIAMETERS_H

// Exact diameter of NNI tree space on 0..N_EXACT leaves: the largest NNI
// distance between any two trees with that many leaves.
// Values for <= 11 leaves from Li et al. (1996, table 1); 12 leaves from
// Fack et al. (2002, table 2, n = 10 internal nodes).
// A diameter bounds the distance between a given pair from above, not below.
constexpr int32_t N_EXACT = 12;
constexpr int32_t exact_diameter[] = {0, 0, 0, 0, 1, 3, 5, 7, 10, 12, 15, 18, 21};

#endif
