#include <stdint.h>
#include "lap.h" /* for cost */

/*************** TYPES      *******************/

typedef uint64_t splitbit;

/*************** CONSTANTS  *******************/

#define BIN_SIZE 64
#define MAX_BINS 32
#define MAX_TIPS BIN_SIZE * MAX_BINS
#define MAX_SPLITS MAX_TIPS /* -3, but quicker if a power of two? */
#define ALL_ONES (splitbit) ~((splitbit) 0U)

#define BIG cost ((cost) 100000000 * (cost) 10000000) /* TODO Calculate */
#define BIGL (double) BIG

const splitbit right16bits = 65535U;
const uint32_t powers_of_two[16] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512,
                                    1024, 2048, 4096, 8192, 16384, 32768};

/*************** FUNCTIONS  *******************/

extern double lg2_trees_matching_split(int a, int b);

extern int count_bits (splitbit x);
