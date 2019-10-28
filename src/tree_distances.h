#include <stdint.h>
#include "lap.h" /* for cost */

/*************** TYPES      *******************/

typedef uint32_t splitbit;

/*************** CONSTANTS  *******************/

#define BIN_SIZE 32
#define MAX_BINS 100
#define MAX_TIPS BIN_SIZE * MAX_BINS
#define MAX_PARTITIONS MAX_TIPS - 3

#define BIG (score_t) (2 ^ sizeof(score_t)) / MAX_PARTITIONS /* Use max_tips for power of 2? */
#define BIGL (double) BIG

/*************** FUNCTIONS  *******************/

extern double lg2_trees_matching_split(int a, int b);
