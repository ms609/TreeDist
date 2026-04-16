#ifndef _TREEDIST_INT_H
#define _TREEDIST_INT_H

#include <TreeDist/types.h>

// Re-export shared types to global scope for backward compatibility.
using TreeDist::int16;
using TreeDist::int32;
using TreeDist::split_int;

// Types used only within TreeDist's own source.
using uint32 = uint_fast32_t;
using grf_match = std::vector<int>;

#endif
