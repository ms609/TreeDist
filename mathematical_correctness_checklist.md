# Checklist for Demonstrating Mathematical Correctness of HierarchicalMutualInfoDist

This document provides a checklist of steps to verify the mathematical correctness of the hierarchical mutual information distance implementation.

## Basic Mathematical Properties

### 1. Non-negativity
- [ ] **Test**: For any two trees T1 and T2, verify that HierarchicalMutualInfoDist(T1, T2) ≥ 0
- **Implementation**: All distance calculations use `max(0, distance)` to ensure non-negativity
- **Verification**: Run tests with various tree pairs

### 2. Identity of Indiscernibles  
- [ ] **Test**: For any tree T, verify that HierarchicalMutualInfoDist(T, T) = 0
- **Implementation**: Explicit check for identical trees returns 0
- **Verification**: Test with identical trees of various sizes

### 3. Symmetry
- [ ] **Test**: For any two trees T1 and T2, verify that HierarchicalMutualInfoDist(T1, T2) = HierarchicalMutualInfoDist(T2, T1)
- **Implementation**: Algorithm treats both trees symmetrically
- **Verification**: Test with multiple tree pairs both ways

## Hierarchical Information Theory Properties

### 4. Information Content Calculation
- [ ] **Test**: Verify that individual split information is calculated correctly using Shannon entropy
- **Formula**: For a split with n_in tips on one side and n_out on the other: H = -(p_in * log2(p_in) + p_out * log2(p_out))
- **Verification**: Manual calculation for simple splits

### 5. Hierarchical Weighting
- [ ] **Test**: Verify that more balanced splits receive higher hierarchical weights
- **Implementation**: Weight = (min(n_in, n_out) / max(n_in, n_out)) * log2(n_total)
- **Verification**: Compare weights for different split balances

### 6. Split Compatibility Detection
- [ ] **Test**: Verify that compatible splits are correctly identified
- **Cases to test**:
  - Identical splits (compatibility = 1.0)
  - Complementary splits (compatibility = 1.0) 
  - Conflicting splits (compatibility = 0.0)
  - Partially overlapping splits (0 < compatibility < 1)

## Distance Properties Specific to Phylogenetic Trees

### 7. Tree Shape Sensitivity
- [ ] **Test**: Verify that different tree shapes produce different distances
- **Test cases**:
  - Balanced vs Pectinate trees (should have substantial distance)
  - Star tree vs Balanced tree (should have substantial distance)
  - Similar shapes (should have smaller distances)

### 8. Tree Size Scaling
- [ ] **Test**: Verify that distances scale appropriately with tree size
- **Expected behavior**: Larger trees should generally have larger distances for similar shape differences

### 9. Normalization Properties
- [ ] **Test**: When normalize=TRUE, verify that 0 ≤ normalized_distance ≤ 1
- **Implementation**: Normalization divides by total hierarchical information content

## Manual Verification Examples

### 10. Simple 4-tip Trees
Create manual calculations for simple cases:
- [ ] Tree1: ((A,B),(C,D)) vs Tree2: ((A,C),(B,D))
- [ ] Expected: Non-zero distance (different topologies)
- [ ] Manual calculation of split information and compatibility

### 11. 6-tip Tree Examples  
- [ ] Balanced tree: (((A,B),(C,D)),(E,F))
- [ ] Pectinate tree: (A,(B,(C,(D,(E,F)))))
- [ ] Star tree: (A,B,C,D,E,F)
- [ ] Manual verification of relative distances

## Algorithmic Correctness

### 12. Split Processing
- [ ] **Test**: Verify that splits are correctly extracted from phylo objects
- [ ] **Test**: Verify that raw split data is correctly converted to logical vectors
- [ ] **Test**: Verify that split information content is calculated correctly

### 13. Shared Information Calculation
- [ ] **Test**: Verify that shared hierarchical information between trees is calculated correctly
- [ ] **Method**: For each split in tree1, find best matching split in tree2
- [ ] **Verification**: Manual calculation for simple tree pairs

### 14. Edge Cases
- [ ] **Test**: Trees with different tip labels (should reduce to common tips)
- [ ] **Test**: Trees with insufficient tips (< 3) (should return 0)
- [ ] **Test**: Star trees (no internal structure)
- [ ] **Test**: Identical trees with different tip orderings

## Implementation Verification Commands

```r
# Load required libraries
library(TreeTools)
library(ape)
source("R/hierarchical_mutual_info.R")

# Test 1: Non-negativity
tree1 <- BalancedTree(6)
tree2 <- PectinateTree(6)
dist <- HierarchicalMutualInfoDist(tree1, tree2)
stopifnot(dist >= 0)

# Test 2: Identity  
self_dist <- HierarchicalMutualInfoDist(tree1, tree1)
stopifnot(abs(self_dist) < 1e-10)

# Test 3: Symmetry
dist_rev <- HierarchicalMutualInfoDist(tree2, tree1)
stopifnot(abs(dist - dist_rev) < 1e-10)

# Test 4: Normalization
norm_dist <- HierarchicalMutualInfoDist(tree1, tree2, normalize = TRUE)
stopifnot(norm_dist >= 0 && norm_dist <= 1)

# Test 5: Tree shape sensitivity
star6 <- StarTree(6)
star_dist <- HierarchicalMutualInfoDist(tree1, star6)
stopifnot(star_dist > 0)

print("All mathematical correctness tests passed!")
```

## References

The implementation is based on:
- Perotti, J. I., Tessone, C. J., and Caldarelli, G. (2015). Hierarchical mutual information for the comparison of hierarchical community structures in complex networks. Physical Review E, 92(6), 062825.

## Notes

This implementation provides a simplified version of the hierarchical mutual information algorithm adapted for phylogenetic trees. The core mathematical principles are preserved while adapting the method to work with tree topology represented as splits.