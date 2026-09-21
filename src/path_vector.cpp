#include <Rcpp/Lightest>
#include <TreeTools/renumber_tree.h> // for postorder_order
#include <memory> // for make_unique
using namespace Rcpp;

#define PO_PARENT(i) edge(postorder[i] - 1, 0)
#define PO_CHILD(i) edge(postorder[i] - 1, 1)

// The path length between leaves a and b is
//   depth(a) + depth(b) - 2 * depth(LCA(a, b)).
// Each pair is written exactly once, at the node where the leaf sets of two of
// its children are first united, so the work is Theta(n^2) whatever the tree
// shape.  (A per-pair walk down shared ancestry, as used previously, costs
// O(n^2 * depth): its data-dependent inner loop dominated the run time and
// made that time sensitive to where the linker happened to place it.)
// [[Rcpp::export]]
IntegerVector path_vector(IntegerMatrix edge) {

  IntegerVector postorder = TreeTools::postorder_order(edge);

  const int n_edge = edge.nrow();
  const int n_vert = n_edge + 1;
  const int root_node = PO_PARENT(n_edge - 1);
  const int n_tip = root_node - 1;

  // Node depths (root = 0), by preorder traversal.
  auto depth = std::make_unique<int[]>(n_vert + 1);
  for (int i = n_edge; i--; ) {
    depth[PO_CHILD(i)] = depth[PO_PARENT(i)] + 1;
  }

  // Leaves below each node, as singly-linked lists threaded through `next`;
  // 0 terminates a list (leaves are numbered from 1).
  auto head = std::make_unique<int[]>(n_vert + 1);
  auto tail = std::make_unique<int[]>(n_vert + 1);
  auto next = std::make_unique<int[]>(n_tip + 1);
  for (int tip = 1; tip <= n_tip; ++tip) {
    head[tip] = tip;
    tail[tip] = tip;
  }

  // The pair (a, b), a < b, sits at position row[a] + b of the result: the
  // lower triangle of a distance matrix, read column-wise, as in `dist`.
  auto row = std::make_unique<int[]>(n_tip + 1);
  for (int a = 1; a <= n_tip; ++a) {
    row[a] = (a - 1) * n_tip - a * (a - 1) / 2 - a - 1;
  }

  IntegerVector ret(n_tip * (n_tip - 1) / 2);
  int* const out = ret.begin();

  for (int i = 0; i != n_edge; ++i) { // Postorder traversal
    const int parent = PO_PARENT(i);
    const int child = PO_CHILD(i);
    const int twice_lca = depth[parent] << 1;

    // Pair each leaf already gathered at `parent` with each leaf below `child`
    for (int a = head[parent]; a; a = next[a]) {
      const int depth_a = depth[a] - twice_lca;
      for (int b = head[child]; b; b = next[b]) {
        const int idx = a < b ? row[a] + b : row[b] + a;
        out[idx] = depth_a + depth[b];
      }
    }

    // Append child's leaves to parent's list
    if (head[parent]) {
      next[tail[parent]] = head[child];
    } else {
      head[parent] = head[child];
    }
    tail[parent] = tail[child];
  }

  return ret;
}
