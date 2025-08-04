#include <Rcpp/Lightest>
#include <TreeTools/assert.h> // for ASSERT
#include <TreeTools/types.h> // for intx
#include <TreeTools/keep_tip.h> // for keep_tip
#include <TreeTools/root_tree.h> // for root_on_node
#include <TreeTools/SplitList.h> // for SplitList, count_bits
#include "reduce_tree.h"

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector mismatch_size (const RawMatrix x, const RawMatrix y) {
  if (double(x.rows()) > double(std::numeric_limits<int16>::max())) {
    Rcpp::stop("This many splits are not (yet) supported.");
  }
  const int16 n_split = int16(x.rows());
  if (n_split != y.rows()) {
    throw std::invalid_argument("`x` and `y` differ in number of splits.");
  }
  if (!x.hasAttribute("nTip")) {
    Rcpp::stop("`x` lacks nTip attribute");
  }
  if (!y.hasAttribute("nTip")) {
    Rcpp::stop("`y` lacks nTip attribute");
  }
  const int16 n_tip = x.attr("nTip");
  if (n_tip != int16(y.attr("nTip"))) {
    Rcpp::stop("`x` and `y` differ in `nTip`");
  }
  
  const TreeTools::SplitList a(x), b(y);
  const int16
    half_tip = n_tip / 2,
    last_bin = a.n_bins - 1,
    unset_tips = (n_tip % SL_BIN_SIZE) ? SL_BIN_SIZE - n_tip % SL_BIN_SIZE : 0
  ;
  constexpr splitbit all_ones = ~(splitbit(0U));
  const splitbit unset_mask = all_ones >> unset_tips;

  IntegerVector ret(n_split * n_split);
  int *ret_ptr = ret.end();
  for (int16 bi = b.n_splits; bi--; ) {
    // Rcout << "a = " << ai << ".\n";
    for (int16 ai = a.n_splits; ai--; ) {
      // Rcout << "  - b = " << bi << ".\n";
      
      --ret_ptr;
      
      // Rcout << "    - last_bin: " << ((a.state[ai][last_bin] ^ b.state[bi][last_bin]) & unset_mask)
      //       << " = " << TreeTools::count_bits(
      // (a.state[ai][last_bin] ^ b.state[bi][last_bin]) & unset_mask
      //       ) << "\n";
      *ret_ptr = TreeTools::count_bits(
        (a.state[ai][last_bin] ^ b.state[bi][last_bin]) & unset_mask
        );
      for (int16 bin = last_bin; bin--; ) {
        // Rcout << "    - bin = " << bin << ".\n";
        // Rcout << "      " << (a.state[ai][bin]);
        // Rcout << " ^ " << (b.state[bi][bin]);
        // Rcout << " = " << (a.state[ai][bin] ^ b.state[bi][bin]) << std::endl;
        *ret_ptr += TreeTools::count_bits(a.state[ai][bin] ^ b.state[bi][bin]);
        // Rcout << "      ret[" << ret_ptr << "] = " 
        //       << TreeTools::count_bits(a.state[ai][bin] ^ b.state[bi][bin]) 
        //       <<".\n";
      }
      if (*ret_ptr > half_tip) {
        *ret_ptr = n_tip - *ret_ptr;
      }
    }
  }
  return ret;
}

// [[Rcpp::export]]
IntegerVector confusion (const RawMatrix x, const RawMatrix y) {
  if (double(x.rows()) > double(std::numeric_limits<int16>::max())) {
    Rcpp::stop("This many splits are not (yet) supported.");
  }
  const int16 n_split = int16(x.rows());
  if (n_split != y.rows()) {
    throw std::invalid_argument("Input splits contain same number of splits.");
  }
  if (!x.hasAttribute("nTip")) {
    Rcpp::stop("`x` lacks nTip attribute");
  }
  if (!y.hasAttribute("nTip")) {
    Rcpp::stop("`y` lacks nTip attribute");
  }
  const int16 n_tip = x.attr("nTip");
  if (n_tip != int16(y.attr("nTip"))) {
    Rcpp::stop("`x` and `y` differ in `nTip`");
  }
  
  const TreeTools::SplitList a(x), b(y);
  const int16
    n_bin = a.n_bins,
    confusion_size = 4
  ;
  IntegerVector ret(n_split * n_split * confusion_size);
  int *ret_ptr = ret.end();
  for (int16 bi = n_split; bi--; ) {
    const int16
      nb = b.in_split[bi],
      nB = n_tip - nb
    ;
    
    for (int16 ai = n_split; ai--; ) {
      
      // x divides tips into a|A; y divides tips into b|B
      int16 a_and_b = 0;
      for (int16 bin = n_bin; bin--; ) {
        a_and_b += TreeTools::count_bits(a.state[ai][bin] & b.state[bi][bin]);
      }
      
      const int16
        na = a.in_split[ai],
        a_and_B = na - a_and_b,
        A_and_b = nb - a_and_b,
        A_and_B = nB - a_and_B
      ;
      *(--ret_ptr) = A_and_B;
      *(--ret_ptr) = A_and_b;
      *(--ret_ptr) = a_and_B;
      *(--ret_ptr) = a_and_b;
    }
  }
  ret.attr("dim") = Dimension(confusion_size, n_split, n_split);
  return ret;
}

IntegerMatrix reverse (const IntegerMatrix x) {
  if (double(x.nrow()) > double(std::numeric_limits<intx>::max())) {
    Rcpp::stop("This many edges are not (yet) supported.");
  }
  const intx n_edge = intx(x.nrow());
  ASSERT(n_edge % 2 == 0); // Tree is binary
  IntegerMatrix ret(n_edge, 2);
  
  for (intx i = n_edge / 2; i--; ) {
    const intx j = n_edge - i - 1;
    ret[i] = x[j];
    ret[i + n_edge] = x[j + n_edge];
    ret[j] = x[i];
    ret[j + n_edge] = x[i + n_edge];
  }
  return ret;
}

// tree1 and tree2 are binary trees in postorder with identical tip.labels
// [[Rcpp::export]]
List keep_and_reroot(const List tree1,
                     const List tree2,
                     const LogicalVector keep) {
  IntegerMatrix
    postorder1 = tree1["edge"],
    postorder2 = tree2["edge"]
  ;
  
  ASSERT(postorder1.nrow() % 2 == 0); // Tree is binary
  ASSERT(postorder2.nrow() % 2 == 0); // Tree is binary
  
  IntegerMatrix
    pre1 = reverse(postorder1),
    pre2 = reverse(postorder2)
  ;
  
  ASSERT(postorder1.nrow() / 2 + 1 == keep.length());
  // Rcout << "\n \n === Keep & Reroot ===\n";
  // Rcout << " Keeping: ";
  // for (int i = 0; i != keep.size(); i++) Rcout << (keep[i] ? "*" : ".");
  IntegerMatrix ret_edge1 = TreeTools::keep_tip(pre1, keep);
  IntegerMatrix ret_edge2 = TreeTools::keep_tip(pre2, keep);
  
  const intx n_node = ret_edge1.nrow() / 2;
  if (!n_node) {
    List nullTree = List::create(Named("edge") = ret_edge1,
                                 _["Nnode"] = n_node,
                                 _["tip.label"] = CharacterVector(0));
    
    nullTree.attr("class") = "phylo";
    nullTree.attr("order") = "preorder";
    return List::create(nullTree, nullTree);
  }
  
  const intx n_tip = n_node + 1;
  CharacterVector
    old_label = tree1["tip.label"],
    new_labels(n_tip)
  ;
  
  // Rcout << ret_edge1.nrow() << " rows; Kept " << n_tip << " tips and "
  //       << n_node << " nodes.\n";
  
  if (old_label.size() > std::numeric_limits<int16>::max()) {
    Rcpp::stop("This many leaves are not (yet) supported.");
  }
  intx next_tip = n_tip;
  for (intx i = intx(old_label.size()); i--; ) {
    if (keep[i]) {
      --next_tip;
      new_labels[next_tip] = old_label[i];
    }
  }
  
  List ret1 = List::create(Named("edge") = ret_edge1,
                           _["Nnode"] = n_node,
                           _["tip.label"] = new_labels),
       ret2 = List::create(Named("edge") = ret_edge2,
                           _["Nnode"] = n_node,
                           _["tip.label"] = new_labels);
  ret1.attr("class") = "phylo";
  ret1.attr("order") = "preorder";
  ret2.attr("class") = "phylo";
  ret2.attr("order") = "preorder";
  return List::create(
    TreeTools::root_on_node(ret1, 1),
    TreeTools::root_on_node(ret2, 1)
  );
}

// [[Rcpp::export]]
List keep_and_reduce(const List tree1,
                     const List tree2,
                     const LogicalVector keep) {
  List rerooted = keep_and_reroot(tree1, tree2, keep);
  
  List rerooted1 = rerooted[0];
  List rerooted2 = rerooted[1];
  IntegerMatrix edge1 = reverse(rerooted1["edge"]);
  IntegerMatrix edge2 = reverse(rerooted2["edge"]);
  
  if (edge1.nrow() < 1) {
    List nullTree = List::create(Named("edge") = NumericMatrix(0, 2),
                                 _["Nnode"] = 0,
                                 _["tip.label"] = CharacterVector(0));
    
    nullTree.attr("class") = "phylo";
    nullTree.attr("order") = "preorder";
    return List::create(nullTree, nullTree);
  }
  
  CharacterVector tip_label = rerooted1["tip.label"];
  return reduce_trees(edge1, edge2, tip_label);
}
