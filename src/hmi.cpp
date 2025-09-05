#include <Rcpp.h>
#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include <set>
#include <map>
#include <algorithm>

using namespace Rcpp;

// ------------------------
// Partition structure
// ------------------------
struct Partition {
  std::string label;
  std::vector<std::shared_ptr<Partition>> children;
  
  Partition() : label("") {}
  Partition(const std::string &l) : label(l) {}
  bool is_leaf() const { return children.empty(); }
};
using PartitionPtr = std::shared_ptr<Partition>;

// ------------------------
// Flatten leaves of a partition
// ------------------------
void flattenator(const PartitionPtr &node, std::set<std::string> &out) {
  if (node->is_leaf()) { 
    out.insert(node->label); 
    return; 
  }
  for (auto &child : node->children) {
    if (child) flattenator(child, out);
  }
}

// ------------------------
// x * log(x), with 0*log(0)=0
// ------------------------
inline double xlnx(double x) {
  return (x > 0.0) ? x * std::log(x) : 0.0;
}

// ------------------------
// Traceable hierarchical mutual information
// ------------------------
std::pair<double,double> compute_HMI(const PartitionPtr &Ut, const PartitionPtr &Us, int depth=0) {
  std::string indent(depth*2,' ');
  
  if (Ut->is_leaf() && Us->is_leaf()) {
    double n_ts = (Ut->label == Us->label) ? 1.0 : 0.0;
    Rcout << indent << "Both leaves: Ut=" << Ut->label 
          << ", Us=" << Us->label 
          << ", n_ts=" << n_ts << std::endl;
    return {n_ts, 0.0};
  }
  
  std::vector<PartitionPtr> u_children = Ut->is_leaf() ? std::vector<PartitionPtr>{Ut} : Ut->children;
  std::vector<PartitionPtr> v_children = Us->is_leaf() ? std::vector<PartitionPtr>{Us} : Us->children;
  
  std::set<std::string> leaves_Ut, leaves_Us;
  flattenator(Ut, leaves_Ut);
  flattenator(Us, leaves_Us);
  
  std::vector<std::string> intersection;
  std::set_intersection(
    leaves_Ut.begin(), leaves_Ut.end(),
    leaves_Us.begin(), leaves_Us.end(),
    std::back_inserter(intersection)
  );
  double n_ts = static_cast<double>(intersection.size());
  if (n_ts == 0.0) {
    Rcout << indent << "n_ts=0 for Ut/Us leaves intersection" << std::endl;
    return {0.0, 0.0};
  }
  Rcout << indent << "Ut internal, Us internal: n_ts=" << n_ts << std::endl;
  
  double H_uv = 0.0, H_u = 0.0, H_v = 0.0, mean_I = 0.0;
  std::vector<double> n_v(v_children.size(),0.0);
  
  for (size_t i=0; i<u_children.size(); ++i) {
    double n_u = 0.0;
    for (size_t j=0; j<v_children.size(); ++j) {
      auto res = compute_HMI(u_children[i], v_children[j], depth+1);
      double n_uv = res.first;
      double I_uv = res.second;
      
      Rcout << indent << "Child pair: i=" << i << ", j=" << j 
            << ", n_uv=" << n_uv << ", I_uv=" << I_uv << std::endl;
      
      H_uv += xlnx(n_uv);
      n_u += n_uv;
      n_v[j] += n_uv;
      mean_I += n_uv*I_uv;
    }
    H_u += xlnx(n_u);
  }
  for (double nv : n_v) H_v += xlnx(nv);
  
  double local_I = std::log(n_ts) - (H_u + H_v - H_uv)/n_ts;
  double I_ts = local_I + mean_I/n_ts;
  
  Rcout << indent << "Returning: n_ts=" << n_ts << ", I_ts=" << I_ts << std::endl;
  return {n_ts, I_ts};
}

// ------------------------
// Hierarchical entropy, variation of information, distance
// ------------------------
double HH(const PartitionPtr &hp) {
  return compute_HMI(hp,hp).second;
}

double HVI(const PartitionPtr &hp1, const PartitionPtr &hp2) {
  return HH(hp1) + HH(hp2) - 2.0 * compute_HMI(hp1,hp2).second;
}

double d_n(const PartitionPtr &t1, const PartitionPtr &t2, int n=1) {
  double ln2d2 = 0.5 * std::log(2.0);
  return 1.0 - std::exp(-n*ln2d2*HVI(t1,t2));
}

// ------------------------
// Convert R nested list to Partition
// ------------------------
PartitionPtr from_sexp(SEXP x) {
  if (Rf_isString(x)) {
    return std::make_shared<Partition>(as<std::string>(STRING_ELT(x,0)));
  } else if (Rf_isVectorList(x)) {
    PartitionPtr node = std::make_shared<Partition>();
    List lst(x);
    for (int i=0; i<lst.size(); ++i) {
      node->children.push_back(from_sexp(lst[i]));
    }
    return node;
  } else {
    stop("Unsupported type in from_sexp");
  }
  return nullptr;
}

// ------------------------
// Rcpp exports
// ------------------------
// [[Rcpp::export]]
double HMI_nested(SEXP t1, SEXP t2) {
  PartitionPtr p1 = from_sexp(t1);
  PartitionPtr p2 = from_sexp(t2);
  return compute_HMI(p1,p2).second;
}

// [[Rcpp::export]]
double d_n_nested(SEXP t1, SEXP t2, int n=1) {
  PartitionPtr p1 = from_sexp(t1);
  PartitionPtr p2 = from_sexp(t2);
  return d_n(p1,p2,n);
}
