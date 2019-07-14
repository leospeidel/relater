#include <RcppCommon.h>

#include "anc.h"

#include <Rcpp.h>
using namespace Rcpp;

// wrap Node struct
namespace Rcpp {
template <>
SEXP wrap(const Node& x);
}

namespace Rcpp {
template <>
SEXP wrap(const Node& x) {
  return Rcpp::wrap(Rcpp::List::create(Rcpp::Named("label")         = x.label,
                                       Rcpp::Named("branch_length") = x.branch_length,
                                       Rcpp::Named("num_events")    = x.num_events,
                                       Rcpp::Named("SNP_begin")     = x.SNP_begin,
                                       Rcpp::Named("SNP_end")       = x.SNP_end
                                       )
                    );
};
}


// [[Rcpp::export]]
double test_node() {
  Node n;
  n.branch_length = 2;
  return n.branch_length;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
test_node()
*/
