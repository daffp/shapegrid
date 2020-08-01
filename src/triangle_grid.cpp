# include <Rcpp.h>
# include "triangle.h"

using namespace Rcpp;

//' @useDynLib shapegrid
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
NumericVector triangle_gridcpp(int n, Rcpp::NumericVector V){
  
  std::vector<double> t = as<std::vector<double>>(V); 

  double* grid = triangle_grid(n, t);
  
  std::vector<double> values(grid, grid + ((n+1)*(n+2)));
  
  return wrap(values);
}
