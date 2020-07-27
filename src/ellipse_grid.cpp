# include <Rcpp.h>
# include "ellipse.h"

using namespace Rcpp;

//' @useDynLib shapegrid
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
NumericVector ellipse_gridcpp(int n, Rcpp::NumericVector R, Rcpp::NumericVector  C){
  
  std::vector<double> r = as<std::vector<double>>(R); 
  std::vector<double> c = as<std::vector<double>>(C);
  
  int ng = ellipse_grid_count(n, r, c);
  double* grid = ellipse_grid(n, r, c, ng);
  
  std::vector<double> values(grid, grid + 2*ng);

  return wrap(values);
}
