# include <Rcpp.h>
# include "polygon.h"

using namespace Rcpp;

//' @useDynLib shapegrid
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
SEXP polygon_gridcpp(int n, int nv, Rcpp::NumericVector V){
  
  std::vector<double> v = as<std::vector<double>>(V); 
  
  
  int ng = polygon_grid_count(n, nv);
  double* grid = polygon_grid_points(n, nv, v, ng);
  
  std::vector<double> xg(grid, grid + 2*ng);
  
  // from polygon_grid_display 
  NumericMatrix out(ng, 2);
  for ( int j = 0; j < ng; j++ ){
    out(j, 0) = xg[0 + j* 2];
    out(j, 1) = xg[1 + j* 2];
  }

  return wrap(out);

}
