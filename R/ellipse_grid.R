#' Ellipse grid
#' 
#' Create a grid of points within an ellipse.
#'
#' @param N The number of subintervals to divide the smaller radius into.
#' @param radius A vector of length two giving the radius of each length.
#' @param center A vector of length two giving the center of the ellipse.
#' @return Returns a matrix with the x and y positions of each point in the grid.
#' @details Code taken from \url{https://people.sc.fsu.edu/~jburkardt/cpp_src/ellipse_grid/ellipse_grid.html}.
#' @examples
#' p = shapegrid::ellipse_grid(25, c(1,1), c(0,0))
#' plot(p)
#' @export
ellipse_grid = function(N, radius, center=c(0,0)){
  
   pts = ellipse_gridcpp(N, radius, center)
   mat = matrix(pts, ncol=2, byrow=TRUE)
   
   return(mat)
}
    


