#' Triangle grid
#' 
#' Create a grid of points within a triangle
#'
#' @param N The number of subintervals each side should be divided into i.e. there are (N+1) points on a side.
#' @param vertices A vector of length six giving the (x,y) coordinates of each vertex of the triangle or a matrix with two rows and three columns of vertex coordinates.
#' @return Returns a matrix with the x and y positions of each point in the grid.
#' @details Code taken from \url{https://people.sc.fsu.edu/~jburkardt/cpp_src/triangle_grid/triangle_grid.html}. 
#' @examples
#' v = cbind(c(0,0), c(1,1), c(2,0))
#' p = shapegrid::triangle_grid(10, v)
#' plot(p)
#' @export
triangle_grid = function(N, vertices){
  
  v = c(vertices)
  pts = triangle_gridcpp(N, v)
  mat = matrix(pts, ncol=2, byrow=TRUE)
  
  return(mat)
}
