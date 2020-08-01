/* See https://people.sc.fsu.edu/~jburkardt/cpp_src/polygon_grid/polygon_grid.html*/
#ifndef POLYGON_H
#define POLYGON_H
  
# include <cmath>  
  
using namespace std;
  

//****************************************************************************80

int polygon_grid_count ( int n, int nv )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_GRID_COUNT counts the grid points inside a polygon.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of subintervals on a side.
//
//    Input, int NV, the number of vertices.
//    3 <= NV.
//
//    Output, int POLYGON_GRID_COUNT, the number of grid points.
//
{
int ng;

ng = 1 + nv * ( n * ( n + 1 ) ) / 2;

return ng;
}  

//****************************************************************************80

double *polygon_grid_points ( int n, int nv, std::vector<double> v, int ng )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_GRID_POINTS computes points on a polygonal grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of subintervals.
//
//    Input, int NV, the number of vertices in the polygon.
//
//    Input, double V[2*NV], the coordinates of the vertices.
//
//    Input, int NG, the number of grid points.
//
//    Output, double POLYGON_GRID_POINTS[2*NG], the coordinates of the 
//    grid points.
//
//
//   Minor changes  by dp:2020: C style array to std::vector
//
{
  int i;
  int j;
  int k;
  int l;
  int lp1;
  int p;
  double vc[2];
  double *xg;
  
    xg = new double[2*ng];
    p = 0;
  //
  //  Determine the centroid.
  //
  vc[0] = 0.0;
  vc[1] = 0.0;
  for ( j = 0; j < nv; j++ )
  {
    vc[0] = vc[0] + v[0+j*2];
    vc[1] = vc[1] + v[1+j*2];
  }
  vc[0] = vc[0] / ( double ) ( nv );
  vc[1] = vc[1] / ( double ) ( nv );
  //
  //  The centroid is the first point.
  //
  xg[0+p*2] = vc[0];
  xg[1+p*2] = vc[1];
  p = p + 1;
  //
  //  Consider each triangle formed by two consecutive vertices and the centroid,
  //  but skip the first line of points.
  //
  for ( l = 0; l < nv; l++ )
  {
    lp1 = ( ( l + 1 ) % nv );
    for ( i = 1; i <= n; i++ )
    {
      for ( j = 0; j <= n - i; j++ )
      {
        k = n - i - j;
        xg[0+p*2] = ( ( double ) ( i ) * v[0+l*2]   
                    + ( double ) ( j ) * v[0+lp1*2] 
                    + ( double ) ( k ) * vc[0] )  
                    / ( double ) ( n );
        xg[1+p*2] = ( ( double ) ( i ) * v[1+l*2]   
                    + ( double ) ( j ) * v[1+lp1*2] 
                    + ( double ) ( k ) * vc[1] )  
                    / ( double ) ( n );
        p = p + 1;
      }
    }
  }
  
  return xg;
}

#endif
