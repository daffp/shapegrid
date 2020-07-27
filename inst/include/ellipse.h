/* See https://people.sc.fsu.edu/~jburkardt/cpp_src/ellipse_grid/ellipse_grid.html */
#ifndef ELLIPSE_H
#define ELLIPSE_H

# include <cmath>

/******************************************************************************/

int i4_ceiling ( double x )
  
  /******************************************************************************/
  /*
   Purpose:
   
   I4_CEILING rounds an R8 up to the nearest I4.
   
   Discussion:
   
   The "ceiling" of X is the value of X rounded towards plus infinity.
   
   Example:
   
   X        I4_CEILING(X)
   
   -1.1      -1
   -1.0      -1
   -0.9       0
   -0.1       0
   0.0       0
   0.1       1
   0.9       1
   1.0       1
   1.1       2
   2.9       3
   3.0       3
   3.14159   4
   
   Licensing:
   
   This code is distributed under the GNU LGPL license.
   
   Modified:
   
   10 November 2011
   
   Author:
   
   John Burkardt
   
   Parameters:
   
   Input, double X, the number whose ceiling is desired.
   
   Output, int I4_CEILING, the ceiling of X.
   */
{
  int value;
  
  value = ( int ) x;
  
  if ( value < x )
  {
    value = value + 1;
  }
  
  return value;
}

//****************************************************************************80

int ellipse_grid_count ( int n, std::vector<double> r, std::vector<double> c )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPSE_GRID_COUNT counts the grid points inside an ellipse.
//
//  Discussion:
//
//    The ellipse is specified as
//
//      ( ( X - C1 ) / R1 )^2 + ( ( Y - C2 ) / R2 )^2 = 1
//
//    The user supplies a number N.  There will be N+1 grid points along
//    the shorter axis.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of subintervals.
//
//    Input, double R[2], the half axis lengths.
//
//    Input, double C[2], the center of the ellipse.
//
//    Output, int ELLIPSE_GRID)_COUNT, the number of grid points inside 
//    the ellipse.
//
//   Minor changes  by dp:2020: C style array to std::vector

{
  double h;
  int i;
  int j;
  int nj;
  int p;
  double x;
  double y;

  if ( r[0] < r[1] )
  {
    h = 2.0 * r[0] / ( double ) ( 2 * n + 1 );
    nj = i4_ceiling ( r[1] / r[0] ) * ( double ) ( n );
  }
  else
  {
    h = 2.0 * r[1] / ( double ) ( 2 * n + 1 );
    nj = n;
  }

  p = 0;

  for ( j = 0; j <= nj; j++ )
  {
    i = 0;
    x = c[0];
    y = c[1] + ( double ) ( j ) * h;

    p = p + 1;

    if ( 0 < j )
    {
      p = p + 1;
    }

    for ( ; ; )
    {
      i = i + 1;
      x = c[0] + ( double ) ( i ) * h;

      if ( 1.0 < pow ( ( x - c[0] ) / r[0], 2 ) 
               + pow ( ( y - c[1] ) / r[1], 2 ) )
      {
        break;
      }

      p = p + 1;
      p = p + 1;

      if ( 0 < j )
      {
        p = p + 1;
        p = p + 1;
      }
    }
  }

  return p;
}

/******************************************************************************/

double *ellipse_grid ( int n, std::vector<double> r, std::vector<double> c, int ng )
   
/******************************************************************************/
  /*
   Purpose:
   
   ELLIPSE_GRID generates grid points inside an ellipse.
   
   Discussion:
   
   The ellipse is specified as
   
   ( ( X - C1 ) / R1 )^2 + ( ( Y - C2 ) / R2 )^2 = 1
   
   The user supplies a number N.  There will be N+1 grid points along
   the shorter axis.
   
   Licensing:
   
   This code is distributed under the GNU LGPL license.
   
   Modified:
   
   12 November 2011
   
   Author:
   
   John Burkardt
   
   Parameters:
   
   Input, int N, the number of subintervals.
   
   Input, double R[2], the half axis lengths.
   
   Input, double C[2], the center of the ellipse.
   
   Input, int NG, the number of grid points inside the ellipse.
   
   Output, double ELLIPSE_GRID[2*NG], the grid points.
   
   Minor changes  by dp:2020: C style array to std::vector
   */
{
   double h;
   int i;
   int j;
   int nj;
   int p;
   double x;
   double *xy;
   double y;
   
   xy = new double[2*ng];
   
   if ( r[0] < r[1] )
   {
      h = 2.0 * r[0] / ( double ) ( 2 * n + 1 );
      nj = i4_ceiling ( r[1] / r[0] ) * ( double ) ( n );
   }
   else
   {
      h = 2.0 * r[1] / ( double ) ( 2 * n + 1 );
      nj = n;
   }
   
   p = 0;
   
   for ( j = 0; j <= nj; j++ )
   {
      i = 0;
      x = c[0];
      y = c[1] + ( double ) ( j ) * h;
      
      xy[0+p*2] = x;
      xy[1+p*2] = y;
      p = p + 1;
      
      if ( 0 < j )
      {
         xy[0+p*2] = x;
         xy[1+p*2] = 2.0 * c[1] - y;
         p = p + 1;
      }
      
      for ( ; ; )
      {
         i = i + 1;
         x = c[0] + ( double ) ( i ) * h;
         
         if ( 1.0 < pow ( ( x - c[0] ) / r[0], 2 ) 
                 + pow ( ( y - c[1] ) / r[1], 2 ) )
         {
            break;
         }
         
         xy[0+p*2] = x;
         xy[1+p*2] = y;
         p = p + 1;
         xy[0+p*2] = 2.0 * c[0] - x;
         xy[1+p*2] = y;
         p = p + 1;
         
         if ( 0 < j )
         {
            xy[0+p*2] = x;
            xy[1+p*2] = 2.0 * c[1] - y;
            p = p + 1;
            xy[0+p*2] = 2.0 * c[0] - x;
            xy[1+p*2] = 2.0 * c[1] - y;
            p = p + 1;
         }
      }
   }
   return xy;
}

#endif
