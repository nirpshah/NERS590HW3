#include <cmath>
#include <limits>
#include <cassert>

#include "Point.h"
#include "QuadSolver.h"
#include "Surface.h"

// evaluates the surface equation w.r.t. to point p
double plane::eval( point p ) {
  return a * p.x  +  b * p.y  +  c * p.z  - d;
}

// determines the mininum positive distance to intersection for a ray r
// (returns a very large number if no intersection along ray for ease of calculation down the line)
double plane::distance( ray r ) {
  point p = r.pos;
  point u = r.dir;

  double denom = a * u.x  +  b * u.y  +  c * u.z;
  if ( std::fabs( denom ) > 100.0 * std::numeric_limits<double>::epsilon() ) {
    double dist = ( d - a * p.x - b * p.y - c * p.z ) / denom;
    if ( dist > 0.0 ) { return dist; }
    else { return std::numeric_limits<double>::max(); }
  }
  else {
    // moving in a direction that is (or is very close to) parallel to the surface
    return std::numeric_limits<double>::max();
  }
}

// get new reflected direction
point plane::reflect( ray r ) {
   assert( std::fabs( eval( r.pos ) ) < std::numeric_limits<float>::epsilon() );

   point  u = r.dir;
   double t = 2.0 * ( a * u.x  +  b * u.y  +  c * u.z ) / ( a*a + b*b + c*c );
   return point( u.x - a*t, u.y - b*t, u.z - c*t );
}

// for sphere /////////////////////////////////////////////////////////////////////////////

double sphere::eval( point p ) {
  return std::pow( p.x - x0, 2 ) + std::pow( p.y - y0, 2 ) + std::pow( p.z - z0, 2 )  - rad*rad;
}

double sphere::distance( ray r ) {
  point p = r.pos;
  point u = r.dir;

  // difference between each coordinate and current point
  point q( p.x - x0, p.y - y0, p.z - z0 );

  // put into quadratic equation form: a*s^2 + b*s + c = 0, where a = 1
  double b = 2.0 * ( q.x * u.x  +  q.y * u.y  +  q.z * u.z );
  double c = eval( p );

  return quad_solve( 1.0, b, c );
}

point sphere::reflect( ray r ) {
   assert( std::fabs( eval( r.pos ) ) < std::numeric_limits<float>::epsilon() );

   point p = r.pos;
   point u = r.dir;

   point q( p.x - x0, p.y - y0, p.z - z0 );

   double t = 2.0 * ( q.x * u.x  +  q.y * u.y  +  q.z * u.z ) / ( rad*rad );
   return point( u.x - q.x * t,  u.y - q.y * t,  u.z - q.z * t );
}


// for xcylinder /////////////////////////////////////////////////////////////////////////////
double xcylinder::eval( point p ) {
  return std::pow( p.y - y0, 2 ) + std::pow( p.z - z0, 2 )  - rad*rad;
}

double xcylinder::distance( ray r ) {
  point p = r.pos;
  point u = r.dir;

  // difference between each coordinate and current point
  double y_bar   = p.y - y0;
  double z_bar   = p.z - z0;

  // put into quadratic equation form: a*s^2 + b*s + c = 0 
  // where s is distance from point to surface
  
  double a = u.y * u.y + u.z * u.z;
  double b = 2.0 * (y_bar * u.y + z_bar * u.z);
  double c = eval( p );

  return quad_solve( a, b, c );
}

point xcylinder::reflect( ray r ) {
   assert( std::fabs( eval( r.pos ) ) < std::numeric_limits<float>::epsilon() );

   point p = r.pos;
   point u = r.dir;

   point q( 0, p.y - y0, p.z - z0 );

   double t = 2.0 * ( q.x * u.x  +  q.y * u.y  +  q.z * u.z ) / ( rad*rad );
   return point( u.x - q.x * t,  u.y - q.y * t,  u.z - q.z * t );
}

// for zcylinder /////////////////////////////////////////////////////////////////////////////
double zcylinder::eval( point p ) {
  return std::pow( p.x - x0, 2 ) + std::pow( p.y - y0, 2 )  - rad*rad;
}

double zcylinder::distance( ray r ) {
  point p = r.pos;
  point u = r.dir;

  // difference between each coordinate and current point
  double x_bar   = p.x - x0;
  double y_bar   = p.y - y0;

  // put into quadratic equation form: a*s^2 + b*s + c = 0 
  // where s is distance from point to surface
  
  double a = u.x * u.x + u.y * u.y;
  double b = 2.0 * (x_bar * u.x + y_bar * u.y);
  double c = eval( p );

  return quad_solve( a, b, c );
}

point zcylinder::reflect( ray r ) {
   assert( std::fabs( eval( r.pos ) ) < std::numeric_limits<float>::epsilon() );

   point p = r.pos;
   point u = r.dir;

   point q( p.x - x0, p.y - y0,0  );

   double t = 2.0 * ( q.x * u.x  +  q.y * u.y  +  q.z * u.z ) / ( rad*rad );
   return point( u.x - q.x * t,  u.y - q.y * t,  u.z - q.z * t );
}