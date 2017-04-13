#include <string>
#include <cmath>

#include "Distribution.h"
#include "Point.h"
#include "Random.h"

double uniform_distribution::sample() {

	return a + Urand() * ( b - a ); 
	
}

double linear_distribution::sample() {

  double r1 = Urand(), r2 = Urand();
  double p  = 2.0 * std::fmin( fa, fb ) / ( fa + fb );
  if ( r1 < p ) { return a + r2 * ( b - a ); }
  else {
    if ( fb > fa ) { return a + ( b - a ) * std::sqrt( r2 ); }
    else           { return a + ( b - a ) * ( 1.0 - std::sqrt( r2 )); }
  }
  
}

double cubic_distribution::sample(){
	double mu, y, f;  

  // rejection sampling to find mu (with respect to x axis)
  // 1/3 efficient for problem 6
  do {
    mu = 2.0  * Urand() - 1.0;
    y =  fmax* Urand();
    f = c1*mu*mu*mu + c2*mu*mu + c3*mu  + c4;
  } while ( y > f );
  return mu;
}

double exponential_distribution::sample() { return -std::log( Urand() ) / lambda; }

double normal_distribution::sample() 
{ return mu + sigma * std::sqrt( -2.0 * std::log( Urand() ) ) * std::cos( twopi * Urand() ); }

double HenyeyGreenstein_distribution::sample() {
  if ( a != 0.0 ) {
    return ( 1.0 + a*a - pow( ( 1.0 - a*a )/( 1.0 + a*(2.0*Urand() - 1.0) ), 2 ) ) / ( 2.0 * a ); 
  }
  else {
    return 2.0 * Urand() - 1.0;
  }
}

int meanMultiplicity_distribution::sample() {
  return (int) std::floor( nu + Urand() );
}

TerrellFission_distribution::TerrellFission_distribution( std::string label, double p1, double p2, double p3 ) 
    : distribution(label), nubar(p1), sigma(p2), b(p3) {
  double c  = 0.0;
  double nu = 0.0;
  while ( c < 1.0 - 1.0e-12 ) {
    double a  = ( nu - nubar + 0.5 + b ) / sigma;
    c = 0.5 * ( 1 + erf( a / sqrt(2.0) ) ) ;

    cdf.push_back(c);
    nu += 1.0;
  }
  cdf.push_back(1.0);
}

int TerrellFission_distribution::sample() {
  double r  = Urand();
  double nu;
  for ( int i = 0 ; i < cdf.size() ; i++ ) {
    if ( r < cdf[i] ) {
      nu = (double) i;
      break;
    }
  }
 return nu; 
}

point isotropicDirection_distribution::sample() {
  // sample polar cosine and azimuthal angle uniformly
  double mu  = 2.0 * Urand() - 1.0;
  double azi = twopi * Urand();

  // convert to Cartesian coordinates
  double c = std::sqrt( 1.0 - mu * mu );
  point p;
  p.x = std::cos( azi ) * c;
  p.y = std::sin( azi ) * c;
  p.z = mu;

  return p;
}

point anisotropicDirection_distribution::sample() {
  double mu  = dist_mu->sample(); 
  double azi = twopi * Urand();
  double cos_azi = std::cos(azi);
  double sin_azi = std::sin(azi);

  // rotate the local particle coordinate system aligned along the incident direction
  // to the global problem (x,y,z) coordinate system 
  double sin_t0 = std::sqrt( 1.0 - mu * mu );
  double c = sin_t0 / sin_t;

  point p;
  p.x = axis.x * mu + ( axis.x * axis.z * cos_azi - axis.y * sin_azi ) * c;
  p.y = axis.y * mu + ( axis.y * axis.z * cos_azi + axis.x * sin_azi ) * c;
  p.z = axis.z * mu - cos_azi * sin_t0 * sin_t;
  return p;
}

point independentXYZ_distribution::sample() {
	return point( dist_x->sample(), dist_y->sample(), dist_z->sample() );
}

  
point shell_distribution::sample() 
{
	double azi  = twopi * Urand();
	double mu	= 2 * (Urand()) - 1.0;
	double rad 	= pow((std::pow(innerRadius,3) + (std::pow(outerRadius,3) - std::pow(innerRadius,3)) * Urand()),(1.0/3.0));
	
  // convert to Cartesian coordinates
  double c = std::sqrt( 1.0 - mu * mu );
  point p;
  p.x = rad * std::cos( azi ) * c;
  p.y = rad * std::sin( azi ) * c;
  p.z = rad * mu;
  
	return p;
}

point zdisk_distribution::sample(){
	double h 	= radius * std::sqrt(Urand());
	double azi  = twopi * Urand();
	double dx 	= h * std::cos(azi);
	double dy 	= h * std::sin(azi);
	point p(disk_center.x + dx, disk_center.y + dy, disk_center.z);
	return p;
	
}

point xdisk_distribution::sample(){
	double h 	= radius * std::sqrt(Urand());
	double azi  = twopi * Urand();
	double dz 	= h * std::cos(azi);
	double dy 	= h * std::sin(azi);
	point p(disk_center.x, disk_center.y + dy, disk_center.z + dz);
	return p;
}
