#include <cmath>
#include <limits>

#include "Particle.h"
#include "Random.h"

// constructor for a new particle
particle::particle( point p, point d, double e, double t ) : p_pos(p), p_dir(d), p_energy(e), p_current_time(t), p_past_time(0.0) {
  p_dir.normalize();
  exist = true;
  p_wgt = 1.0;
  p_cell = nullptr;
}

double particle::speed()
{
  double v = sqrt(2 * p_energy * 1e6 * 1.60218e-19 / 1.674929e-27) * 100 / 1e9; // cm/ns. sqrt(2E/m) the rest is units
  return(v);
}

// move the particle along its current trajectory and update time
void particle::move( double s ) {
  p_pos.x += s * p_dir.x;
  p_pos.y += s * p_dir.y;
  p_pos.z += s * p_dir.z;
  p_past_time = p_current_time;
  p_current_time  += s / speed();
}

// scatter particle given input direction cosine cos_t0 = mu0
void particle::scatter( double cos_t0, std::shared_ptr< nuclide > N ) {
  
  if (N->A() == 0) // not energy dependent
  {
    // do nothing
  }
  else
  {
    p_energy *= ((1 + N->alpha() + (1 - N->alpha())*cos_t0)/2);
  // must change the direction vector to lab frame	
	cos_t0 = (1+N->A() * cos_t0) / std::pow((1.0 + N->A()* N->A() + 2.0 * N->A() * cos_t0), 0.5);
  }
  
  // sample a random azimuthal angle uniformly
  double azi = 2.0 * std::acos(-1.0) * Urand();
  double cos_azi = std::cos(azi);
  double sin_azi = std::sin(azi);

  // rotate the local particle coordinate system aligned along the incident direction
  // to the global problem (x,y,z) coordinate system 
  double sin_t  = std::sqrt( 1.0 - p_dir.z  * p_dir.z  );
  double sin_t0 = std::sqrt( 1.0 - cos_t0 * cos_t0 );

  point q;

  if ( sin_t > std::numeric_limits<double>::epsilon() * 1000.0 ) {
    double c = sin_t0 / sin_t;
    q.x = p_dir.x * cos_t0 + ( p_dir.x * p_dir.z * cos_azi - p_dir.y * sin_azi ) * c;
    q.y = p_dir.y * cos_t0 + ( p_dir.y * p_dir.z * cos_azi + p_dir.x * sin_azi ) * c;
    q.z = p_dir.z * cos_t0 - cos_azi * sin_t0 * sin_t;
  }
  else {
    // if incident direction along z, reorient axes to avoid division by zero
    sin_t  = std::sqrt( 1.0 -  p_dir.y * p_dir.y );
    double c = sin_t0 / sin_t;
    q.x = p_dir.x * cos_t0 + ( p_dir.x * p_dir.y * cos_azi + p_dir.z * sin_azi ) * c;
    q.y = p_dir.y * cos_t0 - cos_azi * sin_t0 * sin_t;
    q.z = p_dir.z * cos_t0 + ( p_dir.y * p_dir.z * cos_azi - p_dir.x * sin_azi ) * c;
  }

  p_dir.x = q.x;
  p_dir.y = q.y;
  p_dir.z = q.z;
}

// set the particles life flag to false
void particle::kill() {
  exist = false;
}

// set the particle's direction
void particle::setDirection( point p ) {
  p_dir.x = p.x;
  p_dir.y = p.y;
  p_dir.z = p.z;
}

// adjust the weight by a factor f
void particle::adjustWeight( double f ) {
  p_wgt *= f;
}

// set the cell pointer for efficiency
void particle::recordCell( std::shared_ptr< cell > cel ) {
  p_cell = cel;
}

// get the next cell pointer
//void particle::getNextCell( std::shared_ptr< cell > ncel ) {
//  n_cell = ncel;
//}

