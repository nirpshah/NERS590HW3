#ifndef _PARTICLE_HEADER_
#define _PARTICLE_HEADER_

#include <memory>

#include "Point.h"
#include "Nuclide.h"

class cell;
class nuclide;

class particle {
  private:
    point p_pos, p_dir;
    double p_energy; // units MeV
    double p_time;   // units ns
    double p_wgt;
    bool   exist;
    std::shared_ptr< cell > p_cell;
//    std::shared_ptr< cell > n_cell;

  public:
     particle( point p, point d, double e, double t );
    ~particle() {};

    point   pos()     { return p_pos; };    // return particle position
    point   dir()     { return p_dir; };    // return particle direction 
    double  energy()  { return p_energy; };
    double  time()    { return p_time; };
    double  wgt()     { return p_wgt; };   // return particle weight
    bool    alive()   { return exist; };   // return particle state flag
    double  speed();

    ray     getRay()  { return ray( p_pos, p_dir ); }

    std::shared_ptr< cell > cellPointer() { return p_cell; }

    void move( double s );
    void scatter( double mu0, std::shared_ptr< nuclide> N );
    void kill();
    void setDirection( point p );
    void adjustWeight( double f );
    void recordCell( std::shared_ptr< cell > cel );
//    void getNextCell( std::shared_ptr< cell > ncel );
};

#endif
