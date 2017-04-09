#ifndef _SOURCE_HEADER_
#define _SOURCE_HEADER_

#include <stack>
#include <memory>

#include "Point.h"
#include "Distribution.h"
#include "Particle.h"

class source {
  private:
    std::shared_ptr< distribution<point> >  dist_pos;
    std::shared_ptr< distribution<point> >  dist_dir;
    std::shared_ptr< distribution<double> > dist_energy;
    std::shared_ptr< distribution<double> > dist_time;
    
  public:
     source(  std::shared_ptr< distribution<point> > pos, 
              std::shared_ptr< distribution<point> > dir, 
              std::shared_ptr< distribution<double> > energy, 
              std::shared_ptr< distribution<double> > time)
            : dist_pos(pos), 
              dist_dir(dir), 
              dist_energy(energy), 
              dist_time(time) {};
    ~source() {};

    std::stack<particle> sample();
};

#endif
