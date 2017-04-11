#ifndef _SIMULATION_HEADER_
#define _SIMULATION_HEADER_

#include "Cell.h"
#include "Particle.h"
#include "Geometry.h"
#include "Source.h"
#include "Distribution.h"

#include <vector>
#include <stack>
#include <memory>


class simulation
{
  private:
    geometry* the_geometry;
    std::vector<std::shared_ptr<estimator>> estimators;
    std::shared_ptr<source> src;
    unsigned long long num_tracks;
    
    void endSimulation();

  public:
    simulation(geometry* in_geometry, std::vector<std::shared_ptr<estimator>> in_estimators, std::shared_ptr< source > in_src) : the_geometry(in_geometry), estimators(in_estimators), src(in_src) { num_tracks = 0;}
    void transport(unsigned long long N_start, unsigned long long number_of_histories);
};

#endif
