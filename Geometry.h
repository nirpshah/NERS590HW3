#ifndef _GEOMETRY_HEADER_
#define _GEOMETRY_HEADER_

#include <vector>
#include <memory>

#include "Cell.h"
#include "Particle.h"


class geometry
{
  private:
    std::vector<std::shared_ptr<cell>> cells;
    
  public:
    
    geometry(std::vector<std::shared_ptr<cell>> inCells);
    std::shared_ptr<cell> findCell(particle* p );
    void particleSplitRoulette(particle* p, std::stack<particle>* bank, double Inew);
};

#endif
