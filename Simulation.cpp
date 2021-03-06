
#include <iostream>
#include "Simulation.h"


void simulation::transport(unsigned long long N_start, unsigned long long number_of_histories)
{
  // RN_init_particle(&N_start); // I thought this was the time skip but it doesn't work as expected
  double rel_error = 0.0;
  for (unsigned long long history = N_start; history < number_of_histories; history++)
  {
    std::stack<particle> bank = src->sample();
    
    while ( ! bank.empty())
    {
      particle p = bank.top();
      bank.pop();
    
      while (p.alive())
      {
        num_tracks++;
        // find the cell the particle is in. You can use geometry::findCell. Again, not fast, but whatever
		
        p.recordCell(the_geometry->findCell(&p ));
        // Create exp distribution. We should create a function in exp distribution that allows us to update the lambda. This should make the code slightly faster. But then again, eh.

		    exponential_distribution transport_exponential("transport exponential",p.cellPointer()->getMaterial()->macro_xs(&p));
        // sample distance to collision
        
		
		    double distance_to_collision = transport_exponential.sample();
        // find nearest surface 
        std::pair<std::shared_ptr<surface>,double> surface_pair = p.cellPointer()->surfaceIntersect(p.getRay());
        // compare distance to collision to surface intersect
        // if collision first
		
        if (distance_to_collision < surface_pair.second)
        {
		      p.cellPointer()->moveParticle(&p, distance_to_collision); // estimator stuff is done here
          p.cellPointer()->sampleCollision(&p, &bank); // sampling and populating bank is done her
          // std::cout << "Particle Energy: " << p.energy() << std::endl;
        }
        else // if surface intersect
        {
          // move to the surface
          // do a surface cross. I assume the surface cross will take care of reflections. Otherwise move just past the surface. Also estimator stuff should happen here
          
		    p.cellPointer()->moveParticle(&p, surface_pair.second); // estimator stuff is done here
		    surface_pair.first->crossSurface(&p); // scores estimators, reflects if necessary, pushes an epsilon distance
          double Inew = the_geometry->findCell(&p )->getImportance();
          the_geometry->particleSplitRoulette(&p, &bank, Inew);
          p.recordCell(the_geometry->findCell(&p )); // finds the new cell if a surface interaction happened. Not the most efficient implementation, but whatever
			
		    }
        
        // Kill particles that have left our world
        if (p.cellPointer()->getImportance() == 0.0)
        {
          p.kill();
        }
        
        // Kill particles if they are past the time cut
        if (time_cut)
        {
          if (p.current_time() > max_time)
          {
            p.kill();
          }
        }
          // update the cell the particle is in
          // conduct any splitting, rouletting, or importance = 0 cases
          // update the exp distribution
      }
    }
    
    for (auto e : estimators)
    {
      e->endHistory();
    }
  }
  endSimulation();
}

void simulation::endSimulation()
{
  for (auto e : estimators)
  { 
    e->setNTracks(num_tracks);
  }
}

