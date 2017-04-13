#include <string>
#include <memory>
#include <vector>
#include <iostream>
#include <cassert>
#include <limits>

#include "pugixml.hpp"
#include "Distribution.h"
#include "Estimator.h"
#include "Reaction.h"
#include "Nuclide.h"
#include "Material.h"
#include "Surface.h"
#include "Cell.h"
#include "Source.h"
#include "Geometry.h"
#include "Simulation.h"
#include "InputData.h"


int main() {
	
	inputdata input;
	input.ReadXMLFile();
	std::vector< std::shared_ptr<material> > materials = input.getmaterials();
	std::vector< std::shared_ptr< cell > > cells = input.getcells();
	std::vector< std::shared_ptr< estimator > > estimators = input.getestimators();
	std::shared_ptr< source > src = input.getsrc();
	unsigned long long N_start = input.getN_start();
	unsigned long long number_of_histories = input.getnumber_of_histories();
	double             max_time = input.getmax_time();
	    
	  geometry the_geometry(cells);
	  std::cout << "no error." << std::endl;
	  std::cout << estimators.at(0)->name() << std::endl;
	  simulation the_simulation(&the_geometry, estimators, src);
	  the_simulation.set_max_time(max_time);
	  the_simulation.transport(N_start, number_of_histories);
	  
	  for (auto e : estimators)
	  {
		e->report();
	  }
  
	
	return 0;
}

