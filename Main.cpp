#include <string>
#include <memory>
#include <vector>
#include <iostream>
#include <cassert>
#include <limits>

#include "pugixml.hpp"
#include "Distribution.h"
#include "Reaction.h"
#include "Nuclide.h"
#include "Material.h"
#include "Surface.h"
#include "Cell.h"
#include "Source.h"
#include "Geometry.h"
#include "Simulation.h"

// function that returns an item from a vector of objects of type T by name provided
// the object class has a string and a method called name() allowing for it to be returned
template< typename T >
std::shared_ptr< T > findByName( std::vector< std::shared_ptr< T > > vec, std::string name ) {
  for ( auto v : vec ) {
    if ( v->name() == name ) { return v; }
	
  }

  return nullptr;
}

int main() {
	#include "ReadXMLFile.h"
	
	// the useful vectors are std::vector< std::shared_ptr<material> > materials; std::vector< std::shared_ptr< cell > > cells; std::vector< std::shared_ptr< estimator > > estimators; std::shared_ptr< source > src; 
    
  geometry the_geometry(cells);
  
  simulation the_simulation(&the_geometry, estimators, src);
  the_simulation.transport(N_start, number_of_histories);
  
  for (auto e : estimators)
  {
    e->report();
  }
  
	
	return 0;
}

