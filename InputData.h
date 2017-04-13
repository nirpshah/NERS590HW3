#ifndef _INPUTDATA_HEADER_
#define _INPUTDATA_HEADER_

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <typeinfo>
#include <vector>
#include <memory>
#include <string>
#include <iostream>
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



class inputdata
{
protected:
    std::vector< std::shared_ptr< distribution<double> > >  double_distributions;
	std::vector< std::shared_ptr< distribution<int>    > >  int_distributions;
	std::vector< std::shared_ptr< distribution<point>  > >  point_distributions;
	std::vector< std::shared_ptr<nuclide> > nuclides;
	std::vector< std::shared_ptr<material> > materials;
	std::vector< std::shared_ptr< surface > > surfaces;
	std::vector< std::shared_ptr< cell > > cells;
	std::vector< std::shared_ptr< estimator > > estimators;
	std::shared_ptr< source > src;  
	unsigned long long N_start;
	unsigned long long number_of_histories;
	double             max_time;
	

public:
	inputdata() {};
	~inputdata() {};
    void ReadXMLFile();
	std::vector< std::shared_ptr< distribution<double> > >  getdouble_distributions() {return double_distributions;};
	std::vector< std::shared_ptr< distribution<int>    > >  getint_distributions(){return int_distributions;};
	std::vector< std::shared_ptr< distribution<point>  > >  getpoint_distributions() {return point_distributions;};
	std::vector< std::shared_ptr<nuclide> > getnuclides() {return nuclides;};
	std::vector< std::shared_ptr<material> > getmaterials() {return materials;};
	std::vector< std::shared_ptr< surface > > getsurfaces() {return surfaces;};
	std::vector< std::shared_ptr< cell > > getcells() { return cells;};
	std::vector< std::shared_ptr< estimator > > getestimators() {estimators;};
	std::shared_ptr< source > getsrc() {return src;};  
	unsigned long long getN_start(){return N_start;};
	unsigned long long getnumber_of_histories(){return number_of_histories;};
	double             getmax_time(){return max_time;};
	
	
	
};



#endif