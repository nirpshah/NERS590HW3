// user enters the XML file name and pugixml will attempt to load
	std::string input_file_name;
	std::cout << " Enter XML input file name: " << std::endl;
	std::cin  >> input_file_name;

	pugi::xml_document input_file;
	pugi::xml_parse_result load_result = input_file.load_file( input_file_name.c_str() );

	// check to see if result failed and throw an exception if it did
	if ( ! load_result ) {
		std::cout << load_result.description() << std::endl;
		// assert is for code development problems
		// throw is for returning error to user
		throw;
	}

	// distributuions
	std::vector< std::shared_ptr< distribution<double> > >  double_distributions;
	std::vector< std::shared_ptr< distribution<int>    > >  int_distributions;
	std::vector< std::shared_ptr< distribution<point>  > >  point_distributions;
	pugi::xml_node input_distributions = input_file.child("distributions");

	// find total number of distributions
	int num_distributions = 0;
	for ( auto d : input_distributions ) { num_distributions++; }


	// since distributions may depend on other distributions, need to iterate
	int set_distributions = 0;
	while ( set_distributions < num_distributions ) {
		int previous_set_distributions = set_distributions;

		for ( auto d : input_distributions ) {
			std::string type = d.name();
			std::string name = d.attribute("name").value();
			std::string data = d.attribute("datatype").value();

			if ( data == "double" ) {
				// skip rest of loop if distribution already done
				if ( findByName( double_distributions, name ) ) { continue; throw; }

				std::shared_ptr< distribution<double> > Dist;
				if ( type == "delta" ) {
					double a = d.attribute("a").as_double();
					Dist = std::make_shared< arbitraryDelta_distribution< double > > ( name, a );
				}
				else if ( type == "uniform" ) {
					double a = d.attribute("a").as_double();
					double b = d.attribute("b").as_double();
					Dist = std::make_shared< uniform_distribution > ( name, a, b );
				}
				else if ( type == "linear" ) {
					double a  = d.attribute("a").as_double();
					double b  = d.attribute("b").as_double();
					double fa = d.attribute("fa").as_double();
					double fb = d.attribute("fb").as_double();
					Dist = std::make_shared< linear_distribution > ( name, a, b, fa, fb );
				}
				else if ( type == "henyeyGreenstein" ) {
					double a = d.attribute("a").as_double();
					Dist = std::make_shared< HenyeyGreenstein_distribution > ( name, a );
				}
				else {
					std::cout << "unsupported distribution with data type " << data << std::endl;
					throw;
				}
			double_distributions.push_back( Dist );
			}
			// integer-valued distributions
			else if ( data == "int" ) {
				// skip rest of loop if distribution already done
				if ( findByName( int_distributions, name ) ) { continue; throw; }

				std::shared_ptr< distribution<int> > Dist;
				if ( type == "delta" ) {
					double a = d.attribute("a").as_int();
					Dist = std::make_shared< arbitraryDelta_distribution< int > > ( name, a );
				}
				else if ( type == "meanMultiplicity" ) {
					double nubar = d.attribute("nubar").as_double();
					Dist = std::make_shared< meanMultiplicity_distribution > ( name, nubar );
				}
				else if ( type == "terrellFission" ) {
					double nubar = d.attribute("nubar").as_double();
					double sigma = d.attribute("sigma").as_double();
					double b     = d.attribute("b").as_double();
					Dist = std::make_shared< TerrellFission_distribution > ( name, nubar, sigma, b );
				}
				else {
					std::cout << "unsupported distribution with data type " << data << std::endl;
					throw;
				}
				int_distributions.push_back( Dist );
			}
			
			// point value distributions
			else if ( data == "point" ) {
				// skip rest of loop if distribution already done
				if ( findByName( point_distributions, name ) ) { continue; throw; }
				

				std::shared_ptr< distribution< point > > Dist;
				if ( type == "delta" ) {
					double x = d.attribute("x").as_double(); 
					double y = d.attribute("y").as_double(); 
					double z = d.attribute("z").as_double();         
					Dist = std::make_shared< arbitraryDelta_distribution< point > > ( name, point( x, y, z ) );
				}
				else if ( type == "isotropic" ) {
					
					Dist = std::make_shared< isotropicDirection_distribution > ( name );
					
				}
				else if ( type == "anisotropic" ) {
					double u = d.attribute("u").as_double(); 
					double v = d.attribute("v").as_double(); 
					double w = d.attribute("w").as_double();         
					std::shared_ptr< distribution<double> > angDist = 
						findByName( double_distributions, d.attribute("distribution").value() );
		  
					// in the angular distribution does not yet, skip to the end of the loop
					if ( ! angDist ) { continue; throw; }

					Dist = std::make_shared< anisotropicDirection_distribution > ( name, point( u, v, w ), angDist );
				}
				else if ( type == "independentXYZ" ) {
					std::shared_ptr< distribution<double> > distX = findByName( double_distributions, d.attribute("x").value() ); 
					std::shared_ptr< distribution<double> > distY = findByName( double_distributions, d.attribute("y").value() ); 
					std::shared_ptr< distribution<double> > distZ = findByName( double_distributions, d.attribute("z").value() ); 

					// if any of these distributions have not yet been resolved, skip to the end of the loop
					if ( !distX || !distY || !distZ ) { continue; throw; }

					Dist = std::make_shared< independentXYZ_distribution > ( name, distX, distY, distZ );
				}
				else if (type == "discrete"){
					std::shared_ptr< distribution<point> > dist1 = findByName( point_distributions, d.attribute("x").value() ); 
					std::shared_ptr< distribution<point> > dist2 = findByName( point_distributions, d.attribute("y").value() ); 
					std::shared_ptr< distribution<point> > dist3 = findByName( point_distributions, d.attribute("z").value() ); 
					
					// if any of these distributions have not yet been resolved, skip to the end of the loop
					if ( !dist1 || !dist2 || !dist3 ) { continue; throw; }

					double probability1 = d.attribute("prob1").as_double();
					double probability2 = d.attribute("prob2").as_double();
					double probability3 = d.attribute("prob3").as_double();
					
					
					
					std::vector< std::pair<point, double >> discrete_vector;
					discrete_vector.push_back(std::make_pair(dist1->sample(),probability1));
					discrete_vector.push_back(std::make_pair(dist2->sample(),probability2));
					discrete_vector.push_back(std::make_pair(dist3->sample(),probability3));
					
					Dist = std::make_shared< arbitraryDiscrete_distribution<point> > (name, discrete_vector ); 
				}
				else if (type == "shell"){
					double r1 = d.attribute("r1").as_double(); 
					double r2 = d.attribute("r2").as_double();       
					Dist = std::make_shared< shell_distribution> ( name, r1, r2 );
				}
				else if (type == "zdisk"){
					double x_center  = d.attribute("x").as_double(); 
					double y_center  = d.attribute("y").as_double();   
					double z_center  = d.attribute("z").as_double(); 
					point p_center(x_center,y_center,z_center);
					double disk_rad	 = d.attribute("rad").as_double();
					Dist = std::make_shared< zdisk_distribution> ( name, p_center, disk_rad );
				}
				else if (type == "xdisk"){
					double x_center  = d.attribute("x").as_double(); 
					double y_center  = d.attribute("y").as_double();   
					double z_center  = d.attribute("z").as_double(); 
					point p_center(x_center,y_center,z_center);
					double disk_rad	 = d.attribute("rad").as_double();
					Dist = std::make_shared< xdisk_distribution> ( name, p_center, disk_rad );
				}
				else if (type == "forwardPeak"){
					Dist = std::make_shared< forwardpeak_distribution> ( name);
				}
				else {
					std::cout << "unsupported " << data << " distribution of type " << type << std::endl;
					throw;
				}
				point_distributions.push_back( Dist );

			}
			else {
				std::cout << "unsupported distribution with data type " << data << std::endl;
				throw;
			}
			// if we reach here, assume distribution has been set
			set_distributions++;
			
		}
		// check to see if number of distributions has increased, if not, caught in an infinite loop
		if ( previous_set_distributions == set_distributions ) { 
			std::cout << "distributions could not be resolved. " << std::endl;
			throw;
		}
	}

		

	// iterate over nuclides
	std::vector< std::shared_ptr<nuclide> > nuclides;
	pugi::xml_node input_nuclides = input_file.child("nuclides");
	for ( auto n : input_nuclides ) {
		std::string name = n.attribute("name").value();
		double A = n.attribute("A").as_double();

		std::shared_ptr< nuclide > Nuc = std::make_shared< nuclide > ( n.attribute("name").value(), A );
		nuclides.push_back( Nuc );

		// iterate over its reactions
		for ( auto r : n.children() ) {
			std::shared_ptr< reaction > Rxn;
			std::string rxn_type = r.name();

			if ( rxn_type == "capture" ) {
			  std::string xs_type = r.attribute("xstype").value();
			  double xs = 0;
			  double a = 0;
			  double b = 0;
			  if (xs_type == "constant")
			  {
			    xs = r.attribute("xs").as_double();
			    a = 0;
			    b = 0;
			  }
			  else if (xs_type == "equation")
			  {
			    xs = 0;
			    a = r.attribute("a").as_double();
			    b = r.attribute("b").as_double();
			  }
			  else
			  {
			    std::cout << "unknown xs type " << xs_type << " in nuclide " << name << std::endl;
			    throw;
			  }
				Nuc->addReaction( std::make_shared< capture_reaction > ( xs_type, xs, a, b ) );
			}
			else if ( rxn_type == "scatter" ) {
				std::string dist_name = r.attribute("distribution").value();
				std::shared_ptr< distribution<double> > scatterDist = findByName( double_distributions, dist_name );
				std::string xs_type = r.attribute("xstype").value();
				double xs = 0;
				double a = 0;
				double b = 0;
			  if (xs_type == "constant")
			  {
			    xs = r.attribute("xs").as_double();
			    a = 0;
			    b = 0;
			  }
			  else if (xs_type == "equation")
			  {
			    xs = 0;
			    a = r.attribute("a").as_double();
			    b = r.attribute("b").as_double();
			  }
			  else
			  {
			    std::cout << "unknown xs type " << xs_type << " in nuclide " << name << std::endl;
			    throw;
			  }
				if ( scatterDist ) {
					Nuc->addReaction( std::make_shared< scatter_reaction > ( xs_type, xs, a, b, scatterDist ) );
				}
				else {
					std::cout << " unknown scattering distribution " << dist_name << " in nuclide " << name << std::endl;
					throw;
				}
			}
			else if ( rxn_type == "fission" ) {
				std::string mult_dist_name = r.attribute("multiplicity").value();
				std::string ener_dist_name = r.attribute("energy").value();
				std::shared_ptr< distribution<int> > multDist = findByName( int_distributions, mult_dist_name );
				std::shared_ptr< distribution<double> > enerDist = findByName( double_distributions, ener_dist_name);
				std::string xs_type = r.attribute("xstype").value();
				double xs = 0;
				double a = 0;
				double b = 0;
			  if (xs_type == "constant")
			  {
			    xs = r.attribute("xs").as_double();
			    a = 0;
			    b = 0;
			  }
			  else if (xs_type == "equation")
			  {
			    xs = 0;
			    a = r.attribute("a").as_double();
			    b = r.attribute("b").as_double();
			  }
			  else
			  {
			    std::cout << "unknown xs type " << xs_type << " in nuclide " << name << std::endl;
			    throw;
			  }
				if ( multDist && enerDist) {
					  Nuc->addReaction( std::make_shared< fission_reaction > ( xs_type, xs, a, b, multDist, enerDist ) );
				}				
				else {
				  if ( ! multDist) { std::cout << " unknown multiplicity distribution " << mult_dist_name << " in nuclide " << name << std::endl; }
				  if ( ! enerDist) { std::cout << " unknown fission energy distribution " << ener_dist_name << " in nuclide " << name << std::endl; }
				  throw;
			  }
			}
			else {
				std::cout << "unknown reaction type " << rxn_type << std::endl;
				throw;
			}	
		}
	} 

	// iterate over materials
	std::vector< std::shared_ptr<material> > materials;
	pugi::xml_node input_materials = input_file.child("materials");
	for ( auto m : input_materials ) {
		std::string name = m.attribute("name").value();
		double      aden = m.attribute("density").as_double();
    
		std::shared_ptr< material > Mat = std::make_shared< material > ( name, aden );    
		materials.push_back( Mat );

		// iterate over nuclides
		for ( auto n : m.children() ) {
			if ( (std::string) n.name() == "nuclide" ) {
				std::string nuclide_name = n.attribute("name").value();
				double      frac         = n.attribute("frac").as_double();
        
				Mat->addNuclide( findByName( nuclides, nuclide_name ), frac );
			}
		}
	}

	// iterate over surfaces
	std::vector< std::shared_ptr< surface > > surfaces;
	pugi::xml_node input_surfaces = input_file.child("surfaces");
	for ( auto s : input_surfaces ) {
		std::string type = s.name();

		std::shared_ptr< surface > S;
		if ( type == "plane") {
			std::string name = s.attribute("name").value();
			double      a    = s.attribute("a").as_double();
			double      b    = s.attribute("b").as_double();
			double      c    = s.attribute("c").as_double();
			double      d    = s.attribute("d").as_double();
			S = std::make_shared< plane > ( name, a, b, c, d );
		} else if (type == "sphere" ) {
			std::string name = s.attribute("name").value();
			double      a    = s.attribute("a").as_double();
			double      b    = s.attribute("b").as_double();
			double      c    = s.attribute("c").as_double();
			double      d    = s.attribute("d").as_double();
			S = std::make_shared< sphere > ( name, a, b, c, d );
		} else if (type == "xcylinder"){
			std::string name = s.attribute("name").value();
			double      a    = s.attribute("a").as_double();
			double      b    = s.attribute("b").as_double();
			double      c    = s.attribute("c").as_double();
			S = std::make_shared< xcylinder > ( name, a, b, c);
		} else if (type == "zcylinder"){
			std::string name = s.attribute("name").value();
			double      a    = s.attribute("a").as_double();
			double      b    = s.attribute("b").as_double();
			double      c    = s.attribute("c").as_double();
			S = std::make_shared< zcylinder > ( name, a, b, c);
		} else if (type == "xcone"){
			std::string name = s.attribute("name").value();
			double      x_0  = s.attribute("x0").as_double();
			double      y_0  = s.attribute("y0").as_double();
			double      z_0  = s.attribute("z0").as_double();
			double      R_in = s.attribute("R").as_double();
			S = std::make_shared< xcone > ( name, x_0, y_0, z_0, R_in);
		}
		else {
			std::cout << " unkown surface type " << type << std::endl;
			throw;
		}

		if ( (std::string) s.attribute("bc").value() == "reflect" ) {
			S->makeReflecting();
		}
		surfaces.push_back( S );
	}

	// iterate over cells
	std::vector< std::shared_ptr< cell > > cells;
	pugi::xml_node input_cells = input_file.child("cells");
	for ( auto c : input_cells ) {
		std::string name = c.attribute("name").value();

		std::shared_ptr< cell > Cel = std::make_shared< cell > ( name );
		

		// cell material
		if ( c.attribute("material") ) {
			std::shared_ptr< material > matPtr = findByName( materials, c.attribute("material").value() );
			if ( matPtr ) {
				Cel->setMaterial( matPtr );
			}
			else {
				std::cout << " unknown material " << c.attribute("material").value() << " in cell " << name << std::endl;
				throw;
			} 
		}

		// cell importance
		if ( c.attribute("importance") ) {
			Cel->setImportance( c.attribute("importance").as_double() );
		}
   
		// iterate over surfaces
		for ( auto s : c.children() ) {
			if ( (std::string) s.name() == "surface" ) {
				std::string name  = s.attribute("name").value();
				int         sense = s.attribute("sense").as_int();

				std::shared_ptr< surface > SurfPtr = findByName( surfaces, name );
				if ( SurfPtr ) {
					Cel->addSurface( findByName( surfaces, name ), sense );
				}
				else {
					std::cout << " unknown surface with name " << name << std::endl;
					throw;
				}
			}
			else {
				std::cout << " unknown data type " << s.name() << " in cell " << name << std::endl;
				throw;
			}
		} 
		cells.push_back( Cel );
	}

	// iterate over estimatators
	std::vector< std::shared_ptr< estimator > > estimators;
	pugi::xml_node input_estimators = input_file.child("estimators");
	for ( auto e : input_estimators ) {
		std::string type = e.name();
		std::string name = e.attribute("name").value();
    
		std::shared_ptr< estimator > Est;
		if ( type == "current" ) {
			Est = std::make_shared< surface_current_estimator > ( name );

			// get the surfaces
			for ( auto s : e.children() ) {
			if ( (std::string) s.name() == "surface" ) {
				std::string name = s.attribute("name").value();
				std::shared_ptr< surface > SurfPtr = findByName( surfaces, name );
				if ( SurfPtr ) {
					SurfPtr->attachEstimator( Est );
				}
				else {
					std::cout << " unknown surface label " << name << " in estimator " << e.attribute("name").value() << std::endl;
				}		
			}
			} 
		}
		else if ( type == "countingSurface" ) {
			Est = std::make_shared< counting_estimator > ( name );

			// get the surfaces
			for ( auto s : e.children() ) {
				if ( (std::string) s.name() == "surface" ) {
					std::string name = s.attribute("name").value();
					std::shared_ptr< surface > SurfPtr = findByName( surfaces, name );
					if ( SurfPtr ) {
						SurfPtr->attachEstimator( Est );
					}
					else {
						std::cout << " unknown surface label " << name << " in estimator " << e.attribute("name").value() << std::endl;
					}
				}
			} 
		}
    else if ( type == "pathLengthReactionRate" ) 
    {
		  double volume = e.attribute("volume").as_double();
		  std::string reaction_name = e.attribute("reactionname").value();
     	Est = std::make_shared< cell_pathLengthReactionRate_estimator > ( name, volume, reaction_name );

			// get the cell
			for ( auto s: e.children()) {
				if ( ( std::string) s.name() == "cell" ) {
					std::string name = s.attribute("name").value();
					std::shared_ptr< cell > CellPtr = findByName( cells, name );
					if ( CellPtr ) {
						CellPtr->attachEstimator( Est );
					}
					else {
						std::cout << " unknown cell label " << name << " in estimator " << e.attribute("name").value() << std::endl;
					}
				}
			}
		}
		else if ( type == "pathLengthFlux" ) {
			double volume = e.attribute("volume").as_double();
			Est = std::make_shared< cell_pathLengthFlux_estimator > ( name, volume );

			// get the cell
			for ( auto s: e.children()) {
				if ( ( std::string) s.name() == "cell" ) {
					std::string name = s.attribute("name").value();
					std::shared_ptr< cell > CellPtr = findByName( cells, name );
					if ( CellPtr ) {
						CellPtr->attachEstimator( Est );
					}
					else {
						std::cout << " unknown cell label " << name << " in estimator " << e.attribute("name").value() << std::endl;
					}
				}
			}
		}
		else if ( type == "pathLengthTimeBin" ) {
			std::string reaction_name = e.attribute("reactionname").value();
			double binnum = e.attribute("binnum").as_double();
			double binmin = e.attribute("binmin").as_double();
			double binmax = e.attribute("binmax").as_double();
			Est = std::make_shared< cell_pathLengthTimeBin_estimator > ( name, reaction_name, binnum, binmin, binmax );

			// get the cell
			for ( auto s: e.children() ) {
				if ( ( std::string ) s.name() == "cell" ) {
					std::string name = s.attribute("name").value();
					std::shared_ptr< cell > CellPtr = findByName ( cells, name );
					if ( CellPtr ) {
						CellPtr->attachEstimator ( Est );
					}
					else {
						std::cout << " uhnknown cell label " << name << " in estimator " << e.attribute("name").value() << std::endl;
					}
				}
			}
		}
		else {
			std::cout << "unknown estimator type " << name << std::endl;
			throw;
		}
    
		estimators.push_back( Est );
	}

	// create source
	pugi::xml_node input_source = input_file.child("source");
	pugi::xml_node input_source_position  = input_source.child("position");
	pugi::xml_node input_source_direction = input_source.child("direction");
	pugi::xml_node input_source_energy    = input_source.child("energy");
	pugi::xml_node input_source_time      = input_source.child("time");

	std::string pos_dist_name   = input_source_position.attribute("distribution").value();
	std::string dir_dist_name   = input_source_direction.attribute("distribution").value();
	std::string ener_dist_name  = input_source_energy.attribute("distribution").value();
	std::string time_dist_name  = input_source_time.attribute("distribution").value();

	std::shared_ptr< distribution< point > > posDist = findByName( point_distributions, pos_dist_name );
	std::shared_ptr< distribution< point > > dirDist = findByName( point_distributions, dir_dist_name );
	std::shared_ptr< distribution< double> > enerDist= findByName( double_distributions, ener_dist_name);
	std::shared_ptr< distribution< double> > timeDist= findByName( double_distributions, time_dist_name);

	std::shared_ptr< source > src;  
	if ( posDist && dirDist && enerDist && timeDist) {
		src = std::make_shared< source > ( posDist, dirDist, enerDist, timeDist);  
	}	
	else {
		if ( ! posDist ) { std::cout << " unknown position distribution "  << pos_dist_name << " in source " << std::endl; }
		if ( ! dirDist ) { std::cout << " unknown direction distribution " << dir_dist_name << " in source " << std::endl; }
		if ( ! enerDist) { std::cout << " unknown energy distribution " << ener_dist_name << " in source " << std::endl; }
		if ( ! timeDist) { std::cout << " unknown time distribution " << time_dist_name << " in source " << std::endl; }
		throw;
	
	}
	
	// create simulation
	pugi::xml_node input_simulation 		= input_file.child("simulation");
	pugi::xml_node input_histories  		= input_simulation.child("histories");
	unsigned long long N_start   			= input_histories.attribute("start").as_int();
	unsigned long long number_of_histories 	= input_histories.attribute("end").as_int();
	double             max_time             = input_simulation.child("timecut").attribute("maxtime").as_double();
	
