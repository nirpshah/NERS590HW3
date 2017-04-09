#ifndef _REACTION_HEADER_
#define _REACTION_HEADER_

#include <vector>
#include <iostream>
#include <string>
#include <memory>
#include <stack>
#include <utility>

#include "Particle.h"
#include "Distribution.h"
#include "Nuclide.h"

class particle;
class nuclide;

class reaction {
  private:
    std::string xs_type;
    double rxn_xs;
    double a;
    double b;
    
  protected:
    std::string rxn_name;
  public:
     reaction( std::string xs_type_in, double rxn_xs_in, double a_in, double b_in ) : xs_type(xs_type_in), rxn_xs(rxn_xs_in), a(a_in), b(b_in) {};
    ~reaction() {};

    virtual std::string name() final { return rxn_name; };
    virtual double xs(particle* p) final;
    virtual void sample( particle* p, std::stack<particle>* bank, std::shared_ptr< nuclide > N ) = 0;
};

class capture_reaction : public reaction {
  private:
 
  public:
     capture_reaction( std::string xs_type_in, double xs_in, double a_in, double b_in ) : reaction(xs_type_in, xs_in, a_in, b_in) { rxn_name = "capture"; };
    ~capture_reaction() {};

    void sample( particle* p, std::stack<particle>* bank, std::shared_ptr< nuclide > N );
};

class scatter_reaction : public reaction {
  private:
    std::shared_ptr< distribution<double> > scatter_dist; 
  public:
     scatter_reaction( std::string xs_type_in, double xs_in, double a_in, double b_in, std::shared_ptr< distribution<double> > D ) :
       reaction(xs_type_in, xs_in, a_in, b_in), scatter_dist(D) { rxn_name = "scatter"; };
    ~scatter_reaction() {};

    void sample( particle* p, std::stack<particle>* bank, std::shared_ptr< nuclide > N );
};

class fission_reaction : public reaction {
  private:
    std::shared_ptr< distribution<int> >    multiplicity_dist; 
    std::shared_ptr< distribution<point> >  isotropic;
    std::shared_ptr< distribution<double> > energy_dist;
    
  public:
     fission_reaction(  std::string xs_type_in, double xs_in, double a_in, double b_in, 
                        std::shared_ptr< distribution<int> > D1, 
                        std::shared_ptr< distribution<double> > D2 ) 
                      : reaction(xs_type_in, xs_in, a_in, b_in),
                        multiplicity_dist(D1), 
                        energy_dist(D2)  
                      { rxn_name = "fission";
                        isotropic = std::make_shared< isotropicDirection_distribution > ( "isotropic" ); 
                      };
    ~fission_reaction() {};

    void sample( particle* p, std::stack<particle>* bank, std::shared_ptr< nuclide > N );
};

#endif
