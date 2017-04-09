#ifndef _ESTIMATOR_HEADER_
#define _ESTIMATOR_HEADER_

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <typeinfo>

#include "Particle.h"
#include "Material.h"
#include "Reaction.h"

// base estimator class
class estimator {
  private:
    std::string estimator_name;
  protected:
    unsigned long long nhist;
  public:
     estimator( std::string label ) : estimator_name(label) {};
    ~estimator() {};

    virtual std::string name() final { return estimator_name; };

    virtual void score( particle* , double  ) = 0; // score at events

    template< typename T >
    void score( particle*, T ) { assert(false); };

    void score( particle*, double, std::shared_ptr< material > ) {};

    virtual void endHistory()       = 0; // closeout history
    virtual void report()           = 0; // write results
    virtual double relError()       = 0; // gets relative error for figure of merit
};

// derived class for simple estimators like current or scalar flux
// future estimators could get spectra or flux on a mesh
class single_valued_estimator : public estimator {
  private:

  protected:
    double tally_hist, tally_sum, tally_squared, var, mean;
  public:
     using estimator::score;

     single_valued_estimator(std::string label ) : estimator(label) { 
       nhist         = 0;
       tally_hist    = 0.0;   
       tally_sum     = 0.0; 
       tally_squared = 0.0;
     };
    ~single_valued_estimator() {};

     virtual void endHistory()    final { 
       nhist++;
       tally_sum     += tally_hist;
       tally_squared += tally_hist * tally_hist;
       tally_hist = 0.0; }

     virtual void score( particle*, double ) = 0;

     virtual void report() final {
       mean = tally_sum / nhist;
       var  = ( tally_squared / nhist - mean*mean ) / nhist;
       std::cout << name() << "   " << mean << "   " << std::sqrt( var ) / mean << std::endl;  
     };
     virtual double relError() final{
       mean = tally_sum / nhist;
       var = ( tally_squared / nhist - mean*mean ) / nhist;
       return std::sqrt( var ) / mean;
     };
};

// current crossing a surface
class surface_current_estimator : public single_valued_estimator {
  private:

  public:
     surface_current_estimator( std::string label ) : single_valued_estimator(label) {};
    ~surface_current_estimator() {};

    void score( particle*, double );
};

class counting_estimator : public estimator {
  private:
    int count_hist;
    std::vector< double > tally;
  public:
     counting_estimator( std::string label ) : estimator(label) { count_hist = 0; };
    ~counting_estimator() {};

    void score( particle*, double  );
    void endHistory();
    void report();
double relError();
};

// volume averaged scalar flux in a cell
class cell_pathLengthFlux_estimator : public single_valued_estimator {
  private:
    double volume;
  public:
    cell_pathLengthFlux_estimator( std::string label, double vol ) : single_valued_estimator(label), volume(vol) {};
   ~cell_pathLengthFlux_estimator() {};
    
    void score( particle*, double );
};

// volume-averaged reaction rate in a cell
// NOTE: CURRENT IMPLEMENTATION IS FOR ONLY THE FIRST NUCLIDE IN A MATERIAL, AND THE FIRST REACTION RATE IN THE NUCLIDE. 
// NEEDS TO BE CHANGED TO BE MORE FLEXIBLE
class cell_pathLengthReactionRate_estimator : public single_valued_estimator {
  private:
    double volume;
    std::string reaction_name;
  public: 
    cell_pathLengthReactionRate_estimator( std::string label, double vol, std::string reaction_name_in ) :  single_valued_estimator(label), volume(vol), reaction_name(reaction_name_in) {};
   ~cell_pathLengthReactionRate_estimator() {};
    
    void score( particle*, double );
};

#endif
