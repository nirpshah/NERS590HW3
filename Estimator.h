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
#include "Utility.h"

// base estimator class
class estimator {
  private:
    std::string estimator_name;
  protected:
    unsigned long long nhist;
    unsigned long long ntracks;
  public:
     estimator( std::string label ) : estimator_name(label) { nhist = 0; ntracks = 0;};
    ~estimator() {};

    virtual std::string name() final { return estimator_name; };
    virtual void setNTracks(unsigned long long num_tracks) final 
    { 
      ntracks = num_tracks; 
    } 

    virtual void score( particle* , double  ) = 0; // score at events

    template< typename T >
    void score( particle*, T ) { assert(false); };

    void score( particle*, double, std::shared_ptr< material > ) {};

    virtual void    endHistory()       = 0; // closeout history
    virtual void    report()           = 0; // write results
    virtual double  relError()         = 0; // gets relative error for figure of merit
    
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
       double FOM  = 1.0/(relError()*relError()*ntracks);
       std::cout << name() << "   " << mean << "   " << std::sqrt( var ) / mean << std::endl;  
       std::cout << "FOM: " << FOM << std::endl;
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
class cell_pathLengthReactionRate_estimator : public single_valued_estimator {
  private:
    double volume;
    std::string reaction_name;
  public: 
    cell_pathLengthReactionRate_estimator( std::string label, double vol, std::string reaction_name_in ) :  single_valued_estimator(label), volume(vol), reaction_name(reaction_name_in) {};
   ~cell_pathLengthReactionRate_estimator() {};
    
    void score( particle*, double );
};

// time binned detector response
class cell_pathLengthTimeBin_estimator : public estimator {
  private:
    std::string reaction_name;
    double binnum, binmin, binmax;
  protected:
    double binmesh;
    std::vector<double> binpoints, tally_hist, tally_sum, tally_squared, mean, var;
  public:
    cell_pathLengthTimeBin_estimator( std::string label, std::string reaction_name_in, double nbin, double minbin, double maxbin ) : estimator(label), reaction_name(reaction_name_in), binnum(nbin), binmin(minbin), binmax(maxbin) {
      binmesh = ( binmax - binmin ) / binnum;
      binpoints.resize(binnum);
      tally_hist.resize(binnum);
      tally_sum.resize(binnum);
      tally_squared.resize(binnum);
      mean.resize(binnum);
      var.resize(binnum);
    for ( int i=0; i <= binnum-1; i++ ) {
      binpoints[i] = ( i + 1 ) * binmesh; // creates an array of the uppoer bounds of each bin for sorting
      tally_hist[i]    = 0.0;
      tally_sum[i]     = 0.0;
      tally_squared[i] = 0.0;
      mean[i]          = 0.0;
      var[i]           = 0.0;
    }
  };
 ~cell_pathLengthTimeBin_estimator() {};
  
  void score( particle*, double );
  void endHistory();
  void report();
  double relError(); // not implemented
};

#endif
