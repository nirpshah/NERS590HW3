#ifndef _NUCLIDE_HEADER_
#define _NUCLIDE_HEADER_

#include <vector>
#include <string>
#include <memory>

#include "Reaction.h"

class reaction;
class particle;

class nuclide {
  private:
    std::string nuclide_name;
    double nuclide_A;
    double nuclide_alpha;
    std::vector< std::shared_ptr< reaction > > rxn;
  public:
     nuclide( std::string label, double A_in ) : nuclide_name(label), nuclide_A(A_in) { nuclide_alpha = ((nuclide_A-1)*(nuclide_A-1))/((nuclide_A+1)*(nuclide_A+1));};
    ~nuclide() {};

    std::string name() { return nuclide_name; }
    std::vector< std::shared_ptr< reaction > > getReactions() {return rxn;} ;
    double alpha() { return nuclide_alpha; };
    double A() {return nuclide_A; };
    
    void        addReaction( std::shared_ptr< reaction > );
    double      total_xs(particle* p);
    double      xs(particle* p, std::string reaction_name);
    

    std::shared_ptr< reaction > sample_reaction(particle* p);
};


#endif
