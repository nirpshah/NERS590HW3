#include <cmath>
#include <iostream>

#include "Estimator.h"
#include "Cell.h"
#include "Material.h"
#include "Nuclide.h"
#include "Reaction.h"
#include "Particle.h"

void surface_current_estimator::score( particle* p, double null ) { tally_hist += p->wgt();}

void counting_estimator::score( particle* p, double null ) { count_hist++; }

void counting_estimator::endHistory() {
  if ( tally.size() < count_hist + 1 ) { tally.resize( count_hist + 1, 0.0 ); }
  tally[ count_hist ] += 1.0;
  nhist++;
  count_hist = 0;
}

double  counting_estimator::relError(){
	return 0;
}

void counting_estimator::report() {
  std::cout << name() << std::endl;
  double s1 = 0.0, s2 = 0.0;
  for ( int i = 0 ; i < tally.size() ; i++ ) 
  {
    double p = tally[i] / nhist;
    std::cout << i << " " << p << "   " << std::sqrt( p * ( 1.0 - p ) / nhist ) / p <<  std::endl;
    s1 += p * i;
    s2 += p * i * i;
  }
  std::cout << "   mean = " << s1 << std::endl;
  std::cout << "   var  = " << s2 - s1*s1 << std::endl;
}

void cell_pathLengthFlux_estimator::score( particle* p, double path_length ) { tally_hist += p->wgt() * path_length; }

void cell_pathLengthReactionRate_estimator::score( particle* p, double path_length ) {
  tally_hist += p->wgt() * path_length * p->cellPointer()->getMaterial()->macro_xs(p, reaction_name);
} 

  
  
  
  
  
  
