#include <cmath>
#include <iostream>

#include "Estimator.h"
#include "Cell.h"
#include "Material.h"
#include "Nuclide.h"
#include "Reaction.h"
#include "Particle.h"
#include "Utility.h"

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

  
void cell_pathLengthTimeBin_estimator::score( particle* p, double path_length ) {

	
	int binstart, binend;
	binstart = bin_search( binpoints, p->past_time() );
	binend   = bin_search( binpoints, p->current_time() );
	double delt = p->current_time() - p->past_time();
	double sigS = p->cellPointer()->getMaterial()->macro_xs(p, reaction_name);
	double tcf = 1.0; // use tcf (time correction factor) to truncate the score for when 
						// the time extends beyond 100 ns
	
	
	if ( binstart - binend == 0 ) { 
		if (p->current_time() > binmax){
			tcf = (binmax - p->past_time())/(p->current_time() - p->past_time());
		}
	tally_hist[binend] += p->wgt() * path_length * sigS * tcf; 
	} //only one bin

	else if ( binstart - binend == 1 ) { // part of one bin and part of the next bin only
		if (p->current_time() > binmax){
			tcf = (binmax - binpoints[binnum - 2])/(p->current_time() - binpoints[binnum - 2]);
		}
		tally_hist[binstart] += ( binpoints[binstart] - p->past_time() ) * p->wgt() * path_length * sigS / delt;
		tally_hist[binend] += ( p->current_time() - binpoints[binstart] ) * p->wgt() * path_length * sigS * tcf / delt;
	}

	else { // multiple bins; loop through interior bins as normal
		if (p->current_time() > binmax){
			tcf = (binmax - binpoints[binnum - 2])/(p->current_time() - binpoints[binnum - 2]);
		}
	
		tally_hist[binstart] += ( binpoints[binstart] - p->past_time() ) * p->wgt() * path_length * sigS / delt;
		tally_hist[binend] += (p->current_time() - binpoints[binend-1] ) * p->wgt() * path_length * sigS * tcf / delt;
		for ( int i = binstart+1; i < binend; i++ ) {
			tally_hist[i] += binmesh * p->wgt() * path_length * sigS / delt;
		}		
	}


}

void cell_pathLengthTimeBin_estimator::endHistory() {
  nhist++;
  for ( int i=0; i<binnum; i++ ) {
    tally_sum[i] += tally_hist[i];
    tally_squared[i] += tally_hist[i] * tally_hist[i];
    tally_hist[i] = 0.0;
  }
}

void cell_pathLengthTimeBin_estimator::report() {
  for ( int i=0; i < binnum; i++ ) {
    mean[i] = tally_sum[i]/nhist;
    var[i] = ( tally_squared[i] / nhist - mean[i] * mean[i] ) / nhist;
    std::cout << name() << " bin upper bound: " << binpoints[i] << " value " << mean[i] << " " << std::sqrt( var[i] ) / mean[i] << std::endl;
  }
}

double cell_pathLengthTimeBin_estimator::relError() { // not implemented
  return 0.0;
}  
  
  
  
  
