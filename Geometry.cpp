#include "Geometry.h"
#include "Particle.h"

geometry::geometry(std::vector<std::shared_ptr<cell>> inCells)
{
  cells = inCells;
}

std::shared_ptr<cell> geometry::findCell(particle* p )
{
  int i = 0;
  int size = cells.size();
  bool foundParticle = false;
  while ((!(foundParticle)) && (i < size))
  {
    foundParticle = cells[i]->testPoint(p->pos());
	if(!foundParticle){i++;}
    
  }
  if (!(foundParticle))
  {
    std::cout << "Lost Particle" << std::endl;
    std::cout << "Position: " << p->pos().x << " " << p->pos().y << " " << p->pos().z << std::endl;
    return(nullptr);
  }
  else
  {
	// if(i == 1){
		// std::cout << i << " Hit Detector!" <<  std::endl;
		// std::cout << cells[i]->name() << std::endl;
		// std::cout << "Position: " << p->pos().x << " " << p->pos().y << " " << p->pos().z << std::endl;
	// }
	
    return(cells[i]);
  }
}

// Inew is the importance of the new cell, or the cell the particle moves into

void geometry::particleSplitRoulette(particle* p, std::stack<particle>* bank, double Inew)
{
  double v = Inew/p->cellPointer()->getImportance();
  if (v>1.0) { //start splitting routine
    double vbar = std::floor(Urand()+v);
    p->adjustWeight(1/vbar);
    for ( int i = 0; i < (vbar-1); i++ ) { bank->push( *p ); }
  }
  else if (v<1.0) {//begin rouletting
    if (Urand()<1-v) {p->kill();}
    else {p->adjustWeight(1/v);}
  }

}
