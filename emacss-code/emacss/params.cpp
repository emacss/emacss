#include "../emacss.h"

  /*Contains 'fixed' data for the simulation, ie, the optimised parameters that
    are used by EMACSS to mimic the evolution of a star cluster. These are 
    designed to be user-adjustable, but are set here as a reference, and to 
    represent optimal values recovered in AG2012 and AGLB2013.*/

node::node(){    //Iterator for a cluster at any given time (initialisation)

  frac = 0.01;
  gamma = 0.02;
  T_SE = 2; 

}

stellar_evo::stellar_evo(){}

void stellar_evo::load(node *innode, dynamics *indyn){
    
  mynode = innode;
  mydyn = indyn;

  nu = 0.072;
  MS_1 = 4.0;
  y = -0.3;
  z0 = 2.0;
}

dynamics::dynamics(){}

void dynamics::load(node *innode, stellar_evo *inse){
    
  mynode = innode;
  myse = inse;

  R1 = 0.22;
  N1 = 1000;
  x = 0.75;
  z = 2.0;
  f = 0.3;
  k_1 = 0.25;

  
  if (mynode->s == 1) {F = 0.2; xi0 = 0.0142; mynode->T_DYN = 1.7;}
  else {F = 0; xi0 = 0.0142; mynode->T_DYN = 25.0;}
}
