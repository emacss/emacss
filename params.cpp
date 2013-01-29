#include "../emacss.h"

  /*Contains 'fixed' data for the simulation, ie, the optimised parameters that
    are used by EMACSS to mimic the evolution of a star cluster. These are 
    designed to be user-adjustable, but are set here as a reference, and to 
    represent optimal values recovered in AG2012 and AGLB2013.*/

node::node(){    //Iterator for a cluster at any given time (initialisation)

  frac = 0.1;
  gamma = 0.02;

}

stellar_evo::stellar_evo(node *innode) : mynode(innode) {

  nu = 0.07;
  T_SE = 2.5;
}

dynamics::dynamics(node *innode, stellar_evo *inse) : mynode(innode),
						       myse(inse){

  R1 = 0.145;
  N1 = 38252;
  x = 0.75;
  z = 1.61;

  if (mynode->s == 1) {
      a = 0.255; F = 0.245; xi0 = 0; T_DYN = 0.2;
  }
  else {
      a = 0; F = 1; xi0 = 0.0142; T_DYN = 0;
  }

}
