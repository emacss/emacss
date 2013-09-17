#include "../emacss.h"

  /*Contains 'fixed' data for the simulation, ie, the optimised parameters that
    are used by EMACSS to mimic the evolution of a star cluster. These are 
    designed to be user-adjustable, but are set here as a reference, and to 
    represent optimal values recovered in AG2012 and AGLB2013.*/

node::node(){    //Iterator for a cluster at any given time (initialisation)

  gamma = 0.11;
  nvar = 11;                           //Number of parameters
  k1 = 0.295;
  Rch0 = 0.1;

 }

dynamics::dynamics(node *innode) : mynode(innode){
  //From AG2012
  R1 = 0.145;
  N1 = 1.5e4;                          //Changed from 38252 (AG2012)in GALB2013
  x = 0.75;
  z = 1.61;
  xi1 = 0.0142;  
   
  N2 = 12.0;
  N3 = N1;
  f = 0.3;
  delta_1 = -0.9*mynode->E.zeta;
  delta_2 = -0.02*mynode->E.zeta;
  }  
