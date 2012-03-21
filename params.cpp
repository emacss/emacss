
#include "emacss.h"

void get_params(parameters *params){
  /*Contains 'fixed' data for the simulation, ie, the optimised parameters that
    are used by EMACSS to mimic the evolution of a star cluster. These are 
    designed to be user-adjustable, but are set here as a reference, and to 
    represent optimal values recovered in AG2012.*/

  params->xi0 = 0.0142;
  params->rhrj1 = 0.145;
  params->gamma = 0.11;
  params->N1 = 38252;
  params->x = 0.75;
  params->z = 1.61;
  params->frac = 0.1;
}
