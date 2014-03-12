#include "../emacss.h"

  /*Contains 'fixed' data for the simulation, ie, the optimised parameters that
    are used by EMACSS to mimic the evolution of a star cluster. These are 
    designed to be user-adjustable, but are set here as a reference, and to 
    represent optimal values recovered in AG2012, GALB2013 and AGLB2014.*/

node::node(){    //Iterator for a cluster at any given time (initialisation)

  //From AG2012 (modified for AGLB2014)
  nvar = 13;                           //Number of parameters
  gamma = 0.02;                        //Modified in input if required
  
  //From GALB2013
  Rch0 = 0.1;
  frac = 1e-4;
 }

stellar_evo::stellar_evo(){}

dynamics::dynamics(){}

void stellar_evo::load(node *innode, dynamics *indyn){
    
  mynode = innode;
  mydyn = indyn;
  
  m_ref = 100;          //Sets stellar Evo properties, in Myr & Msun
  t_ref = 3.3;
  m_ns = 1.2;
  t_se = 3.3;
  
  nu = 0.07;
  MS_1 = 4.0;
  y = -0.3;
   
  //From AGLB2014
  psi1 = 8;
  psi0 = 1.6;
}

void node::load(dynamics *indyn, stellar_evo *inse){

  //Initialises stellar evolution
  inse->tse(M_star,T_star);
  //Initialises dynamical friction
  if (galaxy.f == 1) indyn->tdf(M_star, R_star, pcMyr, T_star);
}

void dynamics::load(node *innode, stellar_evo *inse){
    
  mynode = innode;
  myse = inse;  
    
  //General
  x = 0.75;
  f = 0.3;
   
  if (mynode->s == 0){         //From AG2012  (modified for GALB2013)      
  R1 = 0.145; N1 = 1.5e4; z = 1.61; xi1 = 0.0142; Fej = 0; nc = 0;
  mynode->k1 = 0.295;      //Set in node for easier calls
  }
  else{                        //From AGLB2014
  R1 = 0.22; N1 = 1000; z = 2.0; xi1 = 0.0075; Fej = 0.55; nc = 12.5;  
  mynode->k1 = 0.24;       //Set in node for easier calls
  }
  
  //From GALB2013
  N2 = 12.0;
  N3 = N1;
  delta_1 = -0.09;
  delta_2 = -0.002;
  
  //From AGLB2013
  Y = 100;
  b = 1.35;
  q = 2;
  }  
