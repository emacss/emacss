/***********************************************************************/
//  Evolve Me A Cluster of StarS (EMACSS)                              //
//                                                                     //
//  By: Poul Alexander (University of Cambridge)                       //
//      Mark Gieles (University of Cambridge, University of Surrey)    //
//      Henny Lamers (University of Amsterdam)                         //
//      Holger Baumgardt (Univesity of Brisbane)                       //
//                                                                     //
//  EMACSS is a numerical integrator that solves for the differential  //
//  equations:                                                         //
//
//  dE/dr = e(N,r,p) = -epsilon/t_rh
//                                                                     //
//  dN/dt = f(N,r,p) = - xi*N/t_rh                                     // 
//                                                                     //
//  dr/dt = g(N,r,p) = (zeta - 2*xi) * r/t_rh                          //
//                                                                     // 
//  where                                                              // 
//                                                                     // 
//  epsilon = f(N,r,rt,t)                                              //
//  xi = f(N,r,rt)                                                     //
//  p is an array of parameters/variables, which define the properties //
//  of star clusters at any given time in the life-cycle               //
//                                                                     //
// Equations referenced AG2012 refer to Alexander & Gieles 2012        //
// Equations referenced AGLB2013 refer to Alexander et al. 2013        //
/***********************************************************************/

#include "../emacss.h"

int main(int argc, char* argv[]){

  node cluster;
  cluster.input(argc,argv);
  
  stellar_evo se_module(&cluster);
  dynamics dynamics_module(&cluster,&se_module);
  
  //Sets adiabatic expansion term and core collapse time for dynamics
  se_module.set_X();
  dynamics_module.set_tcc();
  
  //Writes initial conditions and t=0 output.
  cluster.output(se_module,dynamics_module); 
  
  while (cluster.N > 100 && cluster.time.nbody < cluster.out_time.nbody){
    cluster.evolve(se_module,dynamics_module);	
    cluster.output(se_module,dynamics_module);
  }
  return 0;
}

